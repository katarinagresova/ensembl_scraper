import pandas as pd
from tqdm import tqdm
import logging
from pathlib import Path
from utils import download_file, get_2bit_genome_file, prepare_temp_directory, prepare_data_directory, config, save_to_fasta, save_test_to_fasta
from random_negatives import generate_negatives
from preprocessing import remove_low_quality, split_train_val_test
import pyfiglet
import yaml

logging.basicConfig(level=logging.INFO)


def get_supported_organisms() -> list:
    """Get list of all supported organisms specified in a config file.

    Returns
    -------
    list
        a list of supported organisms
    """
    return list(config['organisms'].keys())


def get_fasta_path(organism: str) -> str:
    """Get Ensembl URL for fasta file of given organism

    Parameters
    ----------
    organism : str
        name of organism

    Returns
    -------
    str
        URL for fasta file
    """
    ensembl_ftp = config['global']['ensembl_ftp']
    release = str(config['global']['release'])
    fasta_file = config['organisms'][organism]['fasta_file']

    return ensembl_ftp + 'release-' + release + '/fasta/' + organism + '/dna/' + fasta_file


def get_supported_features(organism: str) -> list:
    """Get list of all supported features for given organism specified in a config file.

    Features correspond to different classes of regulatory data stored in Ensembl.
    In generals, this data contain regions that are predicted to regulate gene expression.
    More at http://www.ensembl.org/info/genome/funcgen/index.html.

    Parameters
    ----------
    organism : str
        name of organism

    Returns
    -------
    list
        a list of supported features
    """
    return list(config['organisms'][organism]['features'].keys())


def get_feature_path(organism: str, feature: str) -> str:
    """Get Ensembl URL for txt.gz file with data about feature class of given organism

    Parameters
    ----------
    organism : str
        name of organism
    feature : str
        name of feature class
    Returns
    -------
    str
        URL for txt.gz file
    """
    ensembl_ftp = config['global']['ensembl_ftp']
    release = str(config['global']['release'])
    feature_file = config['organisms'][organism]['features'][feature]['file']

    return ensembl_ftp + 'release-' + release + '/mysql/regulation_mart_' + release + '/' + feature_file


def get_column_names(feature: str) -> list:
    """Get list column names for given feature

    Parameters
    ----------
    feature : str
        name of feature class

    Returns
    -------
    list
        a list of column names
    """
    return config['global']['features'][feature]['columns']


def get_feature_column_name(columns: list) -> str:
    """Get name of column that contains name of feature type

    Each feature class has s lightly different data structure - different column is storing information about type of
    feature.

    Parameters
    ----------
    columns : list
        list of columns

    Returns
    -------
    str
        column containing name of feature type
    """
    feature_column_name = 'so_name'
    if feature_column_name not in columns:
        feature_column_name = 'feature_type_name'
    return feature_column_name


def parse_feature_file(path: str, feature: str) -> pd.DataFrame:
    """Parse file with with feature class data

    Ensembl feature class file contains many unnecessary columns for this purpose. We only need chromosome, start
    position and end position to retrieve corresponding DNA sequence.

    Parameters
    ----------
    path : str
        local path to downloaded file with feature class data
    feature: str
        name of feature class

    Returns
    -------
    pd.DataFrame
        dataframe with parsed info about positions on chromosomes
    """
    logging.info("parse_feature_file(): Going to parse file {}".format(path))

    header_list = get_column_names(feature)
    feature_column_name = get_feature_column_name(header_list)
    needed_cols = [feature_column_name, 'seq_region_name', 'seq_region_start', 'seq_region_end']
    logging.debug("parse_feature_file(): Using headers: {}".format(header_list))
    df = pd.read_csv(path, sep='\t', names=header_list, usecols=needed_cols)

    logging.debug("parse_feature_file(): Kept columns: {}".format(list(df.keys())))
    logging.info("parse_feature_file(): Done parsing file {}".format(path))

    return df


def find_sequences(organism: str, seqs: pd.DataFrame) -> pd.DataFrame:
    """Find DNA sequences for given positions on chromosomes

    Genome of organism is stored in .2bit file (more at https://genome.ucsc.edu/FAQ/FAQformat.html#format7). Input
    *seqs* dataframe is enriched with 'seq' column containing DNA sequence for each given chromosome and position.

    Parameters
    ----------
    organism : str
        name of organism
    seqs: pd.DataFrame
        dataframe with parsed info about positions on chromosomes

    Returns
    -------
    pd.DataFrame
        dataframe with positions on chromosomes and corresponding DNA sequences
    """
    logging.info('find_sequences_and_save_to_fasta(): Going to find sequences based on genomic loci.')

    genome = get_2bit_genome_file(organism)
    num_seqs = len(seqs)

    seqs['seq'] = ''
    seqs['seq_region_name'] = 'chr' + seqs['seq_region_name']
    for i in range(num_seqs):
        chrom = seqs.iloc[i]['seq_region_name']
        if chrom not in genome:
            logging.debug('Chromosome {} not in 2bit genome. Skipping this line.'.format(chrom))
            continue

        start = seqs.iloc[i]['seq_region_start']
        end = seqs.iloc[i]['seq_region_end']
        seqs.at[i, 'seq'] = genome[chrom][start:end]

    logging.info('find_sequences_and_save_to_fasta(): Done finding sequences.')
    return seqs


def make_dataset_from_loci(feature_loci: pd.DataFrame, organism: str, out_dir: str) -> None:
    """Create dataset based on provided loci

    1) Corresponding DNA sequences are retrieved based on loci information from genome of organism.
    2) Low quality sequences are removed.
    3) Remaining sequences are split into train, validation and test parts.
    4) For each DNA sequence, new random sequence is generated. This sequence is from the same genome, but cannot
        intersect with original sequence. We get negative dataset (dataset containing sequences that do not represent
        selected feature) in this way.
    5) Sequences in negative dataset are also split into train, validation and test parts.
    6) Sequences are saved into fasta files. Six files are created:
         - positive_train.fa
         - positive_valid.fa
         - negative_train.da
         - negative_valid.fa
         - test.fa
         - test_with_labels.fa
        File test.fa is for testing purposes and it doesn't contain labels for sequences. To verify your results, you
        can use file test_with_labels.fa.

    Parameters
    ----------
    feature_loci : pd.DataFrame
        dataframe with loci for one feature
    organism: str
        name of organism
    out_dir: str
        path to directory for storing data
    """
    positive_seqs = find_sequences(organism, feature_loci)
    preprocessed_positive_seqs = remove_low_quality(positive_seqs)
    positive_train, positive_val, positive_test = split_train_val_test(preprocessed_positive_seqs)
    save_to_fasta(Path(out_dir, 'positive_train.fa'), positive_train)
    save_to_fasta(Path(out_dir, 'positive_valid.fa'), positive_val)

    negative_seqs = generate_negatives(organism, preprocessed_positive_seqs)
    # we don't need to preprocess negative sequences since we are generating them
    # to match already preprocessed positive sequences
    negative_train, negative_val, negative_test = split_train_val_test(negative_seqs)
    save_to_fasta(Path(out_dir, 'negative_train.fa'), negative_train)
    save_to_fasta(Path(out_dir, 'negative_valid.fa'), negative_val)

    # save test sequences separately - we don't want to expose information about positive/negative
    save_test_to_fasta(Path(out_dir, 'test.fa'), positive_test, negative_test)


def extract_feature_type_loci(seqs: pd.DataFrame, feature_type: str) -> pd.DataFrame:
    """Get loci of selected feature type from dataframe containing loci of all feature types

    Parameters
    ----------
    seqs : pd.DataFrame
        dataframe with loci for all features of some feature class
    feature_type: str
        name of feature

    Returns
    -------
    pd.DataFrame
        dataframe with feature name, chromosome name, start position and end position
    """
    feature_column_name = get_feature_column_name(list(seqs.keys()))
    # get rows for current feature
    feature_loci = seqs[seqs[feature_column_name] == feature_type].copy()
    # drop column with feature name - there is only one value
    feature_loci.drop(feature_column_name, axis=1, inplace=True)
    # reset indexing to go from 0 with step 1 (filtering on feature type left spaces in indices)
    return feature_loci.reset_index(drop=True)


def get_feature_class_loci(organism: str, feature: str, root_dir: str) -> pd.DataFrame:
    """Get loci of selected feature class for selected organism

    Parameters
    ----------
    organism : str
        name of organism
    feature: str
        name of feature class

    Returns
    -------
    pd.DataFrame
        dataframe with feature name, chromosome name, start position and end position
    """
    feature_path = get_feature_path(organism, feature)

    temp_dir = prepare_temp_directory(root_dir)
    local_feature = temp_dir + organism + '_' + feature + '.txt.gz'
    download_file(feature_path, local_feature)

    seqs = parse_feature_file(local_feature, feature)
    # fix for difference in 0/1-based coordinates in retrieved loci and used genome
    seqs['seq_region_start'] = seqs['seq_region_start'] - 1

    return seqs


def make_feature_dataset_for_organism(organism: str, feature: str, root_dir: str):
    seqs = get_feature_class_loci(organism, feature, root_dir)
    feature_column_name = get_feature_column_name(list(seqs.keys()))
    feature_types = seqs[feature_column_name].unique()
    for feature_type in tqdm(feature_types, desc='Processing feature types'):
        out_dir = prepare_data_directory(root_dir, organism, feature, feature_type)
        feature_loci = extract_feature_type_loci(seqs, feature_type)
        make_dataset_from_loci(feature_loci, organism, out_dir)


def verify_input(user_input: str, supported_names) -> list:

    if user_input.strip() == '*':
        return supported_names

    parts = user_input.strip().split(' ')
    for part in parts:
        if part not in supported_names:
            return list()

    return parts


def main():
    ascii_banner = pyfiglet.figlet_format("Ensembl scraper")
    print(ascii_banner)

    user_input = dict()

    while True:

        print("Please select organisms you are interested in. Supported organisms are: "
              + str(get_supported_organisms()))
        organisms = input("Write names of organisms separated by space. Use '*' for all: ")
        organisms = verify_input(organisms, get_supported_organisms())
        if organisms:
            break
        print("Unrecognized organism. Please try again.")
        print('')

    print('')
    print("Organisms selected: " + str(organisms))
    print("Now you can select feature classes for each organism.")

    for organism in organisms:

        while True:

            print('=========================')
            print(organism)
            print('=========================')
            print("Supported feature classes are: " + str(get_supported_features(organism)))
            features = input("Write names of feature_classes separated by space. Use '*' for all: ")
            features = verify_input(features, get_supported_features(organism))
            if features:
                break
            print("Unrecognized feature class. Please try again.")
            print('')

        user_input[organism] = features

    print('')
    print('=========================')
    print('Selected settings')
    print('=========================')
    print(yaml.dump(user_input))
    print('')

    print('One more thing.')

    while True:
        root_dir = input('Please enter path to directory, where this program can create folders and store data: ')
        if Path(root_dir).exists():
            break
        print("Path doesn't exist. Please try again.")
        print('')

    print('')
    input('All done. Confirm by any input to start.')
    print('')

    for o in tqdm(user_input.keys(), desc='Processing organisms'):
        for f in tqdm(user_input[o], desc='Processing feature files'):
            make_feature_dataset_for_organism(o, f, root_dir)


if __name__ == '__main__':
    main()

