import pandas as pd
from tqdm import tqdm
import logging
from utils import download_file, get_2bit_genome_file, prepare_temp_directory, delete_temp_directory, prepare_data_directory, save_metadata
from config import get_column_names, get_feature_column_name, get_feature_path
from random_negatives import generate_negatives
from preprocessing import remove_low_quality, split_to_csv
from cli import get_user_inputs

logging.basicConfig(level=logging.INFO)


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
    needed_cols = [feature_column_name, 'seq_region_name', 'seq_region_start', 'seq_region_end', 'seq_region_strand']
    logging.debug("parse_feature_file(): Using headers: {}".format(header_list))
    df = pd.read_csv(path, sep='\t', names=header_list, usecols=needed_cols)

    df['seq_region_strand'] = df['seq_region_strand'].map({1: '+', 0: '-'})
    df['seq_region_name'] = 'chr' + df['seq_region_name']

    logging.debug("parse_feature_file(): Kept columns: {}".format(list(df.keys())))
    logging.info("parse_feature_file(): Done parsing file {}".format(path))

    return df


def find_sequences(organism: str, seqs: pd.DataFrame, temp_file:str) -> pd.DataFrame:
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

    genome = get_2bit_genome_file(organism, temp_file)
    num_seqs = len(seqs)

    seqs['seq'] = ''
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


def make_dataset_from_loci(feature_loci: pd.DataFrame, organism: str, out_dir: str, temp_file:str) -> None:
    """Create dataset based on provided loci

    1) Corresponding DNA sequences are retrieved based on loci information from genome of organism.
    2) Low quality sequences are removed.
    3) Remaining sequences are split into train and test parts.
    4) For each DNA sequence, new random sequence is generated. This sequence is from the same genome, but cannot
        intersect with original sequence. We get negative dataset (dataset containing sequences that do not represent
        selected feature) in this way.
    5) Sequences in negative dataset are also split into train and test parts.
    6) Coordinates for sequences are saved into csv files. Folder structure is as follows:
        --- train --- positive.csv
                 |
                  --- negative.csv
        --- test  --- positive.csv
                 |
                  --- negative.csv
        --- metadata.yaml

    Parameters
    ----------
    feature_loci : pd.DataFrame
        dataframe with loci for one feature
    organism: str
        name of organism
    out_dir: str
        path to directory for storing data
    """
    positive_seqs = find_sequences(organism, feature_loci, temp_file)
    preprocessed_positive_seqs = remove_low_quality(positive_seqs)
    split_to_csv(out_dir, "positive", preprocessed_positive_seqs)

    negative_seqs = generate_negatives(organism, preprocessed_positive_seqs)
    # we don't need to preprocess negative sequences since we are generating them
    # to match already preprocessed positive sequences
    split_to_csv(out_dir, "negative", negative_seqs)
    save_metadata(out_dir + 'metadata.yaml', organism)


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


def get_feature_class_loci(organism: str, feature: str, temp_dir: str) -> pd.DataFrame:
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
        dataframe with feature name, chromosome name, start position, end position and strand
    """
    feature_path = get_feature_path(organism, feature)

    local_feature = temp_dir + '/' + organism + '_' + feature + '.txt.gz'
    download_file(feature_path, local_feature)

    seqs = parse_feature_file(local_feature, feature)
    # fix for difference in 0/1-based coordinates in retrieved loci and used genome
    seqs['seq_region_start'] = seqs['seq_region_start'] - 1

    return seqs


def make_feature_dataset_for_organism(organism: str, feature: str, root_dir: str):
    temp_dir = prepare_temp_directory(root_dir)
    seqs = get_feature_class_loci(organism, feature, temp_dir)
    feature_column_name = get_feature_column_name(list(seqs.keys()))
    feature_types = seqs[feature_column_name].unique()
    for feature_type in tqdm(feature_types, desc='Processing feature types'):
        out_dir = prepare_data_directory(root_dir, organism, feature, feature_type)
        feature_loci = extract_feature_type_loci(seqs, feature_type)
        make_dataset_from_loci(feature_loci, organism, out_dir, temp_dir)
    delete_temp_directory(temp_dir)


def main():
    user_input, root_dir = get_user_inputs()

    for o in tqdm(user_input.keys(), desc='Processing organisms'):
        for f in tqdm(user_input[o], desc='Processing feature files'):
            make_feature_dataset_for_organism(o, f, root_dir)


if __name__ == '__main__':
    main()

