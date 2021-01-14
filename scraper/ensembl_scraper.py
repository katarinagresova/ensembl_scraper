import yaml
import pandas as pd
from tqdm import tqdm
from twobitreader import twobit_reader
import logging
from scraper.utils import File, download_file, get_2bit_genome_file

logging.basicConfig(level=logging.DEBUG)
CONFIG_FILE = '../config.yaml'
with open(CONFIG_FILE, "r") as ymlfile:
    config = yaml.load(ymlfile, Loader=yaml.FullLoader)


def get_supported_organisms():
    return list(config['organisms'].keys())


def get_fasta_path(organism):
    ensembl_ftp = config['global']['ensembl_ftp']
    release = str(config['global']['release'])
    fasta_file = config['organisms'][organism]['fasta_file']

    return ensembl_ftp + 'release-' + release + '/fasta/' + organism + '/dna/' + fasta_file


def get_2bit_file_name(organism):
    return config['organisms'][organism]['2bit_file_name']


def get_supported_features(organism):
    return list(config['organisms'][organism]['features'].keys())


def get_feature_path(organism, feature):
    ensembl_ftp = config['global']['ensembl_ftp']
    release = str(config['global']['release'])
    feature_file = config['organisms'][organism]['features'][feature]['file']

    return ensembl_ftp + 'release-' + release + '/mysql/regulation_mart_' + release + '/' + feature_file


def get_feature_column_name(columns):
    feature_column_name = 'so_name'
    if feature_column_name not in columns:
        feature_column_name = 'feature_type_name'
    return feature_column_name


def parse_feature_file(path, feature):
    logging.info("parse_feature_file(): Going to parse file {}".format(path))

    header_list = config['global']['features'][feature]['columns']
    feature_column_name = get_feature_column_name(header_list)
    needed_cols = [feature_column_name, 'seq_region_name', 'seq_region_start', 'seq_region_end']
    logging.debug("parse_feature_file(): Using headers: {}".format(header_list))
    df = pd.read_csv(path, sep='\t', names=header_list, usecols=needed_cols)

    logging.debug("parse_feature_file(): Kept columns: {}".format(list(df.keys())))
    logging.info("parse_feature_file(): Done parsing file {}".format(path))

    return df


def convert_df_to_expected_bed(df):
    seqs_loci = df[['seq_region_name', 'seq_region_start', 'seq_region_end']].values.tolist()
    return ['chr' + ' '.join(str(x) for x in line) for line in seqs_loci]


def find_sequences_and_save_to_fasta(organism, seqs, out_fasta, local_dir='../../ensembl_data/2bit/'):

    feature_column_name = get_feature_column_name(seqs.keys())
    feature_type = seqs[feature_column_name][0]
    logging.info('find_sequences_and_save_to_fasta(): '
                 'Going to find sequences based on genomic loci and save results to fasta file. '
                 'Organism: {}, Feature type: {}'.format(organism, feature_type))

    # fix for difference in 0/1-based coordinates in retrieved loci and used genome
    seqs['seq_region_start'] = seqs['seq_region_start'] - 1

    seqs_loci_list = convert_df_to_expected_bed(seqs)
    fasta_handle = File(out_fasta)
    genome = get_2bit_genome_file(organism, local_dir)
    twobit_reader(genome, seqs_loci_list, fasta_handle.write)

    logging.info('find_sequences_and_save_to_fasta(): Done finding sequences and saving them to fasta file.')


if __name__ == '__main__':

    organisms = get_supported_organisms()
    for o in tqdm(organisms, desc='Processing organisms'):

        #features = get_supported_features(o)
        features = ['regulatory_feature', 'mirna_target_feature', 'external_feature']
        for f in tqdm(features, desc='Processing feature files'):

            feature_path = get_feature_path(o, f)
            local_feature = '../../ensembl_data/feature/' + o + '_' + f + '.txt.gz'
            download_file(feature_path, local_feature)

            seqs = parse_feature_file(local_feature, f)

            feature_column_name = get_feature_column_name(seqs.keys())
            feature_types = seqs[feature_column_name].unique()
            for feature_type in tqdm(feature_types, desc='Processing feature types'):

                feature_seqs = seqs[seqs[feature_column_name] == feature_type].copy()
                feature_seqs = feature_seqs.reset_index(drop=True)

                find_sequences_and_save_to_fasta(o, feature_seqs, '../../ensembl_data/result/' + o + '_' + f + '_' + feature_type + '.fa')
