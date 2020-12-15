import yaml
import urllib.request
import pandas as pd
from tqdm import tqdm
from twobitreader import TwoBitFile, twobit_reader
from twobitreader.download import save_genome
import logging
import os

logging.basicConfig(level=logging.DEBUG)
CONFIG_FILE = '../config.yaml'
with open(CONFIG_FILE, "r") as ymlfile:
    config = yaml.load(ymlfile, Loader=yaml.FullLoader)


class File:
    def __init__(self, name, mode='w'):
        self.f = open(name, mode, buffering=1)

    def write(self, string, newline=True):
        if newline:
            self.f.write(string + '\n')
        else:
            self.f.write(string)


def download_file(url, local_path):
    logging.info("download_file(): Going to download file from path {}".format(url))

    try:
        file = open(local_path)
        file.close()
        logging.info("download_file(): File {} already exists. Not going to download.".format(local_path))
    except FileNotFoundError:
        urllib.request.urlretrieve(url, local_path)
        logging.info("download_file(): File downloaded to path {}.".format(local_path))


def download_2bit_file(genome_name, local_dir):
    logging.info("download_2bit_file(): Going to download 2bit file {}".format(genome_name))

    try:
        file = open(os.path.join(local_dir, genome_name + '.2bit'))
        file.close()
        logging.info("download_2bit_file(): File for {} already exists. Not going to download.".format(genome_name))
    except FileNotFoundError:
        save_genome(genome_name, destdir=local_dir)
        logging.info("download_2bit_file(): File for {} downloaded to path {}.".format(genome_name, os.path.join(local_dir, genome_name + '.2bit')))


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
    logging.debug("parse_feature_file(): Using headers: {}".format(header_list))
    df = pd.read_csv(path, sep='\t', names=header_list)
    feature_column_name = get_feature_column_name(df.keys())
    seqs = df[[feature_column_name, 'seq_region_name', 'seq_region_start', 'seq_region_end']].copy()
    seqs['seq'] = ''
    logging.debug("parse_feature_file(): Kept columns: {}".format(list(seqs.keys())))
    logging.info("parse_feature_file(): Done parsing file {}".format(path))

    return seqs


def convert_df_to_expected_bed(df):
    seqs_loci = df[['seq_region_name', 'seq_region_start', 'seq_region_end']].values.tolist()
    return ['chr' + ' '.join(str(x) for x in line) for line in seqs_loci]


def find_sequences_and_save_to_fasta(organism, seqs, out_fasta, local_dir='../../ensembl_data/2bit/'):

    logging.info('find_sequences_and_save_to_fasta(): Going to find sequences based on genomic loci and save results to fasta file.')

    genome_name = get_2bit_file_name(organism)
    download_2bit_file(genome_name, local_dir)
    twobit_path = os.path.join(local_dir, genome_name + '.2bit')
    seqs_loci_list = convert_df_to_expected_bed(seqs)
    genome = TwoBitFile(twobit_path)
    fasta_handle = File(out_fasta)
    twobit_reader(genome, seqs_loci_list, fasta_handle.write)

    logging.info('find_sequences_and_save_to_fasta(): Done finding sequences and saving them to fasta file.')


if __name__ == '__main__':

    organisms = get_supported_organisms()
    for o in organisms:

        features = get_supported_features(o)
        for f in features:

            feature_path = get_feature_path(o, f)
            local_feature = '../../ensembl_data/feature/' + o + '_' + f + '.txt.gz'
            download_file(feature_path, local_feature)

            seqs = parse_feature_file(local_feature, f)

            feature_column_name = get_feature_column_name(seqs.keys())
            feature_types = seqs[feature_column_name].unique()
            for feature_type in feature_types:

                feature_seqs = seqs[seqs[feature_column_name] == feature_type].copy()
                feature_seqs = feature_seqs.reset_index(drop=True)

                find_sequences_and_save_to_fasta(o, feature_seqs, '../../ensembl_data/result/' + o + '_' + f + '_' + feature_type + '.fa')
