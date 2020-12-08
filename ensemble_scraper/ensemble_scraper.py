import yaml
import urllib.request
import pandas as pd
import gzip
from tqdm import tqdm
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging

logging.basicConfig(level=logging.DEBUG)
CONFIG_FILE = 'config.yaml'
with open(CONFIG_FILE, "r") as ymlfile:
    config = yaml.load(ymlfile, Loader=yaml.FullLoader)


def download_file(url, local_path):
    logging.info("download_file(): Going to download file from path {}".format(url))

    try:
        file = open(local_path)
        file.close()
        logging.info("download_file(): File {} already exists. Not going to download.".format(local_path))
    except FileNotFoundError:
        urllib.request.urlretrieve(url, local_path)
        logging.info("download_file(): File downloaded to path {}.".format(local_path))


def get_supported_organisms():
    return list(config['organisms'].keys())


def get_fasta_path(organism):
    ensembl_ftp = config['global']['ensembl_ftp']
    release = str(config['global']['release'])
    fasta_file = config['organisms'][organism]['fasta_file']

    return ensembl_ftp + 'release-' + release + '/fasta/' + organism + '/dna/' + fasta_file


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


def find_sequences(fasta_file, seqs):
    logging.info("find_sequences(): Going to find sequences in file {} for features {}".format(
        fasta_file,
        seqs[get_feature_column_name(seqs.keys())].unique()
    ))

    def which(self):
        try:
            self = list(iter(self))
        except TypeError as e:
            raise Exception("""'which' method can only be applied to iterables.
            {}""".format(str(e)))
        indices = [i for i, x in enumerate(self) if bool(x) is True]
        return indices

    with gzip.open(fasta_file, "rt") as handle:
        for record in tqdm(SeqIO.parse(handle, "fasta"), total=24):
            sel_seqs = which(seqs.seq_region_name == record.id)
            for i in sel_seqs:
                seqs.loc[i, "seq"] = str(record.seq[(seqs.seq_region_start[i] - 1):seqs.seq_region_end[i]])

            if record.id == "MT":
                # stop, do not read small contigs
                break

    logging.info("find_sequences(): Found sequences from file {}".format(fasta_file))
    return seqs


def save_to_fasta(filename, seq_df):

    logging.info("save_to_fasta(): Going to save sequences to file {}".format(filename))

    with open(filename, 'w') as handle:
        for index, record in tqdm(seq_df.iterrows()):

            # check if record contains sequence
            #  - there might be some loci from scaffolds that was not found in fasta file
            if record.seq != '':
                SeqIO.write(
                    SeqRecord(
                        Seq(record.seq),
                        record.so_name + '_' + record.seq_region_name + ":" + str(record.seq_region_start) + ".." + str(
                            record.seq_region_end),
                        description=""
                    ),
                    handle,
                    'fasta'
                )

    logging.info("save_to_fasta(): Sequences saved to file {}".format(filename))


if __name__ == '__main__':

    organisms = get_supported_organisms()
    for o in organisms:

        fasta_path = get_fasta_path(o)
        local_fasta = '../data/fasta/' + o + '.fa.gz'
        download_file(fasta_path, local_fasta)

        features = get_supported_features(o)
        for f in features:

            feature_path = get_feature_path(o, f)
            local_feature = '../data/feature/' + o + '_' + f + '.txt.gz'
            download_file(feature_path, local_feature)

            seqs = parse_feature_file(local_feature, f)

            feature_column_name = get_feature_column_name(seqs.keys())
            feature_types = seqs[feature_column_name].unique()
            for feature_type in feature_types:

                feature_seqs = seqs[seqs[feature_column_name] == feature_type].copy()
                feature_seqs = feature_seqs.reset_index(drop=True)

                # TODO: looking for sequences might be better to do on whole feature file at once
                feature_seqs = find_sequences(local_fasta, feature_seqs)
                save_to_fasta('../data/result/' + o + '_' + f + '_' + feature_type + '.fa', feature_seqs)
