import os
import urllib.request
import yaml
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging
from tqdm import tqdm
from twobitreader import TwoBitFile
from twobitreader.download import save_genome
from pathlib import Path

CONFIG_FILE = '../config.yaml'
with open(CONFIG_FILE, "r") as ymlfile:
    config = yaml.load(ymlfile, Loader=yaml.FullLoader)


def make_dir(dir):
  Path(dir).mkdir(parents=True, exist_ok=True)


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


class File:
    """Wrapper class around basic python file handling.

    Created to override write() method and add newline at the end of every write()
    """

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


def get_2bit_genome_file(organism, local_dir='../../ensembl_data/2bit/'):
    genome_name = get_2bit_file_name(organism)
    download_2bit_file(genome_name, local_dir)
    twobit_path = os.path.join(local_dir, genome_name + '.2bit')
    return TwoBitFile(twobit_path)


def convert_df_to_format_for_twobitreader(df):
    seqs_loci = df[['seq_region_name', 'seq_region_start', 'seq_region_end']].values.tolist()

    if df['seq_region_name'][0].startswith('chr'):
        return [' '.join(str(x) for x in line) for line in seqs_loci]
    else:
        return ['chr' + ' '.join(str(x) for x in line) for line in seqs_loci]


def get_2bit_file_name(organism):
    return config['organisms'][organism]['2bit_file_name']
