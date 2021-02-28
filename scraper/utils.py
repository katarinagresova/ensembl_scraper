import os
import urllib.request
import yaml
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging
from twobitreader import TwoBitFile
from twobitreader.download import save_genome
from pathlib import Path
import pandas as pd

CONFIG_FILE = '../config.yaml'
with open(CONFIG_FILE, "r") as ymlfile:
    config = yaml.load(ymlfile, Loader=yaml.FullLoader)


def make_dir(dir):
  Path(dir).mkdir(parents=True, exist_ok=True)


def save_test_to_fasta(filename, positives, negatives):
    positives.loc[:, 'positive'] = 1
    negatives.loc[:, 'positive'] = 0
    samples = pd.concat([positives, negatives])

    # shuffle samples
    samples = samples.sample(frac=1).reset_index(drop=True)

    # save positives+negatives without label
    save_to_fasta(filename, samples)

    samples.loc[samples['positive'] == 1, 'seq_region_name'] = 'positive_' + samples['seq_region_name']
    samples.loc[samples['positive'] == 0, 'seq_region_name'] = 'negative_' + samples['seq_region_name']

    # change filename of file with labels
    filename = filename.with_name(filename.stem + '_with_labels.fa')
    # save positives+negatives with label at beginning of ID
    save_to_fasta(filename, samples)


def save_to_fasta(filename, sequences):
    if isinstance(sequences, SeqRecord) or isinstance(sequences, list):
        save_seq_to_fasta(filename, sequences)
    elif isinstance(sequences, pd.DataFrame):
        save_df_to_fasta(filename, sequences)
    else:
        raise TypeError("Unsupported format of sequences. Use SeqRecord, list of SeqRecords or pd.DataFrame.")


def save_seq_to_fasta(filename, seqs):
    with open(filename, 'w') as handle:
        SeqIO.write(seqs, handle, 'fasta')


def save_df_to_fasta(filename, seq_df):

    logging.info("save_to_fasta(): Going to save sequences to file {}".format(filename))

    with open(filename, 'w') as handle:
        for index, record in seq_df.iterrows():

            # check if record contains sequence
            #  - there might be some loci from scaffolds that was not found in fasta file
            if record.seq != '':
                SeqIO.write(
                    SeqRecord(
                        Seq(record.seq),
                            record.seq_region_name + ":" + str(record.seq_region_start) + ".." + str(
                            record.seq_region_end),
                        description=""
                    ),
                    handle,
                    'fasta'
                )

    logging.info("save_to_fasta(): Sequences saved to file {}".format(filename))


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


def get_2bit_file_name(organism):
    return config['organisms'][organism]['2bit_file_name']
