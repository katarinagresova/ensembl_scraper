import os
import urllib.request
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging
from twobitreader import TwoBitFile
from twobitreader.download import save_genome
from pathlib import Path
import pandas as pd
from .config import get_2bit_file_name, get_fasta_path
import yaml


def make_dir(dir):
    Path(dir).mkdir(parents=True, exist_ok=True)


def prepare_data_directory(root_dir: str, organism: str, feature: str, feature_type: str) -> str:
    out_dir = root_dir + organism + '_' + feature + '_' + feature_type + '/'
    make_dir(out_dir + "/train")
    make_dir(out_dir + "/test")
    return out_dir


def delete_temp_directory(temp_dir: str) -> None:
    path = Path(temp_dir) # allow path to be a string
    assert path.is_dir() # make sure it`s a folder
    for p in reversed(list(path.glob('**/*'))): # iterate contents from leaves to root
        if p.is_file():
           p.unlink()
        elif p.is_dir():
            p.rmdir()
    path.rmdir()


def prepare_temp_directory(root_dir: str) -> str:
    temp_dir = root_dir + '/tmp/'
    make_dir(temp_dir)
    return temp_dir


def save_to_csv(path: Path, seqs: pd.DataFrame):
    seqs = seqs.reset_index(drop=True)
    seqs.to_csv(
        path_or_buf=path,
        columns=['seq_region_name', 'seq_region_start', 'seq_region_end', 'seq_region_strand'],
        header=['region', 'start', 'end', 'strand'],
        index_label='id'
    )


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
        save_genome(genome_name, destdir=local_dir, mode='http')
        logging.info("download_2bit_file(): File for {} downloaded to path {}.".format(genome_name, os.path.join(local_dir, genome_name + '.2bit')))


def get_2bit_genome_file(organism, local_dir):
    genome_name = get_2bit_file_name(organism)
    download_2bit_file(genome_name, local_dir)
    twobit_path = os.path.join(local_dir, genome_name + '.2bit')
    return TwoBitFile(twobit_path)


def save_metadata(filename, organism):
    metadata = {
        'version': 0,
        'classes': {
            'positive': {
                'type': 'fa.gz',
                'url': get_fasta_path(organism)
            },
            'negative': {
                'type': 'fa.gz',
                'url': get_fasta_path(organism)
            }
        }
    }

    with open(filename, 'w') as handle:
        yaml.dump(metadata, handle)
