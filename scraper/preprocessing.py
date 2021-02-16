from sklearn.model_selection import train_test_split
from Bio import SeqIO
from pathlib import Path
import numpy as np
from scraper.utils import save_to_fasta
import logging


logging.basicConfig(level=logging.DEBUG)


def split_train_val_test(data, train_ratio = 0.7, validation_ratio = 0.2, test_ratio = 0.1):

    # Produces test split.
    x_remaining, x_test = train_test_split(data, test_size=test_ratio, random_state=42)

    # Adjusts val ratio, w.r.t. remaining dataset.
    ratio_remaining = 1 - test_ratio
    ratio_val_adjusted = validation_ratio / ratio_remaining

    # Produces train and val splits.
    x_train, x_val = train_test_split(x_remaining, test_size=ratio_val_adjusted, random_state=42)

    return x_train, x_val, x_test


def split_fasta(current_dir, prefix, seqs):

    x_train, x_val, x_test = split_train_val_test(seqs)
    save_to_fasta(Path(current_dir, prefix + '_train.fa'), x_train)
    save_to_fasta(Path(current_dir, prefix + '_valid.fa'), x_val)
    save_to_fasta(Path(current_dir, prefix + '_test.fa'), x_test)


def reject_outliers(seqs, m = 3.):
    lengths = np.array([len(seq) for seq in seqs])
    d = np.abs(lengths - np.median(lengths))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return seqs[np.where(s<m)]


def reject_Ns(seqs, threshold = 0.05):

    lengths = np.array([len(seq) for seq in seqs])
    N_counts = np.array([seq.seq.upper().count('N') for seq in seqs])

    return seqs[np.where(N_counts < lengths * threshold)]


def remove_low_quality(file):
    logging.info("remove_low_quality(): Going to preprocess file {}".format(file))

    with open(file, 'r') as handle:
        seqs = np.array([record for record in SeqIO.parse(handle, 'fasta')], dtype=object)

    logging.info("remove_low_quality(): Original number of sequences: {}".format(str(len(seqs))))

    seqs = reject_outliers(seqs)
    logging.info("remove_low_quality(): Number of sequences after outlier rejection: {}".format(str(len(seqs))))

    seqs = reject_Ns(seqs)
    logging.info("remove_low_quality(): Number of sequences after Ns rejection: {}".format(str(len(seqs))))

    logging.info("remove_low_quality(): Done preprocessing file {}".format(file))

    return seqs
