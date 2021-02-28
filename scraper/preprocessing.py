from sklearn.model_selection import train_test_split
from pathlib import Path
import numpy as np
from scraper.utils import save_to_fasta
import logging


logging.basicConfig(level=logging.INFO)


def split_train_val_test(data, train_ratio = 0.7, validation_ratio = 0.2, test_ratio = 0.1):

    if abs((train_ratio + validation_ratio + test_ratio) - 1.0) > 0.00000001:
        raise ValueError("Sum of all ratios must be equal to one.")

    # Produces test split.
    x_remaining, x_test = train_test_split(data, test_size=test_ratio, random_state=42)

    # Adjusts val ratio, w.r.t. remaining dataset.
    ratio_remaining = 1 - test_ratio
    ratio_val_adjusted = validation_ratio / ratio_remaining

    # Produces train and val splits.
    x_train, x_val = train_test_split(x_remaining, test_size=ratio_val_adjusted, random_state=42)

    return x_train, x_val, x_test


def split_to_fasta(current_dir, prefix, seqs):

    x_train, x_val, x_test = split_train_val_test(seqs)
    save_to_fasta(Path(current_dir, prefix + '_train.fa'), x_train)
    save_to_fasta(Path(current_dir, prefix + '_valid.fa'), x_val)


    save_to_fasta(Path(current_dir, prefix + '_test.fa'), x_test)


def reject_outliers(seqs, m = 3.):

    lengths = np.array([len(seq.seq) for index,seq in seqs.iterrows()])
    d = np.abs(lengths - np.median(lengths))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return seqs.loc[np.where(s<m)[0]]


def reject_Ns(seqs, threshold=0.05):

    lengths = np.array([len(seq.seq) for index, seq in seqs.iterrows()])
    N_counts = np.array([seq.seq.upper().count('N') for index, seq in seqs.iterrows()])

    return seqs.loc[np.where(N_counts < lengths * threshold)]


def remove_low_quality(seqs):
    logging.info("remove_low_quality(): Going to preprocess sequences.")
    logging.info("remove_low_quality(): Original number of sequences: {}".format(str(len(seqs))))

    seqs = reject_outliers(seqs)
    seqs = seqs.reset_index(drop=True)
    logging.info("remove_low_quality(): Number of sequences after outlier rejection: {}".format(str(len(seqs))))

    seqs = reject_Ns(seqs)
    seqs = seqs.reset_index(drop=True)
    logging.info("remove_low_quality(): Number of sequences after Ns rejection: {}".format(str(len(seqs))))

    logging.info("remove_low_quality(): Done preprocessing sequences.")

    return seqs
