from .utils import get_2bit_genome_file
import pandas as pd
import numpy as np
import re

MIN_CHROM_LENGTH = 100000
np.random.seed(42)


def get_chr_names_and_lengths(organism, local_dir):

    chr_lengths = {}
    genome = get_2bit_genome_file(organism, local_dir)
    for chromosome in genome.keys():
        # working only with classical chromosomes chr1, chr2, .., chrX, chrY, chrMT
        if re.match(r"^chr([1-9][0-9]*|X|Y|MT)$", chromosome):
            length = len(genome[chromosome])
            if length > MIN_CHROM_LENGTH:
                chr_lengths[chromosome] = len(genome[chromosome])

    # check that all lengths are different from 0
    assert all(x != 0 for x in chr_lengths.values())

    return chr_lengths


def get_random_chr(chr_names_and_lengths):
    chr_lengths = pd.Series(chr_names_and_lengths.values())
    chr_probs = chr_lengths / chr_lengths.sum()
    chr_names = list(chr_names_and_lengths.keys())
    return chr_names[np.argwhere(np.random.multinomial(1, chr_probs))[0][0]]


def is_intersecting(c, pos, df_forbidden):
    intersecting = (df_forbidden.seq_region_name == c) & (df_forbidden.seq_region_start.astype(int) <= pos) & (
                df_forbidden.seq_region_end.astype(int) >= pos)
    return intersecting.any()


def get_random_pos(df_forbidden: pd.DataFrame, chr_names_and_lengths, offset_from_end):
    c = get_random_chr(chr_names_and_lengths)
    c_len = chr_names_and_lengths[c]
    pos = np.random.randint(c_len - offset_from_end) + 1

    while is_intersecting(c, pos, df_forbidden):
        pos = np.random.randint(c_len) + 1

    return c, pos


def generate_negatives(organism, excluded_seqs: pd.DataFrame, local_dir):
    chr_names_and_lengths = get_chr_names_and_lengths(organism, local_dir)
    num_seqs = len(excluded_seqs)

    genome = get_2bit_genome_file(organism, local_dir)
    seqs = [None] * num_seqs
    for i in range(num_seqs):
        while True:
            seq_length = int(excluded_seqs['seq_region_end'][i]) - int(excluded_seqs['seq_region_start'][i])
            chrom, start = get_random_pos(excluded_seqs, chr_names_and_lengths, seq_length)
            end = start + seq_length
            seq = genome[chrom][start:end]
            if 'N' not in seq.upper():
                seqs[i] = [seq, chrom, start, end, '+']
                break
    return pd.DataFrame(seqs, columns=['seq', 'seq_region_name', 'seq_region_start', 'seq_region_end', 'seq_region_strand'])
