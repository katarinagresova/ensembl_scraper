from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from scraper.utils import get_2bit_genome_file
import pandas as pd
import numpy as np
import os

MIN_CHROM_LENGTH = 100000


def get_chr_names_and_lengths(organism):

    chr_lengths = {}
    genome = get_2bit_genome_file(organism)
    for chromosome in genome.keys():
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


def generate_negatives_and_save_to_fasta(organism, excluded_seqs: pd.DataFrame, out_dir, local_dir='../../ensembl_data/2bit/'):
    chr_names_and_lengths = get_chr_names_and_lengths(organism)
    num_seqs = len(excluded_seqs)

    genome = get_2bit_genome_file(organism)
    with open(out_dir + 'negative.fa', 'w') as handle:
        try:
            for i in range(num_seqs):
                while True:
                    seq_length = int(excluded_seqs['seq_region_end'][i]) - int(excluded_seqs['seq_region_start'][i])
                    chrom, start = get_random_pos(excluded_seqs, chr_names_and_lengths, seq_length)
                    end = start + seq_length
                    seq = genome[chrom][start:end]
                    if 'N' not in seq.upper():
                        SeqIO.write(
                            SeqRecord(
                                Seq(seq),
                                chrom + ':' + str(start) + '-' + str(end),
                                description=""
                            ),
                            handle,
                            'fasta'
                        )
                    break
        except:
            # close and delete the file, so we don't have partial results that look like full results
            # the ctx manager will try to close the file again, but that's harmless.
            handle.close()
            os.remove(out_dir + 'negative.fa')
            raise
