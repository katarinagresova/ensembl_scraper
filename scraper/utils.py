from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import logging
from tqdm import tqdm


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