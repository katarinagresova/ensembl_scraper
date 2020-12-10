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