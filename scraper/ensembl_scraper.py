import pandas as pd
from tqdm import tqdm
import logging
from pathlib import Path
from utils import download_file, get_2bit_genome_file, make_dir, config, save_to_fasta, save_test_to_fasta
from random_negatives import generate_negatives
from preprocessing import remove_low_quality, split_train_val_test

logging.basicConfig(level=logging.INFO)


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
    feature_column_name = get_feature_column_name(header_list)
    needed_cols = [feature_column_name, 'seq_region_name', 'seq_region_start', 'seq_region_end']
    logging.debug("parse_feature_file(): Using headers: {}".format(header_list))
    df = pd.read_csv(path, sep='\t', names=header_list, usecols=needed_cols)

    logging.debug("parse_feature_file(): Kept columns: {}".format(list(df.keys())))
    logging.info("parse_feature_file(): Done parsing file {}".format(path))

    return df


def find_sequences(organism, seqs):

    logging.info('find_sequences_and_save_to_fasta(): Going to find sequences based on genomic loci.')

    genome = get_2bit_genome_file(organism)
    num_seqs = len(seqs)

    seqs['seq'] = ''
    seqs['seq_region_name'] = 'chr' + seqs['seq_region_name']
    for i in range(num_seqs):
        chrom = seqs.iloc[i]['seq_region_name']
        if chrom not in genome:
            logging.debug('Chromosome {} not in 2bit genome. Skipping this line.'.format(chrom))
            continue

        start = seqs.iloc[i]['seq_region_start']
        end = seqs.iloc[i]['seq_region_end']
        seqs.at[i, 'seq'] = genome[chrom][start:end]

    logging.info('find_sequences_and_save_to_fasta(): Done finding sequences.')
    return seqs


if __name__ == '__main__':

    organisms = get_supported_organisms()
    for o in tqdm(organisms, desc='Processing organisms'):

        features = ['regulatory_feature', 'mirna_target_feature', 'external_feature']
        for f in tqdm(features, desc='Processing feature files'):

            feature_path = get_feature_path(o, f)
            local_feature = '../../ensembl_data/feature/' + o + '_' + f + '.txt.gz'
            download_file(feature_path, local_feature)

            seqs = parse_feature_file(local_feature, f)
            # fix for difference in 0/1-based coordinates in retrieved loci and used genome
            seqs['seq_region_start'] = seqs['seq_region_start'] - 1

            feature_column_name = get_feature_column_name(seqs.keys())
            feature_types = seqs[feature_column_name].unique()
            for feature_type in tqdm(feature_types, desc='Processing feature types'):

                # get rows for current feature
                feature_seqs = seqs[seqs[feature_column_name] == feature_type].copy()
                # drop column with feature name - there is only one value
                feature_seqs.drop(feature_column_name, axis=1, inplace=True)
                # reset indexing to go from 0 with step 1 (filtering on feature type left spaces in indices)
                feature_seqs = feature_seqs.reset_index(drop=True)

                out_dir = '../../ensembl_data_test/result/' + o + '/' + f + '_' + feature_type + '/'
                make_dir(out_dir)

                positive_seqs = find_sequences(o, feature_seqs)
                preprocessed_positive_seqs = remove_low_quality(positive_seqs)
                positive_train, positive_val, positive_test = split_train_val_test(preprocessed_positive_seqs)
                save_to_fasta(Path(out_dir, 'positive_train.fa'), positive_train)
                save_to_fasta(Path(out_dir, 'positive_valid.fa'), positive_val)

                negative_seqs = generate_negatives(o, preprocessed_positive_seqs)
                # we don't need to preprocess negative sequences since we are generating them
                # to match already preprocessed positive sequences
                negative_train, negative_val, negative_test = split_train_val_test(negative_seqs)
                save_to_fasta(Path(out_dir, 'negative_train.fa'), negative_train)
                save_to_fasta(Path(out_dir, 'negative_valid.fa'), negative_val)

                # save test sequences separately - we don't want to expose information about positive/negative
                save_test_to_fasta(Path(out_dir, 'test.fa'), positive_test, negative_test)
