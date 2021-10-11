CONFIG = {
    "global": {
        "ensembl_ftp": "ftp://ftp.ensembl.org/pub/",
        "release": 100,
        "features": {
            "regulatory_feature": {
                "columns": [
                    "feature_type_description",
                    "bound_seq_region_start",
                    "feature_type_name",
                    "bound_seq_region_end",
                    "so_name",
                    "so_accession",
                    "seq_region_name",
                    "stable_id",
                    "regulatory_feature_id",
                    "seq_region_start",
                    "seq_region_end",
                    "seq_region_strand"
                ]
            },
            "motif_feature": {
                "columns": [
                    "score",
                    "seq_region_name",
                    "binding_matrix_id",
                    "seq_region_end",
                    "feature_type_name",
                    "motif_feature_id",
                    "display_label",
                    "seq_region_start",
                    "seq_region_strand"
                ]
            },
            "mirna_target_feature": {
                "columns": [
                    "display_label",
                    "feature_type_class",
                    "seq_region_end",
                    "accession",
                    "evidence",
                    "so_name",
                    "seq_region_start",
                    "so_accession",
                    "mirna_target_feature_id",
                    "gene_stable_id",
                    "seq_region_name",
                    "seq_region_strand",
                    "feature_type_description"
                ]
            },
            "external_feature": {
                "columns": [
                    "so_name",
                    "feature_type_class",
                    "so_accession",
                    "seq_region_strand",
                    "seq_region_start",
                    "display_label",
                    "seq_region_end",
                    "fs_display_label",
                    "external_feature_id",
                    "feature_type_description",
                    "seq_region_name"
                ]
            }
        }
    },
    "organisms": {
        "mus_musculus": {
            "fasta_file": "Mus_musculus.GRCm38.dna_rm.toplevel.fa.gz",
            "2bit_file_name": "mm10",
            "features": {
                "regulatory_feature": {
                    "file": "mmusculus_regulatory_feature__regulatory_feature__main.txt.gz"
                },
                "motif_feature": {
                    "file": "mmusculus_motif_feature__motif_feature__main.txt.gz"
                },
                "mirna_target_feature": {
                    "file": "mmusculus_mirna_target_feature__mirna_target_feature__main.txt.gz"
                },
                "external_feature": {
                    "file": "mmusculus_external_feature__external_feature__main.txt.gz"
                }
            }
        },
        "homo_sapiens": {
            "fasta_file": "Homo_sapiens.GRCh38.dna.toplevel.fa.gz",
            "2bit_file_name": "hg38",
            "features": {
                "regulatory_feature": {
                    "file": "hsapiens_regulatory_feature__regulatory_feature__main.txt.gz"
                },
                "motif_feature": {
                    "file": "hsapiens_motif_feature__motif_feature__main.txt.gz"
                },
                "mirna_target_feature": {
                    "file": "hsapiens_mirna_target_feature__mirna_target_feature__main.txt.gz"
                },
                "external_feature": {
                    "file": "hsapiens_external_feature__external_feature__main.txt.gz"
                }
            }
        }
    }
}


def get_supported_organisms() -> list:
    """Get list of all supported organisms specified in a config file.

    Returns
    -------
    list
        a list of supported organisms
    """
    return list(CONFIG['organisms'].keys())


def get_fasta_path(organism: str) -> str:
    """Get Ensembl URL for fasta file of given organism

    Parameters
    ----------
    organism : str
        name of organism

    Returns
    -------
    str
        URL for fasta file
    """
    ensembl_ftp = CONFIG['global']['ensembl_ftp']
    release = str(CONFIG['global']['release'])
    fasta_file = CONFIG['organisms'][organism]['fasta_file']

    return ensembl_ftp + 'release-' + release + '/fasta/' + organism + '/dna/' + fasta_file


def get_supported_features(organism: str) -> list:
    """Get list of all supported features for given organism specified in a config file.

    Features correspond to different classes of regulatory data stored in Ensembl.
    In generals, this data contain regions that are predicted to regulate gene expression.
    More at http://www.ensembl.org/info/genome/funcgen/index.html.

    Parameters
    ----------
    organism : str
        name of organism

    Returns
    -------
    list
        a list of supported features
    """
    return list(CONFIG['organisms'][organism]['features'].keys())


def get_feature_path(organism: str, feature: str) -> str:
    """Get Ensembl URL for txt.gz file with data about feature class of given organism

    Parameters
    ----------
    organism : str
        name of organism
    feature : str
        name of feature class
    Returns
    -------
    str
        URL for txt.gz file
    """
    ensembl_ftp = CONFIG['global']['ensembl_ftp']
    release = str(CONFIG['global']['release'])
    feature_file = CONFIG['organisms'][organism]['features'][feature]["file"]

    return ensembl_ftp + 'release-' + release + '/mysql/regulation_mart_' + release + '/' + feature_file


def get_column_names(feature: str) -> list:
    """Get list column names for given feature

    Parameters
    ----------
    feature : str
        name of feature class

    Returns
    -------
    list
        a list of column names
    """
    return CONFIG['global']['features'][feature]['columns']


def get_feature_column_name(columns: list) -> str:
    """Get name of column that contains name of feature type

    Each feature class has s lightly different data structure - different column is storing information about type of
    feature.

    Parameters
    ----------
    columns : list
        list of columns

    Returns
    -------
    str
        column containing name of feature type
    """
    feature_column_name = 'so_name'
    if feature_column_name not in columns:
        feature_column_name = 'feature_type_name'
    return feature_column_name


def get_2bit_file_name(organism):
    return CONFIG['organisms'][organism]['2bit_file_name']
