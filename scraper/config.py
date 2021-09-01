import yaml

CONFIG_FILE = '../config.yaml'
with open(CONFIG_FILE, "r") as ymlfile:
    config = yaml.load(ymlfile, Loader=yaml.FullLoader)

def get_supported_organisms() -> list:
    """Get list of all supported organisms specified in a config file.

    Returns
    -------
    list
        a list of supported organisms
    """
    return list(config['organisms'].keys())


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
    ensembl_ftp = config['global']['ensembl_ftp']
    release = str(config['global']['release'])
    fasta_file = config['organisms'][organism]['fasta_file']

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
    return list(config['organisms'][organism]['features'].keys())


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
    ensembl_ftp = config['global']['ensembl_ftp']
    release = str(config['global']['release'])
    feature_file = config['organisms'][organism]['features'][feature]['file']

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
    return config['global']['features'][feature]['columns']


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
