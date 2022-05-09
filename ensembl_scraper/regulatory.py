from fastdownload import FastDownload
from pathlib import Path
import gzip
import shutil
import sqlparse
import pandas as pd
import re

class Metadata():

    releases = []

    def __init__(self, release):

        self.release = self._verify_release(release)

        # prepare downloader
        self.downloader = FastDownload(base='~/.ensembl_scraper/{}'.format(self.release))
        self._register_gz()

        self._prepare_metadata()

        self.feature_data = {}


    ### GETTERS ###

    @classmethod
    def get_releases(cls):

        if not cls.releases:
            d = FastDownload(base='~/.ensembl_scraper')
            index_path = d.download('http://ftp.ensembl.org/pub/')

            with open(index_path, 'r') as f:
                releases = []
                for line in f.readlines():
                    line = line.strip("\n")
                    r = re.match(r"^.*>release-(?P<release>\d+)/</a>.*$", line)
                    if r:
                        releases.append(r.group('release'))

            cls.releases = releases

        return cls.releases

    def get_release(self):
        return self.release

    def get_organisms(self):

        organisms = []

        for short_name in self.datasets_info['dataset'].apply(lambda x: x.split('_')[0]):
            if short_name in self.organism_info['short_name'].values:
                organisms.append(self.organism_info[self.organism_info['short_name'] == short_name]['sql_name'].values[0])

        return set(organisms)

    def get_feature_classes(self, organism):

        feature_classes = []

        for dataset in self.datasets_info['dataset'].values:
            if dataset.startswith(self.organism_info[self.organism_info['sql_name'] == organism]['short_name'].values[0]):
                feature_classes.append(dataset[dataset.index('_')+1:])

        return set(feature_classes)

    def get_features(self, organism, feature_class):

        path = self._get_feature_class_path(organism, feature_class)

        if not Path(path.parent / organism / feature_class).exists():

            header_list = self._get_column_names(organism, feature_class)
            needed_cols = ['so_term', 'seq_region_name', 'seq_region_start', 'seq_region_end', 'seq_region_strand']
            df = pd.read_csv(path, sep='\t', names=header_list, usecols=needed_cols, low_memory=False)

            # mapping strand symbols to more widely used - (+, -)
            df['seq_region_strand'] = df['seq_region_strand'].map({1: '+', -1: '-'})
            # mapping chromosome names to format chr([1-9][0-9]*|X|Y|MT) - following 2bit file uses this format
            df['seq_region_name'] = df['seq_region_name'].apply(lambda x: 'chr' + str(x))
            # fix for difference in 0/1-based coordinates in retrieved loci and used genome
            df['seq_region_start'] = df['seq_region_start'] - 1

            features = df.unique() 

            for feature in features:
                df[df['so_term'] == feature].to_csv(Path(path.parent / organism / feature_class / feature + '.csv'), index=False)

        else:
            features = [f.stem for f in Path(path.parent / organism / feature_class).glob('*')]

        return features

    def get_full_name(self, organisms):
        """Map organisms sql_name to full_name

        Args:
            organisms (str or List): _description_
        """

        if type(organisms) is list:
            organism_mapping = self._prepare_organisms_info()
            return organism_mapping[organism_mapping['sql_name'].isin(organisms)]['full_name'].values
        elif type(organisms) is str:
            organism_mapping = self._prepare_organisms_info()
            return organism_mapping[organism_mapping['sql_name'] == organisms]['full_name'].values[0]
        else:
            raise TypeError("organisms must be either a list or a string")

    ### PRIVATE METHODS ###

    def _verify_release(self, release):

        if release == 'latest':
            return max([int(rel) for rel in self.get_releases()])

        if release not in self.get_releases():
            raise ValueError("Release {} is not available".format(release))
        return str(release)

    def _prepare_metadata(self):

        self.tables = self._prepare_tables()
        self.organism_info = self._prepare_organisms_info()
        self.datasets_info = self._prepare_datasets_info()

    def _register_gz(self):
            """Register gz files in the downloader"""

            def gunzip_something(gzipped_file_name, work_dir):
                """gunzip the given gzipped fil

                Args:
                    gzipped_file_name (str): path to the gzipped file
                    work_dir (str): path to the directory where the file will be unzipped
                """

                filename = Path(gzipped_file_name).stem

                with gzip.open(gzipped_file_name, 'rb') as f_in:
                    with open(Path(work_dir, filename), 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)



            if 'gz' not in shutil._UNPACK_FORMATS.keys():
                shutil.register_unpack_format(
                    name='gz',
                    extensions=['.gz'],
                    function=gunzip_something,
                    description='Gzipped file'
                )

    def _prepare_tables(self):

        tables_url = "ftp://ftp.ensembl.org/pub/release-{}/mysql/regulation_mart_{}/regulation_mart_{}.sql.gz".format(self.release, self.release, self.release)
        
        tables_path = self.downloader.get(tables_url)
        tables = self._parse_tables_file(tables_path)

        return tables


    def _parse_tables_file(self, tables_path):

        def get_table_name(tokens):
            for token in reversed(tokens):
                if token.ttype is None:
                    return token.value
            return " "

        def get_column_names(token):

            txt = token.value
            columns = txt[1:txt.rfind(")")].split(",\n")
            column_names = []
            for column in columns:
                c = ' '.join(column.split()).split()
                c_name = c[0].replace('\"',"").strip('`')
                if c_name != "KEY" and c_name != "UNIQUE":
                    m = re.search(r"\d", c_name)
                    if m:
                        c_name = c_name[:m.start() - 1]
                    column_names.append(c_name)
            return column_names

        with open(tables_path) as f:
            ddl = f.read()
        parse = sqlparse.parse(ddl)

        tables = {}
        for stmt in parse:

            # Get all the tokens except whitespaces
            tokens = [t for t in sqlparse.sql.TokenList(stmt.tokens) if t.ttype != sqlparse.tokens.Whitespace]
            is_create_stmt = False
            for i, token in enumerate(tokens):
                # Is it a create statements ?
                if token.match(sqlparse.tokens.DDL, 'CREATE'):
                    is_create_stmt = True
                    continue
                
                # If it was a create statement and the current token starts with "("
                if is_create_stmt and token.value.startswith("("):
                    # Get the table name by looking at the tokens in reverse order till you find
                    # a token with None type
                    table = get_table_name(tokens[:i]).strip('`')
                    tables[table] = get_column_names(token)
                    break

        return tables

    def _get_fasta_path(self, organism, force=False):
            
        assembly = self.organism_info[self.organism_info['sql_name'] == organism]['assembly'].values[0]
        # small hack for human assembly because name can contain patch number (like GRCh38.p13) but fasta files in Ensembl are without patch number
        if organism == 'homo_sapiens':
            assembly = assembly.split('.')[0]
        fasta_url = "ftp://ftp.ensembl.org/pub/release-{}/fasta/{}/dna/{}.{}.dna.toplevel.fa.gz".format(self.release, organism, organism.capitalize(), assembly)
        fasta_path = self.downloader.get(fasta_url, force=force)
        return fasta_path

    def _get_feature_class_path(self, organism, feature_class, force=False):

        short_name = self.organism_info[self.organism_info['sql_name'] == organism]['short_name'].values[0]
        feature_class_url = "ftp://ftp.ensembl.org/pub/release-{}/mysql/regulation_mart_{}/{}_{}__{}__main.txt.gz".format(self.release, self.release, short_name, feature_class, feature_class)
        print(feature_class_url)
        feature_class_path = self.downloader.get(feature_class_url, force=force)
        return feature_class_path

    def _get_column_names(self, organism, feature_class):

        return self.tables[self.organism_info[self.organism_info['sql_name'] == organism]['short_name'].values[0] + '_' + feature_class + '__' + feature_class + '__main']

    def _prepare_organisms_info(self):

        info_url = "ftp://ftp.ensembl.org/pub/release-{}/mysql/regulation_mart_{}/dataset_names.txt.gz".format(self.release, self.release)
        info_path = self.downloader.get(info_url)

        return pd.read_csv(info_path, sep='\t', header=None, usecols=[0, 5, 6, 7], names=['short_name', 'full_name', 'sql_name', 'assembly'])

    def _prepare_datasets_info(self):

        info_url = "ftp://ftp.ensembl.org/pub/release-{}/mysql/regulation_mart_{}/meta_conf__dataset__main.txt.gz".format(self.release, self.release)
        info_path = self.downloader.get(info_url)

        return pd.read_csv(info_path, sep='\t', header=None, usecols=[1, 2, 4], names=['dataset', 'display_name', 'type'])

    def _prepare_files_info(self):

        files_path = self.downloader.download('ftp://ftp.ensembl.org/pub/release-104/mysql/regulation_mart_104/')

        with open(files_path, 'r') as f:
            files = []
            for line in f.readlines():
                line = line.strip("\n")
                line = (list(filter(None, line.split(' '))))
                
                files.append([line[4], line[8]])

        files = pd.DataFrame(files, columns=['size', 'file'])
        files['size'] = files['size'].astype(int)
        files = files[files['size'] > 20]

        return files



    