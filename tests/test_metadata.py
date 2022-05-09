from ensembl_scraper.regulatory import Metadata
import pytest
from pathlib import Path

def test_get_organisms():

    organisms_in_releases = {
        '103': ['homo_sapiens'],
        '106': ['mus_musculus', 'homo_sapiens']
    }

    for release, organisms in organisms_in_releases.items():
        metadata = Metadata(release)
        assert set(metadata.get_organisms()) == set(organisms)

def test_get_full_name_exception():

    metadata = Metadata(106)

    with pytest.raises(Exception):
        metadata.get_full_name(set('aotus_nancymaae'))

def test_get_full_name():

    metadata = Metadata(106)

    full_names = {
        'aotus_nancymaae': 'Ma\'s night monkey',
        'bos_taurus': 'Cow',
        'caenorhabditis_elegans': 'Caenorhabditis elegans',
        'cyprinodon_variegatus': 'Sheepshead minnow'
    }

    for sql_name, full_name in full_names.items():
        assert metadata.get_full_name(sql_name) == full_name

def test_get_full_name_list():

    metadata = Metadata(106)

    assert set(metadata.get_full_name(['danio_rerio', 'homo_sapiens', 'ictalurus_punctatus', 'macaca_nemestrina'])) == set(['Zebrafish', 'Human', 'Channel catfish', 'Pig-tailed macaque'])

def test_get_feature_classes():

    feature_classes_for_organisms = {
        'homo_sapiens': ['external_feature', 'regulatory_feature', 'mirna_target_feature', 'peak'],
        'mus_musculus': ['external_feature', 'regulatory_feature', 'mirna_target_feature', 'peak'],
        'drosophila_melanogaster': [],
    }

    metadata = Metadata(106)

    for organism, feature_classes in feature_classes_for_organisms.items():
        assert set(metadata.get_feature_classes(organism)) == set(feature_classes)

def test_get_fasta_path():

    home = Path.home()
    fasta_paths_for_organisms = {
        'homo_sapiens': Path(home, '.ensembl_scraper/106/data/Homo_sapiens.GRCh38.dna.toplevel.fa')
    }

    metadata = Metadata(106)

    for organism, fasta_path in fasta_paths_for_organisms.items():
        assert metadata.get_fasta_path(organism).samefile(fasta_path)

def test_get_column_names():

    column_names_for_organisms_and_features = {
        ('mus_musculus', 'external_feature'): [
            'feature_type_class_1049', 
            'so_accession_1024', 
            'seq_region_strand_1049', 
            'seq_region_start_1049', 
            'display_label_1049', 
            'seq_region_end_1049', 
            'fs_display_label_1049', 
            'so_term_1024', 
            'external_feature_id_1021_key', 
            'feature_type_description_1049', 
            'seq_region_name_1049'
        ],
        ('mus_musculus', 'regulatory_feature'): [
            'feature_type_description_1051',
            'bound_seq_region_start',
            'feature_type_name_1051',
            'bound_seq_region_end',
            'so_accession_1024',
            'seq_region_name_1051',
            'stable_id_1051',
            'so_term_1024',
            'regulatory_feature_id_1036_key',
            'seq_region_start_1051',
            'seq_region_end_1051',
            'seq_region_strand_1051'
        ]
    }

    metadata = Metadata(106)

    for (organism, feature), column_names in column_names_for_organisms_and_features.items():
        assert metadata.get_column_names(organism, feature) == column_names