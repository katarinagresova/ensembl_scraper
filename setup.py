import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

requires = [
    'bio'
    'biopython'
    'certifi'
    'charset-normalizer'
    'idna'
    'joblib'
    'numpy'
    'pandas'
    'plac'
    'pyfiglet'
    'python-dateutil'
    'pytz'
    'PyYAML'
    'request'
    'scikit-learn'
    'scipy'
    'six'
    'threadpoolctl'
    'tqdm'
    'twobitreader'
    'urllib3'
]

setuptools.setup(
    name="scraper",
    version="0.0.1",
    author="Katarina Gresova",
    author_email="gresova11@gmail.com",
    description="Create all functional elements datasets you ever wanted.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/katarinagresova/ensembl_scraper",
    packages=setuptools.find_packages(include=['scraper']),
    python_requires='>=3.6',
    install_requires=requires
)