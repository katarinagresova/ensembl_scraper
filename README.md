# Ensemble scraper

> Create all functional elements datasets you ever wanted.

Ensemble scraper is command-line tool for accessing data from Ensemble and creating classification datasets from them. 

[Ensembl](https://www.ensembl.org/index.html) is a genome browser for vertebrate genomes that supports research in comparative genomics, evolution, sequence variation and transcriptional regulation. Ensembl annotate genes, computes multiple alignments, predicts regulatory function and collects disease data.

## Usage

You can install this package using pip and github repository as follows:

```bash
pip install git+https://github.com/katarinagresova/ensembl_scraper
```

## Local development

If you want to run experiments from this repository or contribute to the package, use following commands to clone the repository and install the package into virtual environment.

```bash
git clone git@github.com:katarinagresova/ensembl_scraper.git
cd ensembl_scraper

virtualenv venv --python=python3.8
source venv/bin/activate

pip install -e .
```

## Running tests

For testing, you need to install `pytest` and `pytest-cov` packages.

To run a specific test

```bash
    pytest -v ./tests/test_specific_file.py
```

To run all tests

```bash
    pytest -v tests/
```

To get a test coverage
```bash
    pytest --cov=ensembl_scraper/ tests/ 
```

## Features

 - downloading loci of functional elements for organisms of interest
 - converting data to specified format
 - preprocessing to remove low quality data
 - generating negative class for classification dataset

