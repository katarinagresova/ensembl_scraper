# Ensemble scraper

> Create all functional elements datasets you ever wanted.

Ensemble scraper is command-line tool for accessing data from Ensemble and creating classification datasets from them. 

[Ensembl](https://www.ensembl.org/index.html) is a genome browser for vertebrate genomes that supports research in comparative genomics, evolution, sequence variation and transcriptional regulation. Ensembl annotate genes, computes multiple alignments, predicts regulatory function and collects disease data.

## Instalation
```
pip install git+https://github.com/katarinagresova/ensembl_scraper.git
wget https://raw.githubusercontent.com/katarinagresova/ensembl_scraper/main/requirements.txt
pip install -r requirements.txt
```


## Features

 - downloading loci of functional elements for organisms of interest
 - converting data to specified format
 - preprocessing to remove low quality data
 - generating negative class for classification dataset

