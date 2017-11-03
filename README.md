# Review of asexual genomes

This repository serves for supplementary analyses performed for review of asexual genomes.

The idea of review is to put together all full genomes are of asexual animals and compare patterns observed, other eukaryotes will be only discussed.
One of difficulties is to compare different genomics projects that are based on different inference methods and focus on different aspects.
This review aims to put them in a line taking this bias into account and unify analysis whenever it's easy enough.

Potential analysis :

- estimates of heterozygosity between haplotypes
- classification of repeats
- analysis of palindromes

## What should be in this repository

- scripts for downloading, processing and analyzing asexual genomes
- a small table of analyzed asexual genomes, their code names and urls for downloading
- one big table -> an overview of all the asexual genomes
- other small summary tables of computationally intensive tasks
- the paper

## Development

The analysis will be automated using `snakemake`, mainly because I would like to try it.
I will probably use combination of `python`, `bash` and `R`.
Details about all used software and versions should be recored here.

I want to pull data from NCBI and `snakemake` seems to have a module exactly for this. There is a function `snakemake.remote.NCBI` for pulling data from NCBI (here is its [documentation](http://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html#genbank-ncbi-entrez)).