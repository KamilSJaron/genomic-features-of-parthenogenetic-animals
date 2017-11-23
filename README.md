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

## Sample labels

The labels of genomes are composed of genus and species name `Gspe`. Some of genomes were sequenced multiple times therefore every sample label is followed by index of the seuqenced individual. For for instance second sequencing of _Meloidogyne incognita_ wouble be `Minc2`.


## Development

The analysis will be automated using `snakemake`, mainly because I would like to try it.
I will probably use combination of `python`, `bash` and `R`.
Details about all used software and versions should be recored here.

I want to pull data from NCBI and `snakemake` seems to have a module exactly for this. There is a function `snakemake.remote.NCBI` for pulling data from NCBI (here is its [documentation](http://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html#genbank-ncbi-entrez)).

## Vital-it execution

I will try to use already installed version. However it's not the most recent one.

```
module add Utility/snakemake/3.5.4
```

seems to be working. However, one have to pay attention to - be in a directory that is accessible from everywhere, no scratch local!

cluster vital-it command :

```
./snakemake_clust.sh {target} {other_flags} {other_flags} ...
```

to run default with other flags you can run

```
./snakemake_clust.sh " " {other_flags} {other_flags} ...
```

other flags to consider / test :

- `--cluster-status 'bjobs'` : allow snakemake to look at the status of jobs; this is not working on version of snakemake on Vital-it.
- `--jobscript cluster_wrapper.sh`
- `--keep-remote`

## other notes

apparently I can produce a graph of the workflow :

```
snakemake --forceall --dag | dot -Tpng > dag1.png
```

Snakemake version of `make clean` is

```
rm $(snakemake --summary | tail -n+2 | cut -f1)
```

Snakemake version of GNU make dry run is

```
snakemake -p --quiet -n
```
