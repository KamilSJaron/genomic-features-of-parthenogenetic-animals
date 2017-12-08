# Review of asexual genomes

This repository serves for supplementary analyses performed for review of asexual genomes.

The idea of review is to put together all full genomes are of asexual animals and compare patterns observed, other eukaryotes will be only discussed.
One of difficulties is to compare different genomics projects that are based on different inference methods and focus on different aspects.
This review aims to put them in a line taking this bias into account and unify analysis whenever it's easy enough.

List of performed analysis :

- estimates of heterozygosity between haplotypes using substitution model (atlas)

List of potential analysis :

- BUSCO
- Blobology
- classification of repeats
- analysis of palindromes

## What should be in this repository

- scripts for downloading, processing and analyzing asexual genomes
- a small table of analyzed asexual genomes, their code names and urls for downloading
- one big table -> an overview of all the asexual genomes
- other small summary tables of computationally intensive tasks
- the paper

## Sample labels

The labels of genomes are composed of **G**enus and **spe**cies name `Gspe`. Some of genomes were sequenced multiple times therefore every sample label is followed by index of the seuqenced individual. For for instance second sequencing of _Meloidogyne incognita_ wouble be `Minc2`.

***

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

### Execution of different cluster

`Snakefile` has no hardcoded any cluster-specific parameters. The resources should be accessed as `{resources.mem}` for memory in kilobytes, `{resources.tmp}` for needed local storage in megabytes and `{threads}` for number of used cores. The command used for cluster execution is stored in a bash wrapper `snakemake_clust.sh`. Modify this script as needed to work with syntax of your cluster. It uses environmental variable `USE_LOCAL` to access if computations should be performed on local disks of computational nodes or not (the job wrapper is `scripts/use_local.sh` and it might need to be adjusted to different cluster settings, now it's set for lsf).

The very last think to look at are program dependencies. These are actually hardcoded in job scripts as `module add UHTS/Analysis/sratoolkit/2.8.0;` for instance. These lines got to be deleted, but all the software got to be available.

- this should be extracted to loading script, so it get's separated from execution scripts

I know it sounds that there is a lot of things to do for reproduction of this analysis. However, I tried my best to combine reproducibility and good computational practice (like using local storage). I will wonk on that.

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
