# Review of asexual genomes

This repository serves for supplementary analyses performed for review of asexual genomes.

The idea of review is to put together all full genomes are of asexual animals and compare patterns observed, other eukaryotes will be only discussed.
One of difficulties is to compare different genomics projects that are based on different inference methods and focus on different aspects.
This review aims to put them in a line taking this bias into account and unify analysis whenever it's easy enough.

List of performed analysis :

- [GenomeScope](https://github.com/schatzlab/genomescope) `dev` - genome profiling from kmer spectra of sequencing reads - estimates of genome size, heterozygosity and repetitive content
- [KAT](https://github.com/TGAC/KAT) `2.4.1` - kmer spectra analysis in context of genome assembly - resolving how much are individual haplotypes collapsed in the assembly
- [MUMmer](https://github.com/mummer4/mummer/blob/master/MANUAL.md) `v4.0.0beta2` - genome self-alignment - evaluation of genome structure
- [dnaPipeTE](https://github.com/clemgoub/dnaPipeTE) `v1.2` - evaluation of repetitive content using sequencing reads
- [BUSCO](https://busco.ezlab.org/) `v3` - benchmarking using single copy orthologs - evaluation of unusual conserved gene content

## Sample labels

The labels of genomes are composed of **G**enus and **spe**cies name `Gspe`. Some of genomes were sequenced multiple times therefore every sample label is followed by index of the seuqenced individual. For for instance second sequencing of _Meloidogyne incognita_ wouble be `Minc2`.

***

## Development

**What should be in this repository:**

- scripts for downloading, processing and analyzing asexual genomes
- a small table of analyzed asexual genomes, their code names and urls for downloading
- one big table -> an overview of all the asexual genomes
- other small summary tables of computationally intensive tasks
- the paper

The analysis is automated using [snakemake](https://snakemake.readthedocs.io/en/stable/), tested with version `4.8.0`.
The scripts for analysis are combination of `bash`, `R` anf `python`.

## Vital-it execution

I wrote a [wrapper](snakemake_clust.sh) around snakemake so it understands our cluster. Then the execution is

```
./snakemake_clust.sh {target} {other_flags} {other_flags} ...
```

to run default with other flags you can run

```
./snakemake_clust.sh " " {other_flags} {other_flags} ...
```

### Execution of different cluster

`Snakefile` has no hardcoded any cluster-specific parameters. The resources should be accessed as `{resources.mem}` for memory in kilobytes, `{resources.tmp}` for needed local storage in megabytes and `{threads}` for number of used cores. The command used for cluster execution is stored in a bash wrapper `snakemake_clust.sh`. Modify this script as needed to work with syntax of your cluster. It uses environmental variable `USE_LOCAL` to access if computations should be performed on local disks of computational nodes or not (the job wrapper is `scripts/use_local.sh` and it might need to be adjusted to different cluster settings, now it's set for lsf).

The very last think to look at are program dependencies. These are actually hardcoded in job scripts as `module add UHTS/Analysis/sratoolkit/2.8.0;` for instance. These lines got to be deleted, but all the software got to be available.

- this should be extracted to loading script, so it get's separated from execution scripts

I know it sounds that there is a lot of things to do for reproduction of this analysis. However, I tried my best to combine reproducibility and good computational practice (like using local storage). I will wonk on that.

## Snakemake hints

apparently I can produce a graph of the workflow (it's really pretty) :

```
snakemake --forceall --dag | dot -Tpng > dag1.png
```

Snakemake version of `make clean` (removed all downloaded and coputed data) is

```
rm $(snakemake --summary | tail -n+2 | cut -f1)
```

Snakemake version of GNU make dry run (show commands, don't execute them) is

```
snakemake -p --quiet -n
```
