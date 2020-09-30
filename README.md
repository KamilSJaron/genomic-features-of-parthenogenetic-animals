# Genomic features of parthenogenetic animals

This repository is archived with [![DOI](https://zenodo.org/badge/106611468.svg)](https://zenodo.org/badge/latestdoi/106611468) and it contrains the analyses performed for a peer-reviewed aricle [Genomic features of parthenogenetic animals](https://doi.org/10.1093/jhered/esaa031).

The idea is to review all genomes of asexual animals and compare patterns observed.
One of difficulties is to compare different genomics projects that are based on different inference methods and focus on different aspects. Therefore we estimate all the genomic features possible using unified methodology.

### Regenerating figures

The main figures are be plotted by following R scripts

```
Rscript scripts/plot_figure_1_questions.R --tricolor
# generates figures/fig1_genomic_studies.pdf
Rscript scripts/plot_figure_2_heterozygosity.R --split_axis --homoeolog --rm_boxes
# generates figures/fig2_heterozygosity_split_axis_homoeolog.pdf
Rscript scripts/plot_figure_3_heterozygosity_structure.R
# generates figures/fig3_heterozygosity_of_tetraploids.pdf
Rscript scripts/plot_figure_4_TEs.R
# generates figures/fig4_TEs.pdf
```

The supplementary figures

```
Rscript scripts/plot_figure_1_questions.R --refs --both --tricolor
# supp figure 3 is a derivative of figure 1
# generates figures/SM_Figure_3_genomic_studies.pdf
Rscript scripts/plot_figure_S6_expected_heterozygosity_structure.R
Rscript scripts/plot_figure_S8_TEs_vs_mode_and_origin.R
# figures/SM_Figure_8_TEs
```

### Supplementary tables

[Supplementary Table 1](LaTeX/SM_table_1_reproduction_modes.pdf): Overview of analysed species. This information was collected directly from the cited literature. References include information regarding cellular mode of reproduction, origin and/or the age of parthenogenesis.


[Supplementary Table 2](tables/genome_table_infered_from_reads.tsv): Genomic features calculated from raw data. We used unified methods to estimate basic genomic properties directly from sequencing reads. Ploidy was estimated using smudgeplot for all species but A. vaga (see section Heterozygosity structure in polyploids for details). Genome size, heterozygosity and repeats were estimated using GenomeScope. Repeats denote the fraction of the genome occurring in more than one copy. The classified repeats, TEs and other types of classified repeats, were estimated using DnaPipeTE.


[Supplementary Table 3](tables/assembly_table.tsv): genome assemblies: size, number of scaffolds, N50, BUSCO, number of annotated genes. Statistics were calculated from the published genome assemblies and genome annotations shared by authors. BUSCO genes were searched using the metazoan database for all the non-nematode species. Nematodes are notoriously known for the high turnover of genes and we therefore used nematode specific BUSCO genes. The number of annotated genes were calculated as the number of lines in the annotation with the tag “gene”. The number of genes was extracted using the tag “mRNA” since the keyword “gene” was not in the annotation file of Diploscapter coronatus.

Supplementary Table 4: Horizontal gene transfer analysis.
HGT candidate genes identified from comparisons to [UniRef90](tables/JOH-2020-024.S4Table.HGT_sheet1_uniref.tsv) and [UniProtKB/Swissprot](tables/JOH-2020-024.S4Table.HGT_sheet2_uniprot.tsv) databases. Column TaxID is NCBI TaxID for focal species; Num.genes is the Number of protein-coding genes in annotation; HGTc is Horizontal gene transfer candidates (i.e. putative foreign gene); Columns E-H: Phylum/Class/Order/Family indicates the taxonomic level at which hits to the focal animal's lineage were discounted; Columns I-L: HGTc expressed as a percentage of total CDS encoded in focal genome (Column D). The supplementary table in the paper is an excel sheet.

### List of performed analysis

- [GenomeScope](https://github.com/tbenavi1/genomescope) `v2 dev` - genome profiling from kmer spectra of sequencing reads. Estimate of genome size, heterozygosity and repetitive content.
- [smudgeplot](https://github.com/tbenavi1/smudgeplot) `v0.1.3` - Estimation of ploidy and heterozygosity structure, helps in interpretation of GenomeScope kmer spectra.
- [MUMmer](https://github.com/mummer4/mummer/blob/master/MANUAL.md) `v4.0.0beta2` - Genome self-alignment. Evaluation of genome structure.
- [dnaPipeTE](https://github.com/clemgoub/dnaPipeTE) `v1.2` - Evaluation of repetitive content using sequencing reads.
- [BUSCO](https://busco.ezlab.org/) `v3` - Benchmarking using single copy orthologs. Evaluation of unusual conserved gene content.
- [MCScanX](http://chibba.pgml.uga.edu/mcscan2/) `untagged version` released 28.3.2013 - Collinearity and palindrome analysis. Evaluation of abundance of palindromes in asexual genomes.
- [hgt](https://github.com/reubwn/hgt) the latest commit `917ce74` - Analysis of horizontal gene transfer.
- [KAT](https://github.com/TGAC/KAT) `2.4.1` - Kmer spectra analysis in context of genome assembly. Resolving how much are individual haplotypes collapsed in the assembly. (performed, but the analysis was not used for any interpretation in the manuscript)

## Sample labels

The labels of genomes are composed of **G**enus and **spe**cies name `Gspe` and an index, which serves only as a distinction of the different sequencing projects. The full list of genomes considered in this study is in the table [tables/download_table.tsv](tables/download_table.tsv).

**What is in this repository:**

- scripts for downloading, processing and analyzing asexual genomes
- a table of analyzed asexual genomes, their code names and urls for downloading
- summary tables of the results

The analysis is automated using [snakemake](https://snakemake.readthedocs.io/en/stable/), tested with version `4.8.0`.
The scripts for analysis are combination of `bash`, `R` anf `python`.

## Cluster execution

I wrote a [wrapper](snakemake_clust.sh) around snakemake so it understands our cluster. Then the execution is

```
./snakemake_clust.sh {target} {other_flags} {other_flags} ...
```

to run default with other flags you can run

```
./snakemake_clust.sh " " {other_flags} {other_flags} ...
```

### Execution on a different cluster

`Snakefile` has no hardcoded any cluster-specific parameters. The resources should be accessed as `{resources.mem}` for memory in kilobytes, `{resources.tmp}` for needed local storage in megabytes and `{threads}` for number of used cores. The command used for cluster execution is stored in a bash wrapper `snakemake_clust.sh`. Modify this script as needed to work with syntax of your cluster. It uses environmental variable `USE_LOCAL` to access if computations should be performed on local disks of computational nodes or not (the job wrapper is `scripts/use_local.sh` and it might need to be adjusted to different cluster settings, now it's set for lsf).

The very last think to look at are program dependencies. These are actually hardcoded in job scripts as `module add UHTS/Analysis/sratoolkit/2.8.0;` for instance. These lines got to be deleted, but all the software got to be available.

- this should be extracted to loading script, so it get's separated from execution scripts

I know it sounds that there is a lot of things to do for reproduction of this analysis. However, I tried my best to combine reproducibility and good computational practice (like using local storage). I will wonk on that.

## Snakemake tips

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
