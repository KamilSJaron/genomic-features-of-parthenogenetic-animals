### Tables in the manuscript

- [palindrome_table.tsv](palindrome_table.tsv) - **Table 1**
- **Supplementary Table 1**, overview of analyzed species. This table with typeset in LaTeX; the source can be found [here](LaTeX/SM_table_1_reproduction_modes.tex).
- [genome_table_infered_from_reads.tsv](genome_table_infered_from_reads.tsv) - **Supplementary Table 2**, genomic features derived from sequencing reads
- [assembly_table.tsv](assembly_table.tsv) - **Supplementary Table 3**, a table of genome assemblies; statistics, BUSCO and gene annotation stats
- [HGT_table.tsv](assembly_table.tsv) - **Supplementary Table 4**,

### The source tables

- [download_table.tsv](download_table.tsv) - the **master source table**; a table of species codes, species names, links to corresponding genome assemblies, IDs of sequencing reads and links to genome annotations.
- [genome_table.tsv](genome_table.tsv) - the **master result table**, a table of all estimated genomic properties; this table is used for nearly all figures and Supplementary tables 2 and 3.
- [Google doc table](https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing) - the **master literature table**; a table of all the collected published information about the genomes in this study

- [download_table_sexuals.tsv](download_table_sexuals.tsv) - analogical table to [download_table.tsv](download_table.tsv) containing sexual sister species when available; this data were used for HGT analysis and for genome comparison of the crayfish to its sister species.
- [triploid_heterozygosity.tsv](triploid_heterozygosity.tsv) - the estimates of the heterozygosity distribution in triploid genomes; used to generate Figure 3
- [tetraploid_heterozygosity.tsv](tetraploid_heterozygosity.tsv) - the estimates of the heterozygosity distribution in tetraploid genomes; used to generate Figure 3
- [gene_annotations.tsv](gene_annotations.tsv) - a summary table of tags in `gtf` and `gff` files referred in [download_table.tsv](download_table.tsv); used to extract protein sequences form genomes and annotation files for palindromes and HGT analyses.

### Other published data

- [Dcor1_meiosis_genes.tsv](Dcor1_meiosis_genes.tsv) - repertoire of meiosis genes in _D. coronatus_, taken from https://doi.org/10.1186/s12864-017-3860-x
- [nt_diversity.tsv](nt_diversity.tsv) - the nucleotide diversities extracted from the supplementary table 1 of  https://doi.org/10.1371/journal.pbio.1001388
