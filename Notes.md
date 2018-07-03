# Notes

## Reference ploidy and heterozygosity

We use kmer-profiles to estimate genome sizes, and heterozygosity estimates (genomescope). The model already works for diploids, the polyploid solution is currently developed and done by collaborators.

An alternative for heterozygosity estimate is bellow

## heterozygosity by atlas

- [atlas](https://bitbucket.org/phaentu/atlas/wiki/Home) - estimates of heterozygosity between haplotypes using substitution model

The elegant way how to estimate heterozygosity is maximum likelihood, implemented for instance in [atlas package](https://bitbucket.org/phaentu/atlas).
To estimate heterozygosity one needs to map reads to reference genome and then estimate heterozygosity of the mapped individual.
The method is not using sequence of the reference, reference is needed only for the alignment of the reads.
The likelihood is calculated using quality scores of individual bases that were sequenced. The model is currently diploid.

We did not used this method in the end.
The reason is because there are couple of potential problem if the reference is not of high quality:

 - some of the paralogs will be collapsed into one sequence (creating artificially heterozygous region)
 - some of the alternative haplotypes get separately assembled (creating two artificially homozygous regions)
 - some of the genomes in this review are not diploid
 - Illumina is overestimating quality scores a lot and it's not streightforward to recalibrate them without other knowledge of genomes (location of conserved loci, known haploid regions like sex chromosomes or mtDNA)

 However, we spent some time thinking about this problems therefore here are some thoughts we had:

### Softmask references

This should be problem of masking tandem replications and stuff like that. Stuff that have been assembled separately.

Hence for my genomes I should

- verify that they are softmasked

OR

- Quick & Dirty window modeler?

OR

- use Jens' TE prediction with Repeat Masker

### Coverage solution

One think we can look at is correlation of coverage and estimated heterozygosity to figure out if high heterozygosity locations are separated regions collapsed during genome assembly.

### Recalibration

The quality of SNP calls will be always reflecting accuracy of errors. Therefore it would be worth to try to recalibrate if possible. If not, we need to rely on the sufficient amount of data we have and just assume that it works.

The possible "haploid" locations could be

- conserved genes (CEGMA, BUSCO)
- mtDNA

Other possibility would be Lancer, the other non parametric recalibration tool.

### Triploid model

We can estimate two thetas instead of one (one corresponding to divergence between two more closely related haplotypes and the second to divergence of these and the more distantly related haplotype). Then on every single position we know exact probabilities of two read comparison being an comparison of theta for the first case (2/9) for the second case (4/9) and coming from the same haplotype (3/9).

- This model require trustworthy triploid mapping
- should not be constructed if problems above won't be resolved

### Marking duplicates

I tried `samblaster` on `Lcla1`, resulting in brutal :

```
samblaster: Marked 6122352 of 12722150 (48.12%) read ids as duplicates using 142804k memory in 49.693S CPU seconds and 22M13S(1333S) wall time.
```

therefore I will try more traditional approach using Picard tools.

```
module add UHTS/Analysis/picard-tools/2.2.1

picard-tools MarkDuplicates INPUT=map_to_Lcla1.bam OUTPUT=map_to_Lcla1_picard_marked.bam METRICS_FILE=picard_metrics.txt
```

The results are comparable to samblaster.  0.486266. I can leave it like this.

# Snakemake notes

I had some decisions and notes I gathered along the way

### Pulling sequencing data using Snakemake

I wrote my own manual puller. However there is an option to use directly `snakemake`. Maybe I should use the native functionality. The problem was that our cluster had an older version without the functionality when I started this project, but it's updated now.

I want to pull data from NCBI and `snakemake` seems to have a module exactly for this. There is a function `snakemake.remote.NCBI` for pulling data from NCBI (here is its [documentation](http://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html#genbank-ncbi-entrez)). One would probably like to specify `--keep-remote` if this option will be used.

### other flags to consider

- `--cluster-status 'bjobs'` : allow snakemake to look at the status of jobs; ~this is not working on version of snakemake on Vital-it~ it was updated, but for some reason cluster status is not allowed even in version 11
- `--jobscript cluster_wrapper.sh`