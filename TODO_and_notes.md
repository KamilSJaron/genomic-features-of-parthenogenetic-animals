# TODO

- Analyze theta outputs OR implement the script form timema project that sub select scaffolds that are at minimum of window size.
- wighted median for representation of genome heterozygosity
- adjust / fix / outsmart *reference ploidy problem* -> look at the subsection bellow

### reference ploidy problem

Haploid / diploid reference will for sure affect heterozygosity calls.

#### Coverage solution

One think we can look at is correlation of coverage and estimated heterozygosity to figure out if high heterozygosity locations are separated regions collapsed during genome assembly.

#### Sexual reference solution

reference assembly is combination of -> decent genomic sequences + contaminants + non-collapsed fragments + tiny contigs. We would like to keep only "decent genomic sequences" for downstream analysis. This idea is based on the same divergence of both asexual haplotypes to sexual sister species.

1. Map asexual genome to sexual genome and Identify contigs that map -> trustworthy contigs
2. Map asexual reads to asexual genome, subselect those reads that map to "trustworthy contigs" and call them "trustworthy reads"
3. Map "trustworthy reads" to sexual assembly
4. Call theta

Also...

- sexual contaminant contigs will have no read coverage of "trustworthy reads"
- asexual that represent non-collapsed haplotypes will map to the same location on sexual genome
- well covered one-to-one mappings represent reliable part of dataset, that can be further supported by coverage analysis
- keep tiny contigs in reference for missmapping errors

# TODO ??

- BUSCO
- palindromes -> Jonathan agreed to gather all the annotations and proteins
- ? Blobology (check amounts of contaminations in assemblies)
- *recalibration* -> look at the subsection bellow
- *triploid theta* -> look at the subsection bellow

### Softmask references

This should be problem of masking tandem replications and stuff like that. Stuff that have been assembled separately.

Hence for my genomes I should

- verify that they are softmasked

OR

- Quick & Dirty window modeler?

OR

- use Jens' TE prediction with Repeat Masker

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

# Notes

- A springtail Holacanthella duospinosa is sequenced. Is it the closed sequenced relative of the asexual springtail genome?

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

Yhe results are comparable to samblaster.  0.486266. I can leave it like this.
