# TODO

- Analyze theta outputs OR implement the script form timema project that sub select scaffolds that are at minimum of window size.
- wighted median for representation of genome heterozygosity
- Haploid / diploid / triploid reference will for sure affect heterozygosity calls. One think we can look at is correlation of coverage and estimated heterozygosity to figure out if high heterozygosity locations are separated regions collapsed during genome assembly.

# TODO ??

- BUSCO
- palindromes
- ? Blobology (check amounts of contaminations in assemblies)
- ? Base recalibration

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
