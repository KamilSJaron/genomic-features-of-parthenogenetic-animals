#!/bin/bash

module add UHTS/Aligner/bwa/0.7.13
module add UHTS/Analysis/samtools/1.3
# samblaster v0.1.24 installed locally

# run bwa mem
bwa mem -M -t 15 -R "@RG\tID:1\tSM:$1\tLB:$1" \
    data/$2/genome.fa.gz \
    data/$1/reads-trimmed-pair1.fastq.gz \
    data/$1/reads-trimmed-pair2.fastq.gz \
    | samblaster -M | samtools sort -O bam - > ${@:$#}
