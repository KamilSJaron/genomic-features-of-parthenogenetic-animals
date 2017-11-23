#!/bin/bash

module add UHTS/Aligner/bwa/0.7.13
module add UHTS/Analysis/samtools/1.3

# run bwa mem
bwa mem -M -t 15 -R "@RG\tID:1\tSM:$1\tLB:$1" \
    data/$2/genome.fa.gz \
    data/$1/reads_R1.fq.gz \
    data/$1/reads_R2.fq.gz \
    | samtools sort -O bam - > ${@:$#}
