#!/bin/bash
# TODO check if samtools view -h -q 20 works

module add UHTS/Aligner/bwa/0.7.15
module add UHTS/Analysis/samtools/1.3
# samblaster v0.1.24 installed locally

READS=$(printf -- '%s\n' "${@}" | grep reads-trimmed)

# run bwa mem
bwa mem -M -t 15 -R "@RG\tID:1\tSM:$1\tLB:$1" \
    data/$2/genome.fa.gz $READS | samtools view -bS - > temp.bam

samtools view -h temp.bam | samblaster -M | samtools view -h -q 20 | samtools sort -@10 -O bam - > ${@:$#}

rm temp.bam
