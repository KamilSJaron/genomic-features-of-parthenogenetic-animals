#!/bin/bash

# get data from NCBI
# 2. arg - what to dl (data/Cbir/genome.fa.gz , data/Avag/reads.fq.gz)
COL=$(head -1 tables/download_table.tsv | tr "\t" "\n" | grep -n "reads" | cut -f 1 -d ':')
SPECIES=$1
ACCESION=$(grep $SPECIES tables/download_table.tsv | cut -f $COL)

if [ -z $ACCESION ]; then
    echo "accession not found"
    exit 0
fi

fastq-dump --accession $ACCESION --outdir data/$SPECIES --split-files --gzip
