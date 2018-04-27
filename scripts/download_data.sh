#!/bin/bash

# get data from NCBI
# 1. arg - sample code
# 2. arg - what to dl (data/Cbir/$DATA.fa.gz , data/Avag/reads.fq.gz) annotation

ROW=$1
DATA=$2

COL=$(head -1 tables/download_table.tsv | tr "\t" "\n" | grep -n $DATA | cut -f 1 -d ':')
URL=$(grep $ROW tables/download_table.tsv | cut -f $COL)

if [[ $DATA == "genome" || $DATA == "proteins" ]]; then
    SUFIX=fa
fi

if [[ $DATA == "annotation" ]]; then
    SUFIX=gff3
fi

if [ -z $SUFIX ]; then
    echo $DATA is not a valid data to download
    exit 0
fi

if [ -z $URL ]; then
    echo $URL is not an adress, means that $ROW is not in the table
    exit 0
fi

if [[ ${URL##*.} == "gz" ]]; then
    echo wget $URL -O data/$ROW/$DATA.$SUFIX.gz
else
    echo wget $URL -O data/$ROW/$DATA.$SUFIX
    echo gzip data/$ROW/$DATA.$SUFIX
fi

touch data/$ROW/$DATA.$SUFIX.gz