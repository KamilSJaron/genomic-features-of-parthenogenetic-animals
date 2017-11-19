#!/bin/bash

# get data from NCBI
# 2. arg - what to dl (data/Cbir/genome.fa.gz , data/Avag/reads.fq.gz)
COL=$(head -1 tables/download_table.tsv | tr "\t" "\n" | grep -n "genome" | cut -f 1 -d ':')

ROW=$(echo $1 | cut -f 2 -d '/')

URL=$(grep $ROW tables/download_table.tsv | cut -f $COL)

if [ -z $URL ]; then
    echo $URL is not an adress
    exit 0
fi

if [ ${URL##*.} -eq gz ]; then
    wget $URL -O $1
else
    wget $URL -O ${1%.gz}
    gzip ${1%.gz}
fi
