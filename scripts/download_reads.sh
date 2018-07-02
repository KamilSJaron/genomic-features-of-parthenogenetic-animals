#!/bin/bash

# for vital it
module add UHTS/Analysis/sratoolkit/2.8.2.1;
module add Utility/aspera_connect/3.7.4.147727;
OPENSSH=/software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh

# get data from NCBI
SPECIES=$1
ACCESION=$2

if [ -z $ACCESION ]; then
    echo "accession not found"
    exit 0
fi

mkdir -p data/"$SPECIES"/raw_reads

URL=anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/${ACCESION::3}/${ACCESION::6}/"$ACCESION"/"$ACCESION".sra

ascp -QT -i $OPENSSH $URL data/"$SPECIES"/raw_reads
fastq-dump data/"$SPECIES"/raw_reads/"$ACCESION".sra --outdir data/$SPECIES/raw_reads --split-files --gzip

rm data/"$SPECIES"/raw_reads/"$ACCESION".sra

# possibly faster
# getFASTQfile( in_acc = c('SRR5115148','SRR5115143'), sra_con, destDir = getwd(), srcType = 'ftp', ascpCMD = NULL )