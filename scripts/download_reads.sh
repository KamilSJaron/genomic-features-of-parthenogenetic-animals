#!/bin/bash

# for vital it
module add UHTS/Analysis/sratoolkit/2.8.2.1;
module add Utility/aspera_connect/3.7.4.147727;
OPENSSH=/software/Utility/aspera_connect/3.7.4.147727/etc/asperaweb_id_dsa.openssh

# get data from NCBI
COL=$(head -1 tables/download_table.tsv | tr "\t" "\n" | grep -n "reads" | cut -f 1 -d ':')
SPECIES=$1
ACCESION=$(grep -w $SPECIES tables/download_table.tsv | cut -f $COL)

if [ -z $ACCESION ]; then
    echo "accession not found"
    exit 0
fi

mkdir -p data/"$SPECIES"

URL=anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/${ACCESION::3}/${ACCESION::6}/"$ACCESION"/"$ACCESION".sra

ascp -QT -i $OPENSSH $URL data/"$SPECIES"

fastq-dump data/"$SPECIES"/"$ACCESION".sra --outdir data/$SPECIES --split-files --gzip

rm data/"$SPECIES"/"$ACCESION".sra
mv data/$SPECIES/"$ACCESION"_1.fastq.gz data/$SPECIES/reads_R1.fq.gz
mv data/$SPECIES/"$ACCESION"_2.fastq.gz data/$SPECIES/reads_R2.fq.gz
