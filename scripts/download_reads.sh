#!/bin/bash

SPECIES=$1
ACCESION=$2
R1_file=data/"$SPECIES"/raw_reads/"$ACCESION"_1.fastq.gz
R2_file=data/"$SPECIES"/raw_reads/"$ACCESION"_2.fastq.gz

# blocker if exist, don't dl again?

if [ -z $ACCESION ]; then
    echo "accession not found"
    exit 0
fi

mkdir -p data/"$SPECIES"/raw_reads

### trying EBI-SRA
# works with SRR610345,SRR610310,SRR801084,SRR7028347,SRR1191749,SRR1192095,SRR1192097,SRR8061554,SRR8061555,SRR8061556
URL_R1=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${ACCESION::6}/"$ACCESION"/"$ACCESION"_1.fastq.gz
URL_R2=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${ACCESION::6}/"$ACCESION"/"$ACCESION"_2.fastq.gz

wget $URL_R1 -O $R1_file
wget $URL_R2 -O $R2_file

if [[ -s $R1_file ]]; then
    echo "reads found at $URL_R1"
    exit 0;
fi

echo "failed to fetch from $URL_R1"

# works with ERR2135445, ERR2135453, ERR2135451
URL_R1=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${ACCESION::6}/00${ACCESION: -1}/"$ACCESION"/"$ACCESION"_1.fastq.gz
URL_R2=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${ACCESION::6}/00${ACCESION: -1}/"$ACCESION"/"$ACCESION"_2.fastq.gz

#### TRY DL
wget $URL_R1 -O $R1_file
wget $URL_R2 -O $R2_file

if [[ -s $R1_file ]]; then
    echo "reads found at $URL_R1"
    exit 0;
fi

echo "failed to fetch from $URL_R1"

### trying NCBI-SRA
# for vital it
module add UHTS/Analysis/sratoolkit/2.8.2.1;

URL=ftp://ftp.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/${ACCESION::3}/${ACCESION::6}/"$ACCESION"/"$ACCESION".sra

wget $URL -O data/"$SPECIES"/raw_reads/"$ACCESION".sra
fastq-dump data/"$SPECIES"/raw_reads/"$ACCESION".sra --outdir data/$SPECIES/raw_reads --split-files --gzip && rm data/"$SPECIES"/raw_reads/"$ACCESION".sra

if [[ -s $R1_file ]]; then
    echo "reads found at $URL"
    exit 0;
fi

echo "reads not found anywhere!"
exit 1;

# possibly faster
# getFASTQfile( in_acc = c('SRR5115148','SRR5115143'), sra_con, destDir = getwd(), srcType = 'ftp', ascpCMD = NULL )

#####
