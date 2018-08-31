#!/bin/bash

module add UHTS/Assembler/cufflinks/2.2.1;
module add Development/java/1.8.0_172;
module add Blast/ncbi-blast/2.7.1+;

SP=$1
GENOME=$2
GFF=$3
MCScanX_DIR=$4
# PROTEINS=CF_000611835.1_CerBir1.0_proteins_cleaned.fa

mkdir -p $MCScanX_DIR

## TODO preparation

##### GENERATE BLAST of ALL PROTEINS vs ALL PROTEINS
# create blast database
makeblastdb -in $PROTEINS -dbtype prot

# blast all proteins vs all proteins
blastp -query $PROTEINS -db $PROTEINS \
    -out "$MCScanX_DIR"/"$SP"_prot.blast \
    -evalue 1e-10 -outfmt 6 -num_alignments 5 -num_threads 4

##### RUN COLINEARITY ANALYSIS

MCScanX "$MCScanX_DIR"/"$SP"_prot

rm $GENOME_UNZIPED $GFF_UNZIPED $PROTEINS*