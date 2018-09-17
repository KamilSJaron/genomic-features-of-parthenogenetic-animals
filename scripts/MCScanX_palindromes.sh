#!/bin/bash

module add Development/java/1.8.0_172;
module add Blast/ncbi-blast/2.7.1+;

SP=$1
GENOME=$2
GFF=$3
MCScanX_DIR=$4
# PROTEINS=CF_000611835.1_CerBir1.0_proteins_cleaned.fa

mkdir -p $MCScanX_DIR

## TODO preparation - extract prot, run blast

##### RUN COLINEARITY ANALYSIS

MCScanX "$MCScanX_DIR"/"$SP"_prot

rm -r "$MCScanX_DIR"/"$SP"_prot.html
mv "$MCScanX_DIR"/"$SP"_prot* /scratch/beegfs/monthly/kjaron/review-of-asexual-genomes/data/$SP/MCScanX/

rm $GENOME_UNZIPED $GFF_UNZIPED $PROTEINS*