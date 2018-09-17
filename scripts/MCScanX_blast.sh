#!/bin/bash

module add Blast/ncbi-blast/2.7.1+;

MCScanX_DIR=data/$SP/MCScanX
PROTEINS=$MCScanX_DIR/annotation_proteins.fa

##### GENERATE BLAST of ALL PROTEINS vs ALL PROTEINS
# create blast database
makeblastdb -in $PROTEINS -dbtype prot

# blast all proteins vs all proteins
blastp -query $PROTEINS -db $PROTEINS \
    -out "$MCScanX_DIR"/"$SP"_prot.blast \
    -evalue 1e-10 -outfmt 6 -num_alignments 5 -num_threads 16

mv "$MCScanX_DIR"/"$SP"_prot* /scratch/beegfs/monthly/kjaron/review-of-asexual-genomes/data/$SP/MCScanX/

rm $PROTEINS*