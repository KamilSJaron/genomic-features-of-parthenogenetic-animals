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

GENOME_UNZIPED=${GENOME%.*}
GFF_UNZIPED=${GFF%.*}
PROTEINS=data/$SP/$(basename $GFF .gff3.gz)_proteins.fa

zcat $GENOME > $GENOME_UNZIPED
zcat $GFF > $GFF_UNZIPED

# extract proteins TODO test if <() works for it
gffread -g $GENOME_UNZIPED -y $PROTEINS $GFF_UNZIPED
sed -i'' -e 's/.$//' $PROTEINS

##### GENERATE BLAST of ALL PROTEINS vs ALL PROTEINS
# create blast database
makeblastdb -in $PROTEINS -dbtype prot

# blast all proteins vs all proteins
blastp -query $PROTEINS -db $PROTEINS \
    -out "$MCScanX_DIR"/"$SP"_prot.blast \
    -evalue 1e-10 -outfmt 6 -num_alignments 5 -num_threads 4

##### GENERATE ANNOTATION FILE UNDERSTANDABLE BY
grep -A 1 "ID=gene" $GFF_UNZIPED | grep "ID=rna" > "$GFF_UNZIPED"_subset
cut -f 1 "$GFF_UNZIPED"_subset > "$GFF_UNZIPED"_scf_names
cut -f 9 "$GFF_UNZIPED"_subset | cut -f 2 -d = | cut -f 1 -d \; > "$GFF_UNZIPED"_transcript_names
cut -f 4,5 "$GFF_UNZIPED"_subset > "$GFF_UNZIPED"_gene_positions

paste "$GFF_UNZIPED"_scf_names "$GFF_UNZIPED"_transcript_names "$GFF_UNZIPED"_gene_positions > "$MCScanX_DIR"/"$SP"_prot.gff
rm "$GFF_UNZIPED"_scf_names "$GFF_UNZIPED"_transcript_names "$GFF_UNZIPED"_gene_positions "$GFF_UNZIPED"_subset

##### RUN COLINEARITY ANALYSIS

MCScanX "$MCScanX_DIR"/"$SP"_prot

rm $GENOME_UNZIPED $GFF_UNZIPED $PROTEINS*