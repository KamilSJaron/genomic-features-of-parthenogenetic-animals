#!/bin/bash

SP=$1

module add Development/java/1.8.0_172;
module add Blast/ncbi-blast/2.7.1+;

MCScanX_DIR=data/$SP/MCScanX
PROTEINS=$MCScanX_DIR/annotation_proteins.fa

##### GENERATE BLAST of ALL PROTEINS vs ALL PROTEINS
# create blast database
makeblastdb -in $PROTEINS -dbtype prot

# blast all proteins vs all proteins
blastp -query $PROTEINS -db $PROTEINS \
    -out "$MCScanX_DIR"/"$SP"_prot.blast \
    -evalue 1e-10 -outfmt 6 -num_alignments 5 -num_threads 32

MCScanX $MCScanX_DIR/"$SP"_prot

grep "# Number of" $MCScanX_DIR/"$SP"_prot.collinearity > $MCScanX_DIR/"$SP"_prot.collinearity_summary.txt;
python3 scripts/colinearity_parser.py $MCScanX_DIR/"$SP"_prot.collinearity >> $MCScanX_DIR/"$SP"_prot.collinearity_summary.txt;
grep "Align" $MCScanX_DIR/"$SP"_prot.collinearity | wc -l >> $MCScanX_DIR/"$SP"_prot.collinearity_summary.txt;

rm -r "$MCScanX_DIR"/"$SP"_prot.html

# when run manually
# cp "$MCScanX_DIR"/"$SP"_prot*.blast /scratch/beegfs/monthly/kjaron/review-of-asexual-genomes/data/$SP/MCScanX/
# mv "$MCScanX_DIR"/"$SP"_prot* /scratch/beegfs/monthly/kjaron/review-of-asexual-genomes/data/$SP/MCScanX/
#
# rm -r data/$SP

# generating summaries aposteriory
# for file in data/*/MCScanX/*_prot.collinearity; do
#     grep "# Number of" $file > "$file"_summary.txt;
#     python3 scripts/colinearity_parser.py $file >> "$file"_summary.txt;
#     grep "Align" $file | wc -l >> "$file"_summary.txt;
# done
