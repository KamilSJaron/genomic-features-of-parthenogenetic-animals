module add Development/java/1.8.0_172;
module add Blast/ncbi-blast/2.7.1+;

PROTEINS=$1
FILTER_SCRIPT=$2
BLASTOUT=$3

##### GENERATE BLAST of ALL PROTEINS vs ALL PROTEINS
# create blast database
makeblastdb -in $PROTEINS -dbtype prot

# blast all proteins vs all proteins
blastp -query $PROTEINS -db $PROTEINS \
       -evalue 1e-10 -outfmt 6 -num_threads 16 | $FILTER_SCRIPT 5 > $BLASTOUT