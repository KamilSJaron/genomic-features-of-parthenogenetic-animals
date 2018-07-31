#/bin/bash
# wild script - used to calculate GC content of kmer peaks

# 1 sp
# 2 jellyfish file
# 3 threshold
# 4 threshold
# 5 (threshold) ...

SP=$1
# SP=Aric1
JK_FILE=$2
# JK_FILE=data/Aric1/genomescope/kmer_counts.jf

shift 2;

LOWER_COUNT=0
OUTPUT=data/$SP/kmer_spectra_GC/summary.tsv
OUTDIR=$(dirname $OUTPUT)

mkdir -p $OUTDIR
calc(){ awk "BEGIN { print "$*" }"; } # function for division
#
# THRESHOLDS=(20 55 90 220 400)
# for threshold in ${THRESHOLDS[*]}; do
for threshold in $@; do
    UPPER_COUNT=$threshold
    jellyfish dump --lower-count=$LOWER_COUNT --upper-count=$UPPER_COUNT $JK_FILE > $OUTDIR/extracted_kmers_"$LOWER_COUNT"-"$UPPER_COUNT".fa
    # do magic $OUTDIR/extracted_kmers.fa >> $OUTPUT
    GC=$(grep -v ">" $OUTDIR/extracted_kmers_"$LOWER_COUNT"-"$UPPER_COUNT".fa | tr -cd GC | wc -c)
    TOTAL=$(grep -v ">" $OUTDIR/extracted_kmers_"$LOWER_COUNT"-"$UPPER_COUNT".fa | tr -d "\n" | wc -c)
    echo $SP $LOWER_COUNT $UPPER_COUNT $(calc '(100*'"$GC"')/'"$TOTAL") $TOTAL >> $OUTPUT
    LOWER_COUNT=$UPPER_COUNT
done;

rm $OUTDIR/extracted_kmers_*