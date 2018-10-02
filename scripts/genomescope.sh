#/bin/bash

# 1 direcotry with kmers
# 2 output dir data/<sp>/genomescope

INDIR=$1
OUTDIR=$2

# name of histogram file should contain all library IDs and kmer size
HIST=$INDIR/*.hist
KMER=$(echo $HIST | cut -f 2 -d '_' | cut -b 2-3)

genomescope.R $HIST $KMER 100 $OUTDIR 10000 verbose
