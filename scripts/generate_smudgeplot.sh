#!/bin/bash

JELLYFISH_DIR=$1
OUTDIR=$2

JF_FILE=$JELLYFISH_DIR/*jf
HISTOGRAM=$JELLYFISH_DIR/*hist
KMER=$(echo $HISTOGRAM | cut -f 2 -d '_' | sed 's/^k//' | sed 's/.hist//')

mkdir $OUTDIR

SP=$(echo $OUTDIR | cut -f 2 -d \/)

if [[ $SP == "Lcla1" || $SP == "Tpre1" || $SP == "Aruf1" || $SP == "Rvar1" ]]; then
    HOMOZYG="--homozygous"
fi

L=$(kmer_cov_cutoff.R $HISTOGRAM L)
U=$(kmer_cov_cutoff.R $HISTOGRAM U)

jellyfish dump -c -L $L -U $U $JF_FILE | hetkmers.py -k $KMER -t 8 -o $OUTDIR/kmer_pairs

smudgeplot.R -i $OUTDIR/kmer_pairs_coverages_2.tsv -o $OUTDIR/$SP -t $SP -L $L -k $KMER $HOMOZYG
