#/bin/bash

# 1 input read directory
# 2 output dir data/<sp>/genomescope

OUTDIR=$2
mkdir -p $OUTDIR

#zcat $1 $2 > $OUTDIR/trimmed_reads.fasta

jellyfish count -C -m 21 -s 1000000000 -t 16 -o $OUTDIR/kmer_counts <(zcat "$1"/*)

if [[ $(echo $OUTDIR/kmer_counts_* | wc -w) -eq 1 ]]
then
        mv $OUTDIR/kmer_counts_0 $OUTDIR/kmer_counts.jf
else
        jellyfish merge $OUTDIR/kmer_counts_* -o $OUTDIR/kmer_counts.jf
fi

jellyfish histo -t 16 $OUTDIR/kmer_counts.jf > $OUTDIR/kmer.hist
genomescope.R $OUTDIR/kmer.hist 21 100 $OUTDIR 1000 verbose

#rm $OUTDIR/kmer_counts.jf
