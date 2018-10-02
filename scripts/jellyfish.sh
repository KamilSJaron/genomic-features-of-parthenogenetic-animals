#/bin/bash

# 1 input read directory
# 2 dl table
# 3 output dir data/<sp>/genomescope

DL_TABLE=$2
OUTDIR=$3

# leaving myself open doors for ther kmer values using environmental variable
if [[ ! -v KMER ]]; then
        KMER=21
fi

# name of histogram file should contain all library IDs and kmer size
HIST=$OUTDIR/$(grep "$sp" "$DL_TABLE" | cut -f 4)_k"$KMER".hist

mkdir -p $OUTDIR

#zcat $1 $2 > $OUTDIR/trimmed_reads.fasta

jellyfish count -C -m $KMER -s 1000000000 -t 16 -o $OUTDIR/kmer_counts <(zcat "$1"/*.fastq.gz)

# this construct makes sure that all the kmer files will end up merged in one file kmer_counts.jf
if [[ $(echo $OUTDIR/kmer_counts_* | wc -w) -eq 1 ]]
then
        mv $OUTDIR/kmer_counts_0 $OUTDIR/kmer_counts.jf
else
        jellyfish merge $OUTDIR/kmer_counts_* -o $OUTDIR/kmer_counts.jf
        rm data/*/genomescope/kmer_counts_*
fi