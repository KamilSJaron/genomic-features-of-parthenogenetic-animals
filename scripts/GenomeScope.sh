#/bin/bash

# 1, 2 input reads
# 3 output dir data/<sp>/genomescope

READS=$(printf -- '%s\n' "${@}" | grep reads-trimmed)
OUTDIR=${@:$#}

jellyfish count -C -m 21 -s 1000000000 -t 16 -o $OUTDIR/kmer_counts <(zcat $READS)

if [[ $(echo $OUTDIR/kmer_counts_* | wc -w) -eq 1 ]]
then
        mv $OUTDIR/kmer_counts_0 $OUTDIR/kmer_counts.jf
else
        jellyfish merge $OUTDIR/kmer_counts_* -o $OUTDIR/kmer_counts.jf
        rm $OUTDIR/kmer_counts_*
fi

jellyfish histo -t 16 $OUTDIR/kmer_counts.jf_0 > $OUTDIR/kmer.hist
genomescope.R $OUTDIR/kmer.hist 21 100 $OUTDIR 1000 verbose

rm $OUTDIR/kmer_counts.jf