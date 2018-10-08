#/bin/bash

# 1 direcotry with kmers
# 2 directory with smudgeplot
# 3 output dir data/<sp>/genomescope

KMERS_HIST=$1/*.hist
SMUDEGEPLOT_SUM=$2/*_verbose_summary.txt
READS=$3/*.fastq.gz
OUTDIR=$4

# name of histogram file should contain all library IDs and kmer size
KMER=$(echo $HIST | cut -f 2 -d '_' | cut -b 2-3)
READ_LEN=$(zcat $READS | head -2 | tail -1 | wc -c)
PLOIDY=$(grep "Estimated ploidy" $SMUDEGEPLOT_SUM | cut -f 2)     # smudgeplot, tried more p when not clear
LAMBDA=$(grep "peak 1n coverage estimate" $SMUDEGEPLOT_SUM | cut -f 2)

SP=$(echo $KMERS_HIST | cut -f 2 -d /)
mkdir -p data/$SP/genomescope_v2

genomescope.R data/$SP/jellyfish/*.hist $KMER $READ_LEN $PLOIDY data/$SP/genomescope_v2 $LAMBDA -1 verbose &> data/$SP/genomescope_v2/log
