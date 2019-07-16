#!/bin/bash

cp -r /scratch/beegfs/monthly/kjaron/genomic-features-of-asexual-animals/data/Pvir1/raw_reads .
cp -r /scratch/beegfs/monthly/kjaron/genomic-features-of-asexual-animals/data/Pvir1/trimmed_reads .

ls raw_reads/*.fastq.gz > FILES
ls trimmed_reads/*.fastq.gz > TRIMMED_FILES
HIST_SUFFIX=_SRR5115143,SRR5115144,SRR5115145,SRR5115146,SRR5115147,SRR5115148_k17.hist

mkdir -p tmp
# 325000001 is max kmer
kmc -k17 -t32 -m64 -ci1 -cs325000005 @FILES data/Pvir1/pvir1_kmc_kmer_counts tmp
kmc_tools transform pvir1_kmc_kmer_counts histogram kmc_$HIST_SUFFIX -cx325000005

kmc -k17 -t60 -m64 -ci1 -cs500000000 @TRIMMED_FILES pvir1_kmc_trimmed_kmer_counts tmp
kmc_tools transform pvir1_kmc_trimmed_kmer_counts histogram temp.hist -cx500000000
grep -v "[[:space:]]0" temp.hist > kmc_trimmed_$HIST_SUFFIX # cleaning up the rows with 0 coverage

jellyfish count -C -m 17 -s 1000000000 -t 32 -o pvir_jelly_kmer_counts <(zcat raw_reads/*.fastq.gz)
jellyfish merge pvir_jelly_kmer_counts_* -o pvir_jelly_kmer_counts.jf
jellyfish histo -h 325000005 pvir_jelly_kmer_counts.jf > jelly_$HIST_SUFFIX

cp -r /scratch/beegfs/monthly/kjaron/genomic-features-of-asexual-animals/data/Pvir1/trimmed_reads .
