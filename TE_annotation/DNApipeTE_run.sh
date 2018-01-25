#!/bin/bash
#script to run dnaPipeTE

#get dnaPipeTE from https://github.com/clemgoub/dnaPipeTE
#follow installation instructions

#here, I used version 1.2, which was setup by Patrick Tran Van
#copy from Patrick's path

#set up environment
#module add Development/java/latest
module add Development/java_jre/1.7.0_17
module add UHTS/Aligner/bowtie2/2.3.0
module add R/latest

#tmp has to be adjusted
mkdir -p /scratch/local/$USER/tmp
export TMPDIR=/scratch/local/$USER/tmp
export _JAVA_OPTIONS="-XX:ParallelGCThreads=6"

clear

#run the dnaPipeTE pipeline
#IMPORTANT: must be run from installed software path !!

#mabye better softlink to program folder if run on local scratch
cd /scratch/beegfs/monthly/jbast/software/dnaPipeTE_old/1.2
#cd /scratch/beegfs/monthly/ptranvan/Software/dnaPipeTE/1.2

#run with single-end reads (can be gzipped)
#IMPORTANT: give genome size of organism
#the genome_coverage and sample_number depends on how many TEs are expected, but the setting here should be generally ok

python3 ./dnaPipeTE.py -input /scratch/local/jbast/review_asex/Avaga/reads-trimmed-pair1.fastq.gz -output /scratch/local/jbast/review_asex/Avaga \
-cpu 12 -genome_size 244000000 -genome_coverage 0.5 -sample_number 3


module rm Development/java_jre/latest
module rm UHTS/Aligner/bowtie2/2.3.0
module rm R/latest

printf "****DONE****\n"
