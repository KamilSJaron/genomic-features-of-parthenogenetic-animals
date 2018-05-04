#/bin/bash

# this part is required only for Vital-it (if you have running BUSCO, delete it)
module add SequenceAnalysis/HMM-Profile/hmmer/3.1b2
module add SequenceAnalysis/GenePrediction/augustus/3.2.3
export AUGUSTUS_CONFIG_PATH=/home/kjaron/src/busco-master/augustus_config

# run busco
# 1 - gzipped genome
# 2 - database with metazons
# 3 - output directory
zcat $1 > genome.fa

run_BUSCO.py -i genome.fa -o busco -m geno -l $2 -c 16

mv run_busco $3

rm genome.fa
