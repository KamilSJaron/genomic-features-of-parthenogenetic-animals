#/bin/bash

# this part is required only for Vital-it (if you have running BUSCO, delete it)
export PATH=\$PATH:/scratch/beegfs/monthly/ptranvan/Software/busco/2.0/
module add Blast/ncbi-blast/2.2.31+
module add SequenceAnalysis/HMM-Profile/hmmer/3.1b2
module add SequenceAnalysis/GenePrediction/augustus/3.2.2
export AUGUSTUS_CONFIG_PATH=/scratch/beegfs/monthly/kjaron/augustus_config

# run busco
# 1 - genome
# 2 - database with metazons
# 3 - output directory
BUSCO.py -i $1 -o $3 -m geno -l $2 -c 32