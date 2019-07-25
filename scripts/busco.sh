#/bin/bash

# this part is required only for Vital-it (if you have running BUSCO, delete it)
module add SequenceAnalysis/HMM-Profile/hmmer/3.1b2
module add SequenceAnalysis/GenePrediction/augustus/3.2.3
module add Blast/ncbi-blast/2.2.31+

export AUGUSTUS_CONFIG_PATH=/home/kjaron/src/busco-master/augustus_config

if [ -z "$TMPDIR" ]; then
	TMPDIR="."
fi

# run busco
# 1 - gzipped genome
# 2 - database with metazons
# 3 - output directory

TMPGENOME=$(mktemp -p $(echo $TMPDIR))
zcat $1 > $TMPGENOME

run_BUSCO.py -i $TMPGENOME -o busco -m geno -l $2 -c 16

mv run_busco $3

rm $TMPGENOME
