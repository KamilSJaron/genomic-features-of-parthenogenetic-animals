#!/bin/bash

sp=$1

export PATH=$PATH:/scratch/beegfs/monthly/ptranvan/Software/KAT/2.4.1/kat_tmp/bin
export PYTHONPATH=$PYTHONPATH:/software/lib64/python3.5/site-packages

mkdir -p data/$sp/KAT

kat comp -o data/$sp/KAT/kmer_profile_in_genome -t 8 data/"$sp"'/trimmed_reads/*q.gz' data/$sp/genome.fa.gz

