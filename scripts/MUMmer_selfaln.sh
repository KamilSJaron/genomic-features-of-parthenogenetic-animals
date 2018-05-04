#!/bin/bash

# MUMmer -ax asm5 genome.fa.gz genome.fa.gz > genome_to_Minc1.sam

mkdir $2
# TODO fix this so genome will contain a species name or it will be saved not in root of comp
zcat $1 > genome.fa

nucmer -t 16 --mum -p $2/selfaln genome.fa genome.fa
show-coords -THrcl $2/selfaln.delta > $2/selfaln.coords
dnadiff -p $2/selfaln -d $2/selfaln.delta