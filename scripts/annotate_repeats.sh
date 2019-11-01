#!/bin/bash

#script to run dnaPipeTE
#set up environment

# reads						data/Avag1/reads-trimmed-pair1.fastq.gz
# sample code				244000000
# output dir				data/Avag1/dnaPipeTE

READ_DIR=data/sPfal1/trimmed_reads
GENOME_SIZE=3528200000
OUTPUT=data/sPfal1/dnaPipeTE
SP_DIR=$(dirname "$OUTPUT")
SHARED_DIR=$(pwd)

#tmp has to be adjusted -> already in the wrapper script, maybe I should veify
mkdir -p "$SHARED_DIR"/"$SP_DIR"/temp
mkdir -p "$SHARED_DIR"/"$SP_DIR"/trimmed_reads
mkdir -p "$SHARED_DIR"/"$OUTPUT"
export TMPDIR="$SHARED_DIR"/"$SP_DIR"/temp
export _JAVA_OPTIONS="-XX:ParallelGCThreads=24"

cat "$READ_DIR"/*fastq.gz > "$SHARED_DIR"/"$READ_DIR"/all_reads.fastq.gz

#run the dnaPipeTE pipeline
#IMPORTANT: must be run from installed software path !!


#mabye better softlink to program folder if run on local scratch
# dnaPipeTE requires be run from it's own dir, uncomment the next line with appropriate path in your cluster
# cd ...

#run with single-end reads (can be gzipped)
#IMPORTANT: give genome size of organism
#the genome_coverage and sample_number depends on how many TEs are expected, but the setting here should be generally ok

python3 ./dnaPipeTE.py -input "$SHARED_DIR"/"$READ_DIR"/all_reads.fastq.gz \
    -output "$SHARED_DIR"/"$OUTPUT" \
    -cpu 32 -genome_size "$GENOME_SIZE" -genome_coverage 0.5 -sample_number 3

cd "$SHARED_DIR"

rm "$READ_DIR"/all_reads.fastq.gz
rm -r "$SP_DIR"/temp
rmdir "$SP_DIR"/trimmed_reads
rmdir "$SP_DIR"

printf "****DONE****\n"
