#!/usr/bin/env bash

TARGET=$1
shift 1

snakemake $TARGET -p --jobs 10 $@ --cluster "bsub "\
"    -J {rule} "\
"    -q bgee "\
"    -n {threads} "\
"    -M {resources.M} "\
"    -R \"span[hosts=1]\" "\
"    -o logs/cluster/{rule}.{wildcards}.out "\
"    -e logs/cluster/{rule}.{wildcards}.err"