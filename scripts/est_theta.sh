#!/bin/bash
#
# is SP
# is a reference genome (VERSION)
# individual
# window_size

SAMPLE=$1
REF=$2
WINDOW=$3

BAM=$4
OUTPUT=$5

atlas task=estimateTheta \
	bam=$BAM \
	window=$WINDOW \
	suppressWarnings verbose

mv map_to_"$REF"_theta_estimates.txt $5