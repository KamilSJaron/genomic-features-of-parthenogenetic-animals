#!/bin/bash

module add UHTS/Analysis/samtools/1.3

samtools sort -n $1 | samtools view -h - | samblaster -M | samtools sort -O bam - > $2
