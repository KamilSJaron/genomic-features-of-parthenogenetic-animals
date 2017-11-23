#!/bin/env bash

wget https://zenodo.org/record/399475/files/MIG-Phylogenomics.zip

unzip MIG-Phylogenomics.zip


### Find assemblies

find MIG-Phylogenomics -name "*.fa*"

# | grep "genome"




rm MIG-Phylogenomics.zip
