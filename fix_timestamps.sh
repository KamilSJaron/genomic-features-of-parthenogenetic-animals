#!/bin/bash

# this script is for cluster that are deleting old data
# the lifespan of data can be extended by renwing timestamps, but that screws up the workflow
# so, if you know that everything is cool, just timestamps are disorganised, you can execute this script
# that will renew timestemps in correct order

find . -name "*bam" -exec touch {} \;

sleep 1

find . -name "*.bai"  -exec touch {} \;

sleep 1

find . -name "*txt" -exec touch {} \;

