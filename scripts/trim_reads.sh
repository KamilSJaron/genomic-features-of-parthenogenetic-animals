#!/bin/bash

# 2 output PATTERN.trimmed
# 1 are R1 and R2 reads ( expanded using [] )

DIR=$(dirname $1)

skewer -z -m pe -n -q 26 -l 48 -t 8 -o "$DIR"/reads $1 $2
