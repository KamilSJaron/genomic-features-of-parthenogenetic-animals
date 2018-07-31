#!/bin/bash

# 3 output PATTERN.trimmed
# 1, 2 are R1 and R2 reads ( expanded using [] )

skewer -z -m pe -n -q 26 -l 21 -t 8 -o $3 $1 $2
