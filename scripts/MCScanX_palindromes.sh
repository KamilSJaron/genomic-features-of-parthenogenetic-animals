#!/bin/bash

# this script was called in a command line
# for blastfile in data/*/MCScanX/*.blast; do
#  SP=$(echo $blastfile | cut -f 2 -d '/');
#  scripts/MCScanX_palindromes.sh $SP;
# done

SP=$1

DIR=data/$SP/MCScanX

##### GENERATE BLAST of ALL PROTEINS vs ALL PROTEINS

MCScanX -s 1 -m 0 $DIR/"$SP"_prot

grep "# Number of" "$DIR"/"$SP"_prot.collinearity > "$DIR"/"$SP"_prot.collinearity_summary.txt;
python3 scripts/colinearity_parser.py $DIR/"$SP"_prot.collinearity >> $DIR/"$SP"_prot.collinearity_summary.txt;
grep "Align" $DIR/"$SP"_prot.collinearity | wc -l >> $DIR/"$SP"_prot.collinearity_summary.txt;

rm -r "$DIR"/"$SP"_prot.html