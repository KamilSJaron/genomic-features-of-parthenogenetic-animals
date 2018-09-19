#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 18:17:38 2017
@author: jklopfen
"""
import sys

# This script will parse the .collinearity file of MCScanX output and look for possible palindromes.
# Palindromes are defined here by having inverted collinear blocks within the same scaffold.

# Test of input file
if ".collinearity" not in sys.argv[1]:
    print("Input file must be a MCScanX .collinearity result")
    sys.exit()

# Generation of the name of the output file
in_name = sys.argv[1].split(".")
out_name = in_name[0] + "_palindromes.collinearity"

# I/O setup
informations = True
write_alignment = False
tmp = []
align_count = 0

# I/O
with open(sys.argv[1], 'r') as in_file,\
 open(out_name, 'w') as out_file:
     for line in in_file:
         if "Alignment" in line:
             write_alignment = False
             informations = False
             tmp = line.rstrip()
             tmp = tmp.split(" ")
             if tmp[-1] == 'minus':
                 tmp = tmp[-2].split("&")
                 if tmp[0] == tmp[1]:
                     align_count += 1
                     write_alignment = True
                 else : write_alignment = False
             else: write_alignment = False
         
         if(informations):
             out_file.write(line)
             
         if(write_alignment):
             out_file.write(line)

print("The .collinearity file was correctly parsed.", align_count, "palindromes have been identified.")
