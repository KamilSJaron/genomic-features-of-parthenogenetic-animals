#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 18:17:38 2017
@author: jklopfen, kjaron
"""
import sys

# data/Hduj1/MCScanX_biggap/Hduj1_prot.collinearity

# This script will parse the .collinearity file of MCScanX output and look for possible palindromes.
# Palindromes are defined here by having inverted collinear blocks within the same scaffold.

# Test of input file
if ".collinearity" not in sys.argv[1]:
    print("Input file must be a MCScanX .collinearity result")
    sys.exit()

# Generation of the name of the output file
in_name = sys.argv[1]
out_name = in_name.split(".")[0] + "_palindromes.collinearity"
gff_name = in_name.split(".")[0] + ".gff"

# I/O setup
write_alignment = False
align_count = 0

# load gff file
annotation_file = open(gff_name, 'r')

annotation = dict()
for line in annotation_file.readlines():
    spitline = line.replace("\n", "").split('\t')
    annotation[spitline[1]] = (int(spitline[2]), int(spitline[3]))

annotation_file.close()

def get_gap_within(genes):
    gaps = 0
    last_gene = genes[0][1]
    for gene in genes[1:]:
        gaps += abs(last_gene - gene[0])
        last_gene = gene[1]
    return(str(gaps))


def parse_alignment(buffer):
    gene_pos_cl1 = []
    gene_pos_cl2 = []
    for line in buffer:
        aln = line.split('\t')
        gene_pos_cl1.append(annotation[aln[1]])
        gene_pos_cl2.append(annotation[aln[2]])
    if min(min(gene_pos_cl1)) < min(min(gene_pos_cl2)):
        gap_between = min(min(gene_pos_cl2)) - max(max(gene_pos_cl1))
    else:
        gap_between = min(min(gene_pos_cl1)) - max(max(gene_pos_cl2))
    return([str(gap_between), get_gap_within(gene_pos_cl1), get_gap_within(gene_pos_cl2)])

def is_palindome(header):
    alignment = header.rstrip().split(" ")
    scfs = alignment[-2].split("&")
    return(alignment[-1] == 'minus' and scfs[0] == scfs[1])

palindrome_distances = []

with open(in_name, 'r') as in_file, open(out_name, 'w') as out_file:
    alignment_header = ""
    while not "Alignment" in alignment_header:
        alignment_header = in_file.readline()
    alignment_buffer = []

    #for line in in_file:
    for line in in_file:
        if "Alignment" in line:
            if is_palindome(alignment_header):
                palindrome_distances.append("\t".join(parse_alignment(alignment_buffer)))
            alignment_header = line
            write_alignment = False
            alignment_buffer = []
            if is_palindome(line):
                align_count += 1
                write_alignment = True
        else:
            alignment_buffer.append(line)
        if write_alignment:
            out_file.write(line)

print("Parsed:\t", in_name)
print("Created:\t", out_name)
print("Palindromes found:\t", align_count)
print("Palindrome table")
print("spacer\tblock1_internal_gaps\tblock2_internal_gaps")
print("\n".join(palindrome_distances))