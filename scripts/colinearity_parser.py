#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 18:17:38 2017
@author: jklopfen, kjaron
"""
import sys
import gzip

# data/Hduj1/MCScanX_biggap/Hduj1_prot.collinearity

# This script will parse the .collinearity file of MCScanX output and look for possible palindromes.
# Palindromes are defined here by having inverted collinear blocks within the same scaffold.

# Test of input file
if ".collinearity" not in sys.argv[1]:
    sys.stderr.write("Input file must be a MCScanX output (.collinearity)")
    sys.exit()

if ".gff3.gz" not in sys.argv[2]:
    sys.stderr.write("Input annotation must be a gzipped gff3 file (.gff3.gz)")
    sys.exit()

# Generation of the name of the output file
# in_name = 'data/Pfor1/MCScanX_biggap/Pfor1_prot.collinearity'
in_name = sys.argv[1]
# gff_name = 'data/Pfor1/annotation.gff3.gz'
gff_name = sys.argv[2]
out_name = in_name.split(".")[0] + "_palindromes.collinearity"

###########
# load gff file into two dictionaries
###########
# - transcript -> gene; which will be used to remove palindromes of alternative transcripts
# - gene -> position; which will be used to retrieve position of genes with respect to strand
annotation_file = gzip.GzipFile(gff_name,'r')
lines = annotation_file.readlines()

rna2gene = dict()
gene2pos = dict()

for line in lines:
    processed_line = line.decode()
    if processed_line[0] == "#":
        continue
    processed_line = processed_line.rstrip().split('\t')
    line_details = processed_line[8].split(';')
    if processed_line[2] == 'gene':
        gene_name = line_details[0].split('=')[1]
        gene2pos[gene_name] = (int(processed_line[3]), int(processed_line[4]), processed_line[6])
    if processed_line[2] == 'mRNA':
        transcript_name = line_details[0].split('=')[1]
        gene_name = line_details[1].split('=')[1]
        rna2gene[transcript_name] = gene_name

annotation_file.close()

### some parsing functions to make code bit cleaner

def get_gap_within(genes):
    gaps = 0
    last_gene = genes[0][1]
    for gene in genes[1:]:
        gaps += abs(last_gene - gene[0])
        last_gene = gene[1]
    return(str(gaps))

# this function should return buffer but just with valid palindrome sequences
def keep_palindomes(header, buffer):
    palindrome_buffer = []
    # if it's reverse tandem it has pontential to be a palindrome as well
    if is_rerverse_tandem(header):
        sys.stderr.write("Parsing: " + header)
        for line in buffer:
            aln = line.split('\t')
            gene_1 = rna2gene[aln[1]]
            gene_2 = rna2gene[aln[2]]
            if gene_1 == gene_2:
                sys.stderr.write("\tDiscarding:" + aln[1] + " and " + aln[2] + " are alternative transcripts of the same gene.\n")
                continue
            if gene2pos[gene_1][2] == gene2pos[gene_2][2]:
                sys.stderr.write("\tDiscarding:" + gene_1 + " and " + gene_2 + " are genes of the same strand (" + gene2pos[gene_2][2] + ")\n")
                continue
            palindrome_buffer.append(line)
    return(palindrome_buffer)


def parse_alignment(buffer):
    gene_pos_cl1 = []
    gene_pos_cl2 = []
    for line in buffer:
        aln = line.split('\t')
        gene_pos_cl1.append(gene2pos[rna2gene[aln[1]]][0:2])
        gene_pos_cl2.append(gene2pos[rna2gene[aln[2]]][0:2])
    if min(min(gene_pos_cl1)) < min(min(gene_pos_cl2)):
        gap_between = min(min(gene_pos_cl2)) - max(max(gene_pos_cl1))
    else:
        gap_between = min(min(gene_pos_cl1)) - max(max(gene_pos_cl2))
    return("\t".join([str(gap_between), get_gap_within(gene_pos_cl1), get_gap_within(gene_pos_cl2)]))

def is_rerverse_tandem(header):
    header = header.rstrip().split(" ")
    scfs = header[-2].split("&")
    return(header[-1] == 'minus' and scfs[0] == scfs[1])

### end functions ###

putative_palindrome_count = 0
palindrome_count = 0
palindrome_distances = []

in_file = open(in_name, 'r')
with open(in_name, 'r') as in_file, open(out_name, 'w') as out_file:
    # removing the initial lines with parameters
    for alignment_header in in_file:
        if "Alignment" in alignment_header:
            break
    alignment_buffer = []
    if is_rerverse_tandem(alignment_header):
        putative_palindrome_count += 1

    # iterating though alignment part of the file
    for line in in_file:
        if "Alignment" in line:
            alignment_buffer = keep_palindomes(alignment_header, alignment_buffer)
            if alignment_buffer:
                palindrome_count += 1
                palindrome_distances.append(parse_alignment(alignment_buffer))
                out_file.write(alignment_header)
                for aln in alignment_buffer:
                    out_file.write(aln)
            alignment_header = line
            if is_rerverse_tandem(alignment_header):
                putative_palindrome_count += 1
            alignment_buffer = []
        else:
            alignment_buffer.append(line)
    # if the last one is palindrome
    alignment_buffer = keep_palindomes(alignment_header, alignment_buffer)
    if alignment_buffer:
        palindrome_distances.append(parse_alignment(alignment_buffer))

print("Parsed:\t", in_name)
print("Created:\t", out_name)
print("Reverse synteny blocks found:\t", putative_palindrome_count)
print("Palindromes found:\t", palindrome_count)
print("Palindrome table")
print("spacer\tblock1_internal_gaps\tblock2_internal_gaps")
print("\n".join(palindrome_distances))