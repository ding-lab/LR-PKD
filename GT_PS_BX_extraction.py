#!/usr/bin/env python3
# Yize Li
import os
import sys
import numpy as np
import pysam
import pandas as pd
from sys import argv

"""
This script is to:
1. Extract the GT, PS, BX information from the VCF file
Usage: python3 <script> <10x VCF file>
E.g. python3 GT_PS_BX_extraction.py input.vcf

Args:
    <script>: the python script (required)
    <readcount file>: the file includes genotype, phase block and barcode information (required)
    Args need to follow this certain order
"""

# Save a new output file to merge with the MAF file
newfile = open(sys.argv[1] + "_GT_PS_BX.txt", "w")

# Make a list including all the existing allele (A, C, G, T etc) across all columns in the current file
chr_list = []
pos_list = []
GT_list = []
PS_list = []
BX_list = []

f = open(sys.argv[1], "r")
for line in f: # Loop through each row (mutation call)
    line = line.strip() # Remove "\n"
    # Skip the header which comment sign (#)
    if not line.startswith("#"):
        for i in range(len(line.split("\t"))):
            if i == 0:
                chr_list.append(line.split("\t")[0])
            elif i == 1:
                pos_list.append(line.split("\t")[1])
            elif i == 8:
                format_list = []
                for j in range(len((line.split("\t")[8]).split(":"))):
                    format_list.append((line.split("\t")[8]).split(":")[j])
                indices_GT = [k for k, s in enumerate(format_list) if 'GT' in s][0]
                indices_PS = [k for k, s in enumerate(format_list) if 'PS' in s][0]
                indices_BX = [k for k, s in enumerate(format_list) if 'BX' in s][0]
            elif i == 9:
                GT_i = (line.split("\t")[9]).split(":")[indices_GT]
                PS_i = (line.split("\t")[9]).split(":")[indices_PS]
                BX_i = (line.split("\t")[9]).split(":")[indices_BX]
                GT_list.append(GT_i)
                PS_list.append(PS_i)
                BX_list.append(BX_i)
            else:
                pass
f.close() # Close file

# Print the new file header
newfile.write("Chromosome" + "\t" + "Start_Position" + "\t" + "GT" + "\t" + "PS" + "\t" + "BX" + "\n")
for l in range(len(GT_list)):
    newfile.write(chr_list[l] + "\t" + pos_list[l] + "\t" + GT_list[l] + "\t" + PS_list[l] + "\t" + BX_list[l] + "\n")

file1 = pd.read_csv('input_d1.vep.maf', na_filter=False, float_precision=None, sep='\t', low_memory=False)
file2 = pd.read_csv('input.vep.vcf_GT_PS_BX.txt', na_filter=False, float_precision=None, sep='\t', low_memory=False)

merge_df = pd.merge(file1, file2, on=['Chromosome', 'Start_Position']).fillna(0) #, how='left'
merge_df.to_csv('input.vep.GTPSBX.temp.maf', sep='\t')
