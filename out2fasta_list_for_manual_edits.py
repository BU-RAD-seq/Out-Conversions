#!/usr/bin/env python3

##################################
##
## out2fastas_for_manual_edits.py
##
## Version 1.00 -- 28 July 2017
##
## Created by Jeffrey DaCosta
## Copyright (c) 2011-2017 Boston University. All rights reserved.
##
## This Python (v3) script produces a fasta file for each cluster in a
## defined list. The resulting fasta files can be imported to Geneious
## or a similar program to manually revise cluster alignments. The two
## alleles for each sample are assigned to "a" and "b" sequences in the
## fasta file. For samples with missing data (i.e., genotype code=0) a
## single "N" is written for each allele. For samples with low depth
## (code=2) a single "N" is written for the second allele.
##
## Output fasta file names begin with the provided "base" name followed
## by _cluster#.fasta
##
## This script is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
##################################

import sys, os, random, argparse
from argparse import RawTextHelpFormatter

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description=
                                 'This Python (v3) script produces a fasta file for each cluster in a\n'+
                                 'defined list. The resulting fasta files can be imported to Geneious\n'+
                                 'or a similar program to manually revise cluster alignments. The two\n'+
                                 'alleles for each sample are assigned to "a" and "b" sequences in the\n'+
                                 'fasta file. For samples with missing data (i.e., genotype code=0) a\n'+
                                 'single "N" is written for each allele. For samples with low depth\n'+
                                 '(code=2) a single "N" is written for the second allele.\n\n'+
                                 'Output fasta file names begin with the provided "base" name followed\n'+
                                 'by _cluster#.fasta\n\n'+
                                 'This script is free and distributed WITHOUT warranty; without even\n'+
                                 'the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR\n'+
                                 'PURPOSE.', formatter_class=RawTextHelpFormatter)

parser.add_argument('-i', type=str, metavar='infile', required=True, help='Name of input .out file')
parser.add_argument('-base', type=str, metavar='outfile_basename', required=True, help='Base name for output fasta files')
parser.add_argument('-ns', type=int, metavar='num_samples', required=True, help='Number of samples in input outfile')
parser.add_argument('-l', type=str, metavar='cluster_list', required=True, help='Name of text file containing list of target clusters')
args = parser.parse_args()

#count number of clusters
print('Counting number of clusters in .out file')
infile = open(args.i,'r')
file = infile.read()
num_clusters = file.count('Clstr')
infile.close()
print('Found '+str(num_clusters)+' clusters\n')

#gather target clusters
print('\nGathering target clusters')
listfile = open(args.l,'r')
targets = []
for line in listfile:
    cluster = line.strip('\n')
    targets.append(cluster)
listfile.close()
print('Found '+str(len(targets))+' target clusters')

print('\nCreating fasta files for target clusters')
infile = open(args.i,'r')
for z in range(num_clusters):
    #define cluster number, skip line 2
    header1 = infile.readline()
    header1 = header1.split()
    cluster = header1[1]
    header2 = infile.readline()

    if cluster in targets:
        outfile = open(args.base+'_'+cluster+'.fasta','w')
        for i in range(args.ns):
            allele1 = infile.readline().split()
            if int(allele1[7]) == 0:                                                        #if genotype code=0 (missing)
                outfile.write('>'+allele1[0]+'a\nN\n>'+allele1[0]+'b\nN\n')
                infile.readline()
            elif int(allele1[7]) == 2:                                                      #if genotype code=2 (low depth)
                outfile.write('>'+allele1[0]+'a\n'+allele1[1]+'\n>'+allele1[0]+'b\nN\n')
                infile.readline()
            else:
                outfile.write('>'+allele1[0]+'a\n'+allele1[1]+'\n')
                allele2 = infile.readline().split()
                outfile.write('>'+allele2[0]+'b\n'+allele2[1]+'\n')
        outfile.close()
    else:
        for i in range(args.ns*2):
            data = infile.readline()

infile.close()
print('\nFinished!!\n')
