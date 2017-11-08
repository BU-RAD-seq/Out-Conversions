#!/usr/bin/env python3

# Stolen from: http://bugs.python.org/issue12806
import argparse
import re
import textwrap
class FlexiFormatter(argparse.RawTextHelpFormatter):
    def _split_lines(self, text, width):
        lines = list()
        main_indent = len(re.match(r'( *)',text).group(1))
        for line in text.splitlines():
            indent = len(re.match(r'( *)',line).group(1))
            list_match = re.match(r'( *)(([*-+>]+|\w+\)|\w+\.) +)',line)
            if(list_match):
                sub_indent = indent + len(list_match.group(2))
            else:
                sub_indent = indent
            line = self._whitespace_matcher.sub(' ', line).strip()
            new_lines = textwrap.wrap(
                text=line,
                width=width,
                initial_indent=' '*(indent-main_indent),
                subsequent_indent=' '*(sub_indent-main_indent),)
            lines.extend(new_lines or [' '])
        return lines

import os, sys, random, argparse
from argparse import RawTextHelpFormatter

print()
#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description=
                                 'This Python (v3) script converts sequences in an out file (with select\n'+
                                 'clusters) to the fasta format. It assumes that all clusters are\n'+
                                 'located on autosomes (hence "A" in script name). A separate fasta\n'+
                                 'file is created for each cluster. Samples with a "-9" in the population\n'+
                                 'column of the sample info file are skipped. Missing data are written as\n'+
                                 'a single N.\n\n'+
                                 'This script is free and distributed WITHOUT warranty; without\n'+
                                 'even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n'+
                                 'PARTICULAR PURPOSE.', formatter_class=FlexiFormatter)

parser.add_argument('-i', type=str, metavar='infile', required=True, help='Name of input .out file with selected clusters.')
parser.add_argument('-base', type=str, metavar='basename', required=True, help='Base name of fasdta files for individual clusters. '+
                    'Names of these files will have the format basename_clstr_#.fasta')
parser.add_argument('-si', type=str, metavar='infofile', required=True, help='Name of sample info file.')
parser.add_argument('-na', type=int, metavar='num_alleles', required=True, help='Number of alleles per sample to write to fasta file. '+
                    '1=major allele for low depth or flagged genotypes, random draw of first or second allele for good genotypes; '+
                    '2=alleles assigned to "a" or "b".')
parser.add_argument('-hemi', type=int, metavar='hemizygous_genotypes', required=True, help='Allow (1) or do not allow (0) hemizygous '+
                    'genotypes in fasta files. Applies only when na=2. If 1 then major allele will be written for low depth and '+
                    'flagged genotypes and second allele will contain a single N characters. If 0 then both alleles will contain a '+
                    'single N character')
args = parser.parse_args()

#count number of clusters
print('\nCounting number of clusters')
infile = open(args.i,'r')
file = infile.read()
num_clusters = file.count('Clstr')
print('Found '+str(num_clusters)+' clusters\n')
infile.close()

#gather sample info from popfile, skip pop=-9
print('\nGathering info from sample info file, skipping samples where population is -9:')
infofile = open(args.si,'r')
header = infofile.readline()
num_samples = 0
incl_samples = 0
samplearray = []
for line in infofile:
    pop = line.split() [4]
    samplearray.append(pop)
    if pop != '-9':
        num_samples += 1
        incl_samples += 1
    else:
        num_samples += 1
infofile.close()
print('Found '+str(num_samples)+' samples, of which '+str(incl_samples)+' will be included in output file')

percent10 = round(num_clusters*0.1,0)
cluster_count = 0
target = percent10
percent = 10

print('\nGathering cluster data, writing fasta file for each cluster\n\nAnalyzed:\n')
infile = open(args.i,'r')
for z in range(num_clusters):
    #define cluster number, skip line 2
    cluster = infile.readline().split() [1]
    header2 = infile.readline()

    locus_array = []
    for k in range(num_samples):
        if samplearray[k] != '-9':
            data = infile.readline().split()
            locus_array.append(data)
            data = infile.readline().split()
            locus_array.append(data)
        else:
            infile.readline()
            infile.readline()

    fasta_file = open(args.base+'_clstr_'+cluster+'.fasta','w')
    
    #randomly define allele 1 and 2, write to individual nexus file or interleave file
    for i in range(incl_samples):
        x = random.randint(1,2)

        if int(locus_array[i*2][7][0]) == 0:
            if args.na == 1:
                fasta_file.write('>'+locus_array[i*2][0]+'\nN\n')
            else: #args.na == 2
                fasta_file.write('>'+locus_array[i*2][0]+'a\nN\n')
                fasta_file.write('>'+locus_array[i*2+1][0]+'b\nN\n')

        elif int(locus_array[i*2][7][0]) == 1:
            if args.na == 1:
                if x ==1:
                    fasta_file.write('>'+locus_array[i*2][0]+'\n'+locus_array[i*2][1]+'\n')
                else:
                    fasta_file.write('>'+locus_array[i*2+1][0]+'\n'+locus_array[i*2+1][1]+'\n')
            else: #args.na == 2
                fasta_file.write('>'+locus_array[i*2][0]+'a\n'+locus_array[i*2][1]+'\n')
                fasta_file.write('>'+locus_array[i*2+1][0]+'b\n'+locus_array[i*2+1][1]+'\n')

        elif int(locus_array[i*2][7][0]) > 1:
            if args.hemi == 1:
                if args.na == 1:
                    fasta_file.write('>'+locus_array[i*2][0]+'\n'+locus_array[i*2][1]+'\n')
                else: #args.na == 2
                    fasta_file.write('>'+locus_array[i*2][0]+'a\n'+locus_array[i*2][1]+'\n')
                    fasta_file.write('>'+locus_array[i*2+1][0]+'b\nN\n')
            else:
                if args.na == 1:
                    fasta_file.write('>'+locus_array[i*2][0]+'\nN\n')
                else: #args.na == 2
                    fasta_file.write('>'+locus_array[i*2][0]+'a\nN\n')
                    fasta_file.write('>'+locus_array[i*2+1][0]+'b\nN\n')
                   
    fasta_file.close()

    cluster_count += 1
    if cluster_count == target:
        print(str(cluster_count)+' clusters ~ '+str(percent)+'%')
        target = target + percent10
        percent = percent + 10

infile.close()
print('\nFinished!!\n')
