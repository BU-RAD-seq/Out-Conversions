#!/usr/bin/env python3

##################################
##
## out2phylipA.py
##
## Version 1.00 -- 28 May 2019
##
## Created by Jeffrey DaCosta
## Copyright (c) 2019 Boston College. All rights reserved.
##
## This Python (v3) script converts sequences in an out file (with select
## clusters) to the phylip format. It assumes that all clusters are
## located on autosomes (hence "A" in script name). A separate phylip
## file is created for each cluster, and a phylip file with concatenated
## sequences is also created. Samples with a "-9" in the population
## column of the sample info file are skipped.
##
## This script is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
##################################

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

import os, sys, math, random, subprocess, argparse
from argparse import RawTextHelpFormatter

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description=
                                 'This Python (v3) script converts sequences in an out file (with select\n'+
                                 'clusters) to the phylip format. It assumes that all clusters are\n'+
                                 'located on autosomes (hence "A" in script name). A separate phylip\n'+
                                 'file is created for each cluster, and a phylip file with concatenated\n'+
                                 'sequences is also created. Samples with a "-9" in the population\n'+
                                 'column of the sample info file are skipped.\n\n'+
                                 'This script is free and distributed WITHOUT warranty; without\n'+
                                 'even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n'+
                                 'PARTICULAR PURPOSE.', formatter_class=FlexiFormatter)

requiredParam = parser.add_argument_group('required parameters')
requiredParam.add_argument('-i', type=str, metavar='infile', required=True, help='Name of input .out file with selected clusters.')
requiredParam.add_argument('-base', type=str, metavar='basename', required=True, help='Base name of phylip files for individual clusters. '+
                    'Names of these files will have the format basename_clstr_#.nex')
requiredParam.add_argument('-cat', type=str, metavar='catfile', required=True, help='Name of phylip file with concatenated sequences.')
requiredParam.add_argument('-si', type=str, metavar='infofile', required=True, help='Name of sample info file.')
requiredParam.add_argument('-na', type=int, metavar='num_alleles', required=True, help='Number of alleles per sample to write to phylip file. '+
                    '1=major allele for low depth or flagged genotypes, random draw of first or second allele for good genotypes; '+
                    '2=alleles  randomly assigned to an "a" or "b".')
requiredParam.add_argument('-hemi', type=int, metavar='hemizygous_genotypes', required=True, help='Allow (1) or do not allow (0) hemizygous '+
                    'genotypes in phylip files. Applies only when na=2. If 1 then major allele will be written for low depth and '+
                    'flagged genotypes and second allele will contain a string of ?. If 0 then both alleles will contain a string '+
                    'of ?')
args = parser.parse_args()
           
#count number of clusters
print('Counting number of clusters:')
infile = open(args.i,'r')
file = infile.read()
num_clusters = file.count('Clstr')
print('Found '+str(num_clusters)+' clusters')
infile.close()

#gather sample info from popfile, skip pop=-9
print('\nGathering info from sample info file, skipping samples where population is -9:')
infofile = open(args.si,'r')
header = infofile.readline()
num_samples = 0
incl_samples = 0
incl_samples_list = []
samplearray = []
for line in infofile:
    pop = line.split() [4]
    samplearray.append(pop)
    if pop != '-9':
        num_samples += 1
        incl_samples += 1
        incl_samples_list.append(line.split() [1])
    else:
        num_samples += 1
infofile.close()
print('Found '+str(num_samples)+' samples, of which '+str(incl_samples)+' will be included in output file')

#set values for printing progress to screen
percent10 = round(num_clusters*0.1,0)
cluster_count = 0
target = percent10
percent = 10

#start count of total length
total_length = 0

#open concatenated phylip file
concatfile = open(args.cat,'w')

#open list for concatenated data
concat = []
if args.na == 1:
    for s in range(incl_samples):
        concat.append('')
else:
    for s in range(incl_samples*2):
        concat.append('')

print('\nGathering cluster data, writing nexus file for each cluster\n\nAnalyzed:\n')
#loop for the range of the number of clusters
infile = open(args.i,'r')
for z in range(int(num_clusters)):
    #define cluster number, skip line 2
    header1 = infile.readline()
    header1 = header1.split()
    cluster = header1[1]
    header2 = infile.readline()

    #create array of data for each locus
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

    max_length = 0
    for y in range(len(locus_array)):
        length = len(locus_array[y][1])
        if length > max_length:
            max_length = length

    #make appropriate sized string of ? for missing
    missing = max_length*'?'
        
    if args.na == 1:
        #open and write header of individual cluster phylip file
        phy_file = open(args.base+'_clstr_'+cluster+'.phy','w')
        phy_file.write(str(incl_samples)+'\t'+str(max_length)+'\n')
        
    elif args.na == 2:
        #open and write header of individual cluster phylip file
        phy_file = open(args.base+'_clstr_'+cluster+'.phy','w')
        phy_file.write(str(incl_samples*2)+'\t'+str(max_length)+'\n')

    else:
        print('\nERROR: number of alleles must equal 1 or 2!\n')

    #randomly define allele 1 and 2, write to individual phylip file or append to concat file
    for i in range(int(incl_samples)):
        x = random.randint(1,2)
        if int(locus_array[i*2][7][0]) == 0:
            if args.na == 1:
                phy_file.write(locus_array[i*2][0]+'\t'+missing+'\n')
                concat[i] += missing
            else: #args.na == 2
                phy_file.write(locus_array[i*2][0]+'a\t'+missing+'\n')
                phy_file.write(locus_array[i*2+1][0]+'b\t'+missing+'\n')
                concat[i*2] += missing
                concat[i*2+1] += missing
        elif int(locus_array[i*2][7][0]) == 1:
            if x == 1:
                if args.na == 1:
                    phy_file.write(locus_array[i*2][0]+'\t'+locus_array[i*2][1]+'\n')
                    concat[i] += locus_array[i*2][1]
                else: #args.na == 2
                    phy_file.write(locus_array[i*2][0]+'a\t'+locus_array[i*2][1]+'\n')
                    phy_file.write(locus_array[i*2+1][0]+'b\t'+locus_array[i*2+1][1]+'\n')
                    concat[i*2] += locus_array[i*2][1]
                    concat[i*2+1] += locus_array[i*2+1][1]
            else:
                if args.na == 1:
                    phy_file.write(locus_array[i*2+1][0]+'\t'+locus_array[i*2+1][1]+'\n')
                    concat[i] += locus_array[i*2+1][1]
                else: #args.na == 2
                    phy_file.write(locus_array[i*2+1][0]+'a\t'+locus_array[i*2+1][1]+'\n')
                    phy_file.write(locus_array[i*2][0]+'b\t'+locus_array[i*2][1]+'\n')
                    concat[i*2+1] += locus_array[i*2+1][1]
                    concat[i*2] += locus_array[i*2][1]
        elif int(locus_array[i*2][7][0]) > 1:
            if x == 1:
                if args.hemi == 0:
                    if args.na == 1:
                        phy_file.write(locus_array[i*2][0]+'\t'+missing+'\n')
                        concat[i] += missing
                    else: #args.na == 2
                        phy_file.write(locus_array[i*2][0]+'a\t'+missing+'\n')
                        phy_file.write(locus_array[i*2+1][0]+'b\t'+missing+'\n')
                        concat[i*2] += missing
                        concat[i*2+1] += missing
                else: #args.hemi == 1
                    if args.na == 1:
                        phy_file.write(locus_array[i*2][0]+'\t'+locus_array[i*2][1]+'\n')
                        concat[i] += locus_array[i*2][1]
                    else: #args.na == 2
                        phy_file.write(locus_array[i*2][0]+'a\t'+locus_array[i*2][1]+'\n')
                        phy_file.write(locus_array[i*2+1][0]+'b\t'+missing+'\n')
                        concat[i*2] += locus_array[i*2][1]
                        concat[i*2+1] += missing
            else:
                if args.hemi == 0:
                    if args.na == 1:
                        phy_file.write(locus_array[i*2][0]+'\t'+missing+'\n')
                        concat[i] += +missing
                    else: #args.na == 2
                        phy_file.write(locus_array[i*2][0]+'a\t'+missing+'\n')
                        phy_file.write(locus_array[i*2+1][0]+'b\t'+missing+'\n')
                        concat[i*2] += missing
                        concat[i*2+1] += missing
                else: #args.hemi == 1
                    if args.na == 1:
                        phy_file.write(locus_array[i*2][0]+'\t'+locus_array[i*2][1]+'\n')
                        concat[i] += locus_array[i*2][1]
                    else: #args.na == 2
                        phy_file.write(locus_array[i*2+1][0]+'a\t'+missing+'\n')
                        phy_file.write(locus_array[i*2][0]+'b\t'+locus_array[i*2][1]+'\n')
                        concat[i*2] += missing
                        concat[i*2+1] += locus_array[i*2][1]
            
    #close individual phylip file
    phy_file.close()

    cluster_count += 1
    if cluster_count == target:
        print(str(cluster_count)+' clusters ~ '+str(percent)+'%')
        target = target + percent10
        percent = percent + 10

total_length = len(concat[0])

#write header of interleave file
if args.na == 1:
    concatfile.write(str(incl_samples)+'\t'+str(total_length)+'\n')
    for x in range(incl_samples):
        concatfile.write(incl_samples_list[x]+'\t'+concat[x]+'\n')
else: #args.na == 2
    concatfile.write(str(incl_samples*2)+'\t'+str(total_length)+'\n')
    for x in range(incl_samples):
        concatfile.write(incl_samples_list[x]+'a\t'+concat[x*2]+'\n')
        concatfile.write(incl_samples_list[x]+'b\t'+concat[x*2+1]+'\n')       
concatfile.close()

print('\nFinished!!\n\n')
