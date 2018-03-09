#!/usr/bin/env python3

##################################
##
## out2fineRADstructure.py
##
## Version 1.00 -- 23 July 2017
##
## Created by Michael Sorenson and Jeffrey DaCosta
## Copyright (c) 2017 Boston University. All rights reserved.
##
## This Python (v3) converts an out file (with select clusters)
## to an input file for fineRADstructure. It assumes that all 
## clusters are located on autosomes.
##
## This script is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
##################################

#format RAD-seq data for fineRADstructure
#assume hemizygous data are not OK
	#thus, use only genotypes == 1
#assume gap chars and indels are not OK
	#recode as pseudo-SNPs

import os, sys, argparse
from argparse import RawTextHelpFormatter

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description=
                                 'This Python (v3) converts an out (with select clusters) \n'+
                                 'to an input file for fineRADstructure. It assumes that all \n'+
                                 'clusters are located on autosomes.\n\n'+
                                 'This script is free and distributed WITHOUT warranty; without\n'+
                                 'even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n'+
                                 'PARTICULAR PURPOSE.', formatter_class=RawTextHelpFormatter)

parser.add_argument('-i', type=str, metavar='infile', required=True, help='Name of input .out file with selected clusters')
parser.add_argument('-o', type=str, metavar='outfile', required=True, help='Name of output fineRADstructure file')
parser.add_argument('-si', type=str, metavar='infofile', required=True, help='Name of sample info file')
args = parser.parse_args()

def transform(genotypes):
    all_genos = [item for sublist in genotypes for item in sublist]
    uni_genos = list(set(all_genos))
    if '' in uni_genos:
        uni_genos.pop(uni_genos.index(''))
    if uni_genos == ['.']:
        genotypes=[]
    else:
        gaps = False
        indels = False
        for i in uni_genos:
            if '-' in i:
                gaps = True
            elif '0' in i:
                indels = True
        if gaps == True:
            #find columns with gaps
            gap_pos = []
            for i in uni_genos:
                temp = list(i)
                if '-' in temp:
                    gap_pos.append([j for j in range(len(temp)) if temp[j] == '-' ])
            gap_pos = [item for sublist in gap_pos for item in sublist]
            gap_pos = list(set(gap_pos))
            gap_pos.sort()
                       
            #find common base in those columns
            acgt = ['A','C','G','T']
            common = []
            for i in gap_pos:
                states = [all_genos[j][i] for j in range(len(all_genos)) if all_genos[j] != '']
                counts = [states.count(j) for j in acgt]
                common.append(acgt[counts.index(max(counts))])
            
            #change gaps in those columns to common base
            for i,j in zip(gap_pos,common):
                for x in range(len(genotypes)):
                    if genotypes[x][0] != '':
                        if genotypes[x][0][i] == '-':
                            genotypes[x][0] = genotypes[x][0][:i] + j + genotypes[x][0][i+1:]
                        if genotypes[x][1][i] == '-':
                            genotypes[x][1] = genotypes[x][1][:i] + j + genotypes[x][1][i+1:]

            
        genotypes = [i[0]+'/'+i[1] if i[0] != i[1] else i[0] for i in genotypes]
        if indels == True:
            genotypes = [i.replace('0','A') for i in genotypes]
            genotypes = [i.replace('1','T') for i in genotypes]
       
    return(genotypes)

outfile = open(args.o,'w')

#count number of clusters in input file
print('Counting number of clusters in input file:')
infile = open(args.i,'r')
file = infile.read()
num_clusters = file.count('Clstr')
print('Found '+str(num_clusters)+' clusters in input .out file')
infile.close()

#gather sample info from popfile, skip pop=-9
print('\nGathering info from sample info file, skipping samples where population is -9:')
pops = []
incl_samples = []
infofile = open(args.si,'r')
header = infofile.readline()
sample_count = 0
for line in infofile:
    pop = line.split() [4]
    pops.append(pop)
    sample_count += 1
    if pop != '-9':
        incl_samples.append(line.split()[1])
infofile.close()
outfile.write('\t'.join(incl_samples)+'\n')
print('Found '+str(sample_count)+', of which '+str(len(incl_samples))+' will be included in output file')

#convert genotypes for fineRADstructure
print('\nConverting genotypes for use in fineRADstructure')
infile = open(args.i,'r')
for line in infile:
    if 'Clstr' in line:
        infile.readline()
        genos = []
        for i in range(sample_count):
            if pops[i] != '-9':
                data1 = infile.readline().split()
                data2 = infile.readline().split()
                genos.append([data1[2],data2[2],data1[7]])
            else:
                infile.readline()
                infile.readline()
        genos = [['',''] if i[2] != '1' else [i[0],i[1]] for i in genos]                
        genos = transform(genos)
        
        if genos != []:
            outfile.write('\t'.join(genos)+'\n')

outfile.close()

print('\nFinished!!\n\n')
            
