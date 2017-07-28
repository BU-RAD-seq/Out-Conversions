#!/usr/bin/env python3

#####################################################################
##
## out2parseclusters.py
##
## Version 1.00 -- 10 July 2017
##
## Created by Jeffrey DaCosta
## Copyright (c) 2011-2017 Boston University. All rights reserved.
##
## This Python (v3) script parses a defined list (or inverse of list)
## of clusters to a separate .out file
##
## This script is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
######################################################################

import sys, os, argparse
from argparse import RawTextHelpFormatter

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description=
                                 'This Python (v3) script parses a defined list (or inverse of list)\n'+
                                 'of clusters to a separate .out file\n\n'+
                                 'This script is free and distributed WITHOUT warranty; without\n'+
                                 'even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n'+
                                 'PARTICULAR PURPOSE.', formatter_class=RawTextHelpFormatter)

parser.add_argument('-i', type=str, metavar='infile', required=True, help='Name of input .out file for cluster parsing')
parser.add_argument('-o', type=str, metavar='outfile', required=True, help='Name of output .out file for parsed clusters')
parser.add_argument('-ns', type=int, metavar='num_samples', required=True, help='Number of samples in input outfile')
parser.add_argument('-l', type=str, metavar='cluster_list', required=True, help='Name of text file containing list of target clusters')
parser.add_argument('-inv',type=str, metavar='inverse_list', default='False', help='True or False: Parse inverse of list [False]')
args = parser.parse_args()

outfile = open(args.o,'w')

#count number of clusters in input file
print('Counting number of clusters in input file')
infile = open(args.i,'r')
file = infile.read()
num_clusters = file.count('Clstr')
print('Found '+str(num_clusters)+' clusters in input .out file')
infile.close()

parse_list = []
#count number of clusters to parse
#add clusters to parse to list
print('\nGathering information on clusters to parse')
line_count = 0
parsefile = open(args.l,'r')
for line in parsefile:
    parse_list.append(str(line.strip('\n')))
    line_count += 1
print('Found '+str(line_count)+' clusters in list')
parsefile.close()

parse_count = 0
inverse_count = 0
print('\nParsing clusters')
#loop for the range of the number of clusters
infile = open(args.i,'r')
for i in range(int(num_clusters)):
    #define cluster number, skip line 2
    header1 = infile.readline()
    header1s = header1.split()
    cluster = header1s[1]

    if str(cluster) in parse_list:
        if args.inv == 'False':
            parse_count += 1
            outfile.write(header1)
            for i in range((args.ns*2)+1):
                data = infile.readline()
                outfile.write(data)
        else:
            for i in range((args.ns*2)+1):
                data = infile.readline()

    else:
        if args.inv == 'True':
            inverse_count += 1
            outfile.write(header1)
            for i in range((args.ns*2)+1):
                data = infile.readline()
                outfile.write(data)
        else:
            for i in range((args.ns*2)+1):
                data = infile.readline()

infile.close()
outfile.close()

if args.inv == 'False':
    print('\n'+str(parse_count)+' clusters written to '+args.o)
else:
    print('\n'+str(inverse_count)+' clusters written to '+args.o)

print('\nFinished!!\n')
