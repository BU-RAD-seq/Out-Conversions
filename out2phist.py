#!/usr/bin/env python3

##################################
##
## out2phist.py
##
## Version 1.03 -- 18 April 2014
##
## Created by Michael Sorenson and Jeffrey DaCosta
## Copyright (c) 2014 Boston University. All rights reserved.
##
## This Python (v3) script calculates locus-by-locus phi-st
## values from an out file containing filtered clusters from our
## ddRAD pipeline. Script assumes all loci are autosomal.
##
## This script is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
##################################

import os, sys, argparse, random
from argparse import RawTextHelpFormatter

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description=
                                 'This Python (v3) script calculates locus-by-locus phi-st\n'+
                                 'values from an out file containing filtered clusters from our\n'+
                                 'ddRAD pipeline. Script assumes all loci are autosomal\n\n'+
                                 'This script is free and distributed WITHOUT warranty; without\n'+
                                 'even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n'+
                                 'PARTICULAR PURPOSE.', formatter_class=RawTextHelpFormatter)

parser.add_argument('-i', type=str, metavar='infile', required=True, help='Name of input .out file with filtered clusters')
parser.add_argument('-o', type=str, metavar='outfile', required=True, help='Name of output SNPs file')
parser.add_argument('-si', type=str, metavar='infofile', required=True, help='Name of sample info file')
args = parser.parse_args()

#first three functions below are not currently implemented
def calculate_pvalue(seq_array,observed):
    count = 0
    for i in range(1000):
        randomize_array(seq_array)
        phist = calculate_phist(seq_array)
        if phist >= observed:
            count += 1
    p = count/1000
    return(p)

def calculate_random(seq_array):
    results =[]
    for i in range(10):
        seq_array = randomize_array(seq_array)
        phist = calculate_phist(seq_array)
        results.append(phist)
    return(results)

def randomize_array(seq_array):
    seqs=[]
    for i in range(len(seq_array)):
        seqs.append(seq_array[i][2])
    random.shuffle(seqs)
    for i in range(len(seq_array)):
        seq_array[i][2] = seqs[i]        
    return(seq_array)

#calculation of phist for each locus
def calculate_phist(seq_array):
    Nsamples = [0] * len(populations) #sample size per population
    SSDin = [0] * len(populations) #sum of squared differences wihtin populations
    SSDtot = 0 #sum of squared differences total

    missing = [] #find missing data
    for i in range(len(seq_array)):
        if seq_array[i][2] == '.' or seq_array[i][1]==-9:
            missing.append(True)
        else:
            missing.append(False)
    
    for i in range(len(seq_array)): #get sequence length
        if missing[i] == False:
            seqlen=len(seq_array[i][2])
            break

    for i in range(len(seq_array)): #sum differences within pops and total
        if missing[i]==False:
            Nsamples[int(seq_array[i][1])] += 1
            for j in range(i+1,len(seq_array)):
                if missing[j] == False:
                    diffs = sum(ch1 != ch2 for ch1,ch2 in zip(seq_array[i][2],seq_array[j][2]))
                    SSDtot += diffs
                    if seq_array[i][1]==seq_array[j][1]:
                        SSDin[int(seq_array[i][1])] += diffs
    #multiply by 2 for square (full) matrix 
    SSDtot = SSDtot*2
    SSDin=[i*2 for i in SSDin] 
                
    #calculate phi-st
    SSD_WP = sum(float(SSDin[i])/(2*Nsamples[i]) for i in range(len(populations)) if Nsamples[i]>0)
    SSD_AP = float(SSDtot)/(2*sum(Nsamples)) - SSD_WP
    Vwithin = SSD_WP/(sum(Nsamples)-sum(i > 0 for i in Nsamples))
    popsum = (sum(pow(i,2) for i in Nsamples))/sum(Nsamples) 
    pops_with_samples = sum(i > 0 for i in Nsamples)
    if pops_with_samples > 1:
        weightedN = (sum(Nsamples) - popsum)/(pops_with_samples-1)
        Vamong = (SSD_AP/(pops_with_samples-1) - Vwithin)/weightedN
        if Vamong+Vwithin > 0:
            phist = Vamong/(Vamong+Vwithin)
        else:
            phist=0
    else:
        phist=1

    nucdiv = [float(SSDin[i])/(Nsamples[i]*Nsamples[i]*seqlen) for i in range(len(populations)) if Nsamples[i]>0]
    nucdiv.append(float(SSDtot)/(sum(Nsamples)*sum(Nsamples)*seqlen))

    return(phist,nucdiv,seqlen,SSDtot,SSDin,Nsamples)

outf = open(args.o,'w')

print('Gathering info from sample info file, skipping samples where population is -9:')
infof = open(args.si,'r')
header = infof.readline()
lines = infof.readlines()
PopVector = [i.split()[4] for i in lines for x in (0,1)]
pops = list(set(PopVector))
populations = [i for i in pops if i != '-9']
populations.sort()

print('Found %.0f' % ((sum(1 for i in PopVector if i != '-9'))/2)+' samples in '+str(len(populations))+' populations')
print(populations)

#count number of clusters
print('\nCounting number of clusters')
inf = open(args.i,'r')
file = inf.read()
num_clusters = file.count('Clstr')
print('Found '+str(num_clusters)+' clusters')
inf.close()

percent10 = num_clusters*0.1
target = percent10
percent = 10

total_length = 0
GT_SSDtot = [0] * int(((len(populations) * (len(populations)-1))/2) + 1)
GT_SSDin = [0 for i in populations]
nd_tot = [[0,0] for i in range(len(populations)+1)]
pop_samples = [0 for i in populations]

outf.write('Cluster\tSeq_len\t')
for i in range(len(populations)):
    for j in range(i+1,len(populations)):
        outf.write(populations[i]+'-'+populations[j]+'\t')
outf.write('All_Pops\t')
for i in range(len(populations)):
    outf.write(populations[i]+'\t')
outf.write('All_Pops\n')

print('\nCalculating phi-st for each cluster\n')
inf = open(args.i,"r")
for x in range(num_clusters):
    if x+1 >= target:
        print(str(x)+' clusters ~ '+str(percent)+'%')
        target += percent10
        percent += 10

    line = inf.readline()
    data = line.split()
    cluster = data[1]

    inf.readline()
    locus_array = []
    for y in PopVector:
        data = inf.readline()
        data = data.split()
        locus_array.append([data[0],y,data[1]])
        
    PHIst_values=[]
    pop_pair = -1
    for i in range(len(populations)):
        for j in range(i+1,len(populations)):
            for k in range(len(locus_array)):
                if PopVector[k]==populations[i]:
                    locus_array[k][1]=0
                elif PopVector[k]==populations[j]:
                    locus_array[k][1]=1
                else:
                    locus_array[k][1]=-9
            pop_pair += 1
            PHIst,nd,length,locSSDtot,locSSDin,locus_samples = calculate_phist(locus_array)
            PHIst_values.append(PHIst) 
            GT_SSDtot[pop_pair] += locSSDtot

    for i in range(len(locus_array)):
        if PopVector[i] in populations:
            locus_array[i][1]=populations.index(PopVector[i])
        else:
            locus_array[i][1]=-9
    PHIst,nd,length,locSSDtot,locSSDin,locus_samples = calculate_phist(locus_array)
    for i in range(len(populations)):
        if locus_samples[i] > 0:
            nd_tot[i][0] += nd[i]*length
            nd_tot[i][1] += length
    nd_tot[-1][0] += nd[-1]*length
    nd_tot[-1][1] += length
    pop_samples = [i+j for i,j in zip (locus_samples,pop_samples)]
    
    PHIst_values.append(PHIst)
    GT_SSDtot[-1] += locSSDtot
    GT_SSDin = [GT_SSDin[i]+locSSDin[i] for i in range(len(locSSDin))]

#        for j in populations:
#            outf2.write(str(locus)+'\t'+str(sum(allele_counts[j-1]))+'\t'+str(len(allele_counts[j-1])))
#            for i in allele_counts[j-1]:
#                outf2.write('\t'+str(i))
#            outf2.write('\n')

    outf.write(str(cluster)+'\t'+str(length))
    for i in PHIst_values:
        outf.write('\t'+str(i))
    for i in nd:
        outf.write('\t'+str(i))
    outf.write('\n')
    total_length += length

###calculate overall values across all loci###
outf.write('ALL'+'\t'+str(total_length)+'\t')
#Nsamples = [PopVector.count(i) for i in populations]
Nsamples = [i/num_clusters for i in pop_samples]
pop_pair = -1
for i in range(len(populations)):
    for j in range(i+1,len(populations)):
        pop_pair += 1
        SSD_WP = 0
        if Nsamples[i]>0 and Nsamples[j]>0:
            SSD_WP = float(GT_SSDin[i])/(2*Nsamples[i]) + float(GT_SSDin[j])/(2*Nsamples[j])
            SSD_AP = float(GT_SSDtot[pop_pair])/(2*(Nsamples[i]+Nsamples[j])) - SSD_WP
            Vwithin = SSD_WP/(Nsamples[i]+Nsamples[j]-2)
            weightedN = (Nsamples[i]+Nsamples[j]) - (pow(Nsamples[i],2)/(Nsamples[i]+Nsamples[j]) + pow(Nsamples[j],2)/(Nsamples[i]+Nsamples[j]))
            Vamong = (SSD_AP - Vwithin)/weightedN
            if Vamong+Vwithin > 0:
                phist = Vamong/(Vamong+Vwithin)
            else:
                phist=0
        else:
            phist=1
        outf.write(str(phist)+'\t')

SSD_WP = 0
for i in range(len(populations)):
    if Nsamples[i]>0:
        SSD_WP += float(GT_SSDin[i])/(2*Nsamples[i])
SSD_AP = float(GT_SSDtot[-1])/(2*sum(Nsamples)) - SSD_WP
Vwithin = SSD_WP/(sum(Nsamples)-len(Nsamples))
popsum = sum(pow(i,2)/sum(Nsamples) for i in Nsamples)
pops_with_samples = sum(i > 0 for i in Nsamples)
if pops_with_samples > 1:
    weightedN = (sum(Nsamples) - popsum)/(pops_with_samples-1)
    Vamong = (SSD_AP/(pops_with_samples-1) - Vwithin)/weightedN
    if Vamong+Vwithin > 0:
        phist = Vamong/(Vamong+Vwithin)
    else:
        phist=0
else:
    phist=1
outf.write(str(phist))
#for i in range(len(populations)):
#    nd = float(GT_SSDin[i])/(Nsamples[i]*Nsamples[i]*total_length)
#    outf.write(str(nd)+'\t')
#nd = float(GT_SSDtot[-1])/(sum(Nsamples)*sum(Nsamples)*total_length)
#outf.write(str(nd)+'\n')
#outf.write('ALL\t')
#for i in GT_SSDtot:
#    outf.write('\t')
for i in nd_tot:
    outf.write('\t'+str(i[0]/i[1])) 
outf.write('\n')
 
print('\n\nFinished!!\n\n')
outf.close()
#outf2.close()
