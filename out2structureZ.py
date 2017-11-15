#!/usr/bin/env python3

###############################################################################
##
## This Python (v3) script converts sequences in an out file (with select
## clusters) to the STRUCTURE format. It can extract either haplotype numbers, all
## SNPs/indels all biallelic SNPs/indels, or one randomly drawn biallelic
## SNP/indel per locus. Samples with a "-9" in the population column of
## the sample info file are skipped.
##
## In contrast to the script out2structureA.py, this script assumes that
## all loci are from the Z (sex) chromosome of birds. The sex of samples
## must be provided in the sample info file. Males are coded as diploid and
## and (heterogametic) females are coded as haploid (only major allele used.
## The script can also be used for X chromosomes in XY systems, but sex
## coding needs to be (pseudo) reversed (i.e., heterogametic sex must be
## female.
##
## This script is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
###############################################################################

import os, sys, random, argparse, subprocess
from argparse import RawTextHelpFormatter

print()

#create variables that can be entered as arguments in command line
parser = argparse.ArgumentParser(description=
                                 'This Python (v3) script converts sequences in an out file (with select\n'+
                                 'clusters) to the STRUCTURE format. It can extract either haplotype\n'+
                                 'numbers, all SNPs/indels all biallelic SNPs/indels, or one randomly drawn\n'+
                                 'biallelic SNP/indel per locus. Samples with a "-9" in the population\n'+
                                 'column of the sample info file are skipped.\n\n'+
                                 'In contrast to the script out2structureA.py, this script assumes that\n'+
                                 'all loci are from the Z (sex) chromosome of birds. The sex of samples\n'+
                                 'must be provided in the sample info file. Males are coded as diploid and\n'+
                                 'and (heterogametic) females are coded as haploid (only major allele used.\n'+
                                 'The script can also be used for X chromosomes in XY systems, but sex\n'+
                                 'coding needs to be (pseudo) reversed (i.e., heterogametic sex must be\n'+
                                 'female.\n\n'+
                                 'This script is free and distributed WITHOUT warranty; without\n'+
                                 'even the implied warranty of MERCHANTABILITY or FITNESS FOR A\n'+
                                 'PARTICULAR PURPOSE.', formatter_class=RawTextHelpFormatter)

parser.add_argument('-i', type=str, metavar='infile', required=True, help='Name of input .out file with selected clusters.')
parser.add_argument('-o', type=str, metavar='outfile', required=True, help='Name of output STRUCTURE file')
parser.add_argument('-si', type=str, metavar='infofile', required=True, help='Name of sample info file.')
parser.add_argument('-ct', type=str, metavar='character_type', required=True, help='Type of characters to use (must be HAP, ALLSNP, ALLBISNP, or 1BISNP)')
parser.add_argument('-min', type=int, metavar='min_freq', default=1, help='Minimum minor allele count when using biallelic snps/indels [1]')
parser.add_argument('-hemi', type=int, metavar='hemizygous_genotypes', default=1, help='Allow (1) or do not allow (0) hemizygous genotypes in STRUCTURE\n'+
                    'file. If 1 then major allele will be written for low depth\n'+
                    'and flagged genotypes and second allele will be scored as -9.\n'+
                    'If 0 then both alleles will be scored as -9 [default: 1]')
args = parser.parse_args()

#check ct parameter
goodct = ['HAP','ALLSNP','ALLBISNP','1BISNP']
if args.ct not in goodct:
    print('ERROR: ct parameter does not match one of four possible options!\n\n')
    quit()
    
#count number of clusters
print('\nCounting number of clusters and gathering sample names')
infile = open(args.i,'r')
file = infile.read()
num_clusters = file.count('Clstr')
print('Found '+str(num_clusters)+' clusters')
infile.close()

#establish cluster count and printing progress every 10% completed
percent10 = round(num_clusters*0.1-1,0)
cluster_count = 0
target = percent10
percent = 10

#gather sample info from popfile, skip pop=-9
print('\nGathering info from sample info file, skipping samples where population is -9:')
infofile = open(args.si,'r')
header = infofile.readline()
num_samples = 0
incl_samples = 0
samplearray = []
inclsamplearray = []
inclsamplearray.append('')
inclsamplesexarray = []
for line in infofile:
    pop = line.split() [4]
    sample = line.split() [1]
    sex = line.split() [5]
    samplearray.append(pop)
    if pop != '-9':
        num_samples += 1
        incl_samples += 1
        inclsamplearray.append(sample+'a')
        inclsamplearray.append(sample+'b')
        inclsamplesexarray.append(sex)
    else:
        num_samples += 1
infofile.close()
outputarray = []
outputarray.append(inclsamplearray)
print('Found '+str(num_samples)+' samples, of which '+str(incl_samples)+' will be included in output file')

if args.ct == 'HAP':
    print('\nGathering cluster data, creating structure file based on haplotype numbers\n')

    #define output array, append included samples
    hap_array = []
    hap_array.append(inclsamplearray)
    
    #loop for the range of the number of clusters
    infile = open(args.i,'r')
    for i in range(int(num_clusters)):
        if i == target:
            print(str(i)+' clusters ~ '+str(percent)+'%')
            target = target + percent10
            percent = percent + 10

        #define cluster number, skip line 2
        cluster = infile.readline().split() [1]
        header2 = infile.readline()
        
        #create array of cluster data for included samples
        locus_array = []
        for j in range(int(num_samples)):
            if samplearray[j] != '-9':
                data = infile.readline().split()
                locus_array.append(data)
                data = infile.readline().split()
                locus_array.append(data)
            else:
                infile.readline()
                infile.readline()

        #check if cluster is variable for included samples
        haps = []
        for x in range(len(locus_array)):
            hap = locus_array[x][3]
            if hap not in haps:
                haps.append(hap)
        
        #if cluster is variable for included samples, proceed with data extraction
        if len(haps) > 1:
    
            locus_hap = []
            locus_hap.append(cluster)

            for k in range(int(incl_samples)):
                if int(locus_array[k*2][7]) == 0:
                    locus_hap.append('-9')
                    locus_hap.append('-9')
                elif int(locus_array[k*2][7]) == 1:
                    if inclsamplesexarray[k] != 'F':
                        locus_hap.append(int(locus_array[k*2][3])+1)
                        locus_hap.append(int(locus_array[k*2+1][3])+1)
                    else:
                        locus_hap.append(int(locus_array[k*2][3])+1)
                        locus_hap.append('-9')
                else: #int(locus_array[i*2][7]) > 1:
                    if args.hemi == 0:
                        locus_hap.append('-9')
                        locus_hap.append('-9')
                    else: #args.hemi == 1
                        locus_hap.append(int(locus_array[k*2][3])+1)
                        locus_hap.append('-9')
                            
            hap_array.append(locus_hap)

    print('\nFound '+str(len(hap_array)-1)+' clusters that were variable for included samples')

    #write output STRUCTURE file
    outfile = open(args.o,'w')
    for i in range((incl_samples*2)+1):
        for j in range(len(hap_array)):
            outfile.write(str(hap_array[j][i])+'\t')
        outfile.write('\n')
    outfile.close()

    print('\nStructure file '+args.o+' created.\n\nFinished!!\n\n')

if args.ct == 'ALLSNP':
    print('\nGathering cluster data, creating structure file based on all SNPs/indels numbers\n')

    #define output array, append included samples
    allsnp_array = []
    allsnp_array.append(inclsamplearray)

    var_cluster = 0
    
    #loop for the range of the number of clusters
    infile = open(args.i,'r')
    for i in range(int(num_clusters)):
        if i == target:
            print(str(i)+' clusters ~ '+str(percent)+'%')
            target = target + percent10
            percent = percent + 10

        #define cluster number, skip line 2
        cluster = infile.readline().split() [1]
        header2 = infile.readline()
        
        #create array of cluster data for included samples
        locus_array = []
        for j in range(int(num_samples)):
            if samplearray[j] != '-9':
                data = infile.readline().split()
                locus_array.append(data)
                data = infile.readline().split()
                locus_array.append(data)
            else:
                infile.readline()
                infile.readline()

        #check if cluster is variable for included samples
        varsites_array = []
        for x in range(len(locus_array)):
            varsites = locus_array[x][2]
            if varsites != '.' and varsites not in varsites_array:
                varsites_array.append(varsites)
        
        #if cluster is variable for included samples, proceed with data extraction
        if len(varsites_array) > 1:
            var_cluster += 1
    
            #get number of snps in locus
            num_snp = 0
            for x in range(len(locus_array)):
                if locus_array[x][2] != '.':
                    for y in locus_array[x][2]:
                        num_snp += 1
                    break

            for s in range(num_snp):
                locus_allsnps = []
                locus_allsnps.append(cluster+'.'+str(s+1))
                
                for k in range(int(incl_samples)):
                    if int(locus_array[k*2][7]) == 0:
                        locus_allsnps.append('-9')
                        locus_allsnps.append('-9')
                    elif int(locus_array[k*2][7]) == 1:
                        if inclsamplesexarray[k] != 'F':
                            if locus_array[k*2][2][s] == 'A':
                                locus_allsnps.append('1')
                            elif locus_array[k*2][2][s] == 'C':
                                locus_allsnps.append('2')
                            elif locus_array[k*2][2][s] == 'G':
                                locus_allsnps.append('3')
                            elif locus_array[k*2][2][s] == 'T':
                                locus_allsnps.append('4')
                            elif locus_array[k*2][2][s] == '0':
                                locus_allsnps.append('5')
                            elif locus_array[k*2][2][s] == '1':
                                locus_allsnps.append('6')
                            else: #locus_array[k*2][2][s] == '-':
                                locus_allsnps.append('-9')

                            if locus_array[k*2+1][2][s] == 'A':
                                locus_allsnps.append('1')
                            elif locus_array[k*2+1][2][s] == 'C':
                                locus_allsnps.append('2')
                            elif locus_array[k*2+1][2][s] == 'G':
                                locus_allsnps.append('3')
                            elif locus_array[k*2+1][2][s] == 'T':
                                locus_allsnps.append('4')
                            elif locus_array[k*2+1][2][s] == '0':
                                locus_allsnps.append('5')
                            elif locus_array[k*2+1][2][s] == '1':
                                locus_allsnps.append('6')
                            else: #locus_array[k*2][2][s] == '-':
                                locus_allsnps.append('-9')
                        else:
                            if locus_array[k*2][2][s] == 'A':
                                locus_allsnps.append('1')
                            elif locus_array[k*2][2][s] == 'C':
                                locus_allsnps.append('2')
                            elif locus_array[k*2][2][s] == 'G':
                                locus_allsnps.append('3')
                            elif locus_array[k*2][2][s] == 'T':
                                locus_allsnps.append('4')
                            elif locus_array[k*2][2][s] == '0':
                                locus_allsnps.append('5')
                            elif locus_array[k*2][2][s] == '1':
                                locus_allsnps.append('6')
                            else: #locus_array[k*2][2][s] == '-':
                                locus_allsnps.append('-9')
                            locus_allsnps.append('-9')
                            
                    else: #int(locus_array[i*2][7]) > 1:
                        if args.hemi == 0:
                            locus_allsnps.append('-9')
                            locus_allsnps.append('-9')
                        else: #args.hemi == 1
                            if locus_array[k*2][2][s] == 'A':
                                locus_allsnps.append('1')
                            elif locus_array[k*2][2][s] == 'C':
                                locus_allsnps.append('2')
                            elif locus_array[k*2][2][s] == 'G':
                                locus_allsnps.append('3')
                            elif locus_array[k*2][2][s] == 'T':
                                locus_allsnps.append('4')
                            elif locus_array[k*2][2][s] == '0':
                                locus_allsnps.append('5')
                            elif locus_array[k*2][2][s] == '1':
                                locus_allsnps.append('6')
                            else: #locus_array[k*2][2][s] == '-':
                                locus_allsnps.append('-9')
                            locus_allsnps.append('-9')
                            
                allsnp_array.append(locus_allsnps)

    print('\nFound '+str(var_cluster)+' variable clusters for included samples, '+
          'which collectively contain '+str(len(allsnp_array)-1)+' SNPs/indels')

    print('\nWriting output file')
    #write output STRUCTURE file
    outfile = open(args.o,'w')
    for i in range((incl_samples*2)+1):
        for j in range(len(allsnp_array)):
            outfile.write(str(allsnp_array[j][i])+'\t')
        outfile.write('\n')
    outfile.close()

    print('\nStructure file '+args.o+' created.\n\nFinished!!\n\n')


if args.ct == 'ALLBISNP':
    print('\nGathering cluster data, creating structure file based on all biallelic SNPs/indels\n')

    #define output array, append included samples
    allbisnp_array = []
    allbisnp_array.append(inclsamplearray)

    var_cluster = 0
    
    #loop for the range of the number of clusters
    infile = open(args.i,'r')
    for i in range(int(num_clusters)):
        if i == target:
            print(str(i)+' clusters ~ '+str(percent)+'%')
            target = target + percent10
            percent = percent + 10

        #define cluster number, skip line 2
        cluster = infile.readline().split() [1]
        header2 = infile.readline()
        
        #create array of cluster data for included samples
        locus_array = []
        for j in range(int(num_samples)):
            if samplearray[j] != '-9':
                data = infile.readline().split()
                locus_array.append(data)
                data = infile.readline().split()
                locus_array.append(data)
            else:
                infile.readline()
                infile.readline()

        #check if cluster is variable for included samples
        varsites_array = []
        for x in range(len(locus_array)):
            varsites = locus_array[x][2]
            if varsites != '.' and varsites not in varsites_array:
                varsites_array.append(varsites)
        
        #if cluster is variable for included samples, proceed with data extraction
        if len(varsites_array) > 1:
            var_cluster += 1
    
            #get number of snps in locus
            num_snp = 0
            for x in range(len(locus_array)):
                if locus_array[x][2] != '.':
                    for y in locus_array[x][2]:
                        num_snp += 1
                    break

            for s in range(num_snp):
                                
                locus_allbisnps = []
                locus_allbisnps.append(cluster+'.'+str(s+1))
                
                for k in range(int(incl_samples)):
                    if int(locus_array[k*2][7]) == 0:
                        locus_allbisnps.append('-9')
                        locus_allbisnps.append('-9')
                    elif int(locus_array[k*2][7]) == 1:
                        if inclsamplesexarray[k] != 'F':
                            if locus_array[k*2][2][s] == 'A':
                                locus_allbisnps.append('1')
                            elif locus_array[k*2][2][s] == 'C':
                                locus_allbisnps.append('2')
                            elif locus_array[k*2][2][s] == 'G':
                                locus_allbisnps.append('3')
                            elif locus_array[k*2][2][s] == 'T':
                                locus_allbisnps.append('4')
                            elif locus_array[k*2][2][s] == '0':
                                locus_allbisnps.append('5')
                            elif locus_array[k*2][2][s] == '1':
                                locus_allbisnps.append('6')
                            else: #locus_array[k*2][2][s] == '-':
                                locus_allbisnps.append('-9')

                            if locus_array[k*2+1][2][s] == 'A':
                                locus_allbisnps.append('1')
                            elif locus_array[k*2+1][2][s] == 'C':
                                locus_allbisnps.append('2')
                            elif locus_array[k*2+1][2][s] == 'G':
                                locus_allbisnps.append('3')
                            elif locus_array[k*2+1][2][s] == 'T':
                                locus_allbisnps.append('4')
                            elif locus_array[k*2+1][2][s] == '0':
                                locus_allbisnps.append('5')
                            elif locus_array[k*2+1][2][s] == '1':
                                locus_allbisnps.append('6')
                            else: #locus_array[k*2][2][s] == '-':
                                locus_allbisnps.append('-9')
                        else:
                            if locus_array[k*2][2][s] == 'A':
                                locus_allbisnps.append('1')
                            elif locus_array[k*2][2][s] == 'C':
                                locus_allbisnps.append('2')
                            elif locus_array[k*2][2][s] == 'G':
                                locus_allbisnps.append('3')
                            elif locus_array[k*2][2][s] == 'T':
                                locus_allbisnps.append('4')
                            elif locus_array[k*2][2][s] == '0':
                                locus_allbisnps.append('5')
                            elif locus_array[k*2][2][s] == '1':
                                locus_allbisnps.append('6')
                            else: #locus_array[k*2][2][s] == '-':
                                locus_allbisnps.append('-9')
                            locus_allbisnps.append('-9')

                    else: #int(locus_array[i*2][7]) > 1:
                        if args.hemi == 0:
                            locus_allbisnps.append('-9')
                            locus_allbisnps.append('-9')
                        else: #args.hemi == 1
                            if locus_array[k*2][2][s] == 'A':
                                locus_allbisnps.append('1')
                            elif locus_array[k*2][2][s] == 'C':
                                locus_allbisnps.append('2')
                            elif locus_array[k*2][2][s] == 'G':
                                locus_allbisnps.append('3')
                            elif locus_array[k*2][2][s] == 'T':
                                locus_allbisnps.append('4')
                            elif locus_array[k*2][2][s] == '0':
                                locus_allbisnps.append('5')
                            elif locus_array[k*2][2][s] == '1':
                                locus_allbisnps.append('6')
                            else: #locus_array[k*2][2][s] == '-':
                                locus_allbisnps.append('-9')
                            locus_allbisnps.append('-9')

                #check if snp is biallelic and passes min freq threshold
                snp_array = []
                for x in range(1,len(locus_allbisnps)):
                    if locus_allbisnps[x] != '-9':
                        if locus_allbisnps[x] not in snp_array:
                            snp_array.append(locus_allbisnps[x])
                if len(set(snp_array)) == 2:
                    bases = list(set(snp_array))
                    if min(locus_allbisnps.count(bases[0]),locus_allbisnps.count(bases[1])) >= args.min:
                        allbisnp_array.append(locus_allbisnps)

    print('\nFound '+str(var_cluster)+' variable clusters for included samples, '+
          'which collectively contain '+str(len(allbisnp_array)-1)+' biallelic SNPs/indels')

    print('\nWriting output file')
    #write output STRUCTURE file
    outfile = open(args.o,'w')
    for i in range((incl_samples*2)+1):
        for j in range(len(allbisnp_array)):
            outfile.write(str(allbisnp_array[j][i])+'\t')
        outfile.write('\n')
    outfile.close()

    print('\nStructure file '+args.o+' created.\n\nFinished!!\n\n')

if args.ct == '1BISNP':
    print('\nGathering cluster data, creating structure file based on one biallelic SNPs/indel per cluster\n')

    #define output array, append included samples
    onebisnp_array = []
    onebisnp_array.append(inclsamplearray)

    var_cluster = 0
    
    #loop for the range of the number of clusters
    infile = open(args.i,'r')
    for i in range(int(num_clusters)):
        if i == target:
            print(str(i)+' clusters ~ '+str(percent)+'%')
            target = target + percent10
            percent = percent + 10

        #define cluster number, skip line 2
        cluster = infile.readline().split() [1]
        header2 = infile.readline()
        
        #create array of cluster data for included samples
        locus_array = []
        for j in range(int(num_samples)):
            if samplearray[j] != '-9':
                data = infile.readline().split()
                locus_array.append(data)
                data = infile.readline().split()
                locus_array.append(data)
            else:
                infile.readline()
                infile.readline()

        #check if cluster is variable for included samples
        varsites_array = []
        for x in range(len(locus_array)):
            varsites = locus_array[x][2]
            if varsites != '.' and varsites not in varsites_array:
                varsites_array.append(varsites)
        
        #if cluster is variable for included samples, proceed with data extraction
        if len(varsites_array) > 1:
            var_cluster += 1
    
            #get number of snps in locus
            num_snp = 0
            for x in range(len(locus_array)):
                if locus_array[x][2] != '.':
                    for y in locus_array[x][2]:
                        num_snp += 1
                    break

            #randomize order of snps
            sites = []
            for z in range(num_snp):
                sites.append(z)
            sites_rand = random.sample(sites,len(sites))

            for s in range(num_snp):
                                
                locus_onebisnp = []
                locus_onebisnp.append(cluster+'.'+str(sites_rand[s]+1))
                
                for k in range(int(incl_samples)):
                    if int(locus_array[k*2][7]) == 0:
                        locus_onebisnp.append('-9')
                        locus_onebisnp.append('-9')
                    elif int(locus_array[k*2][7]) == 1:
                        if inclsamplesexarray[k] != 'F':
                            if locus_array[k*2][2][sites_rand[s]] == 'A':
                                locus_onebisnp.append('1')
                            elif locus_array[k*2][2][sites_rand[s]] == 'C':
                                locus_onebisnp.append('2')
                            elif locus_array[k*2][2][sites_rand[s]] == 'G':
                                locus_onebisnp.append('3')
                            elif locus_array[k*2][2][sites_rand[s]] == 'T':
                                locus_onebisnp.append('4')
                            elif locus_array[k*2][2][sites_rand[s]] == '0':
                                locus_onebisnp.append('5')
                            elif locus_array[k*2][2][sites_rand[s]] == '1':
                                locus_onebisnp.append('6')
                            else: #locus_array[k*2][2][sites_rand[s]] == '-':
                                locus_onebisnp.append('-9')

                            if locus_array[k*2+1][2][sites_rand[s]] == 'A':
                                locus_onebisnp.append('1')
                            elif locus_array[k*2+1][2][sites_rand[s]] == 'C':
                                locus_onebisnp.append('2')
                            elif locus_array[k*2+1][2][sites_rand[s]] == 'G':
                                locus_onebisnp.append('3')
                            elif locus_array[k*2+1][2][sites_rand[s]] == 'T':
                                locus_onebisnp.append('4')
                            elif locus_array[k*2+1][2][sites_rand[s]] == '0':
                                locus_onebisnp.append('5')
                            elif locus_array[k*2+1][2][sites_rand[s]] == '1':
                                locus_onebisnp.append('6')
                            else: #locus_array[k*2][2][sites_rand[s]] == '-':
                                locus_onebisnp.append('-9')
                        else:
                            if locus_array[k*2][2][sites_rand[s]] == 'A':
                                locus_onebisnp.append('1')
                            elif locus_array[k*2][2][sites_rand[s]] == 'C':
                                locus_onebisnp.append('2')
                            elif locus_array[k*2][2][sites_rand[s]] == 'G':
                                locus_onebisnp.append('3')
                            elif locus_array[k*2][2][sites_rand[s]] == 'T':
                                locus_onebisnp.append('4')
                            elif locus_array[k*2][2][sites_rand[s]] == '0':
                                locus_onebisnp.append('5')
                            elif locus_array[k*2][2][sites_rand[s]] == '1':
                                locus_onebisnp.append('6')
                            else: #locus_array[k*2][2][sites_rand[s]] == '-':
                                locus_onebisnp.append('-9')
                            locus_onebisnp.append('-9')

                    else: #int(locus_array[i*2][7]) > 1:
                        if args.hemi == 0:
                            locus_onebisnp.append('-9')
                            locus_onebisnp.append('-9')
                        else: #args.hemi == 1
                            if locus_array[k*2][2][sites_rand[s]] == 'A':
                                locus_onebisnp.append('1')
                            elif locus_array[k*2][2][sites_rand[s]] == 'C':
                                locus_onebisnp.append('2')
                            elif locus_array[k*2][2][sites_rand[s]] == 'G':
                                locus_onebisnp.append('3')
                            elif locus_array[k*2][2][sites_rand[s]] == 'T':
                                locus_onebisnp.append('4')
                            elif locus_array[k*2][2][sites_rand[s]] == '0':
                                locus_onebisnp.append('5')
                            elif locus_array[k*2][2][sites_rand[s]] == '1':
                                locus_onebisnp.append('6')
                            else: #locus_array[k*2][2][sites_rand[s]] == '-':
                                locus_onebisnp.append('-9')
                            locus_onebisnp.append('-9')

                #check if snp is biallelic and passes min freq threshold
                snp_array = []
                for x in range(1,len(locus_onebisnp)):
                    if locus_onebisnp[x] != '-9':
                        if locus_onebisnp[x] not in snp_array:
                            snp_array.append(locus_onebisnp[x])
                if len(set(snp_array)) == 2:
                    bases = list(set(snp_array))
                    if min(locus_onebisnp.count(bases[0]),locus_onebisnp.count(bases[1])) >= args.min:
                        onebisnp_array.append(locus_onebisnp)
                        break

    print('\nFound '+str(var_cluster)+' variable clusters for included samples, '+
          'and '+str(len(onebisnp_array)-1)+' contained at least one biallelic SNP/indel')

    print('\nWriting output file')
    #write output STRUCTURE file
    outfile = open(args.o,'w')
    for i in range((incl_samples*2)+1):
        for j in range(len(onebisnp_array)):
            outfile.write(str(onebisnp_array[j][i])+'\t')
        outfile.write('\n')
    outfile.close()

    print('\nStructure file '+args.o+' created.\n\nFinished!!\n\n')
