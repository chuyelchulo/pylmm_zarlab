#!/usr/bin/python

import sys
import pdb
from optparse import OptionParser,OptionGroup
import numpy as np
from numpy import *
def correlation(a):
#snp is a matrix, each row a snp
	snp=np.array(a)
	snp[np.isnan(snp)]= 0
#each row already normalized	
	return cov(snp,bias=1)	

usage = """usage: %prog [options] --[tfile | bfile] plinkFileBase outfile """

parser = OptionParser(usage=usage)
basicGroup = OptionGroup(parser, "basic options")
basicGroup.add_option("--tfile", dest="tfile", 
			help="the base for a PLINK tped file")
basicGroup.add_option("--bfile", dest="bfile", 
			help="the base for a PLINK binary ped file")
basicGroup.add_option("--emmaSNP", dest="emmaFile", default=None, 
			help="\"EMMA\" file formats. This is just a text file with individuals on the rows and snps on the coumns.")
basicGroup.add_option("--emmaNumSNPs", dest="numSNPs", type="int", default=0,
			help="When providing the emmaSNP file you need to specify how many snps are in the file")
basicGroup.add_option("--neighborNum", dest="neighborNums", type="int",default=100,help="You can specify the number of neighbor snp to choose from to pair with a specific snp")
basicGroup.add_option("--candidateNum",dest="candidateNums",type="int", default=10, help="You can specify the number of candidate snp to look at to impute")

basicGroup.add_option("-v", "--verbose", action="store_true", dest="verbose",
			default=False, help="print extra info")

parser.add_option_group(basicGroup)
(options, args) = parser.parse_args()
import sys
import os
import numpy as np
from pylmm import input
if len(args) != 1:
	parser.print_help()
	sys.exit()

outFile = args[0]


if not options.tfile and not options.bfile and not options.emmaFIle:
	parser.error("You mustprovide at least one PLINK input file base (--tfile or --bfile) or an emma formatted file (--emmaSNP).")

if options.verbose: sys.stderr.write("Reading PLINK input...\n")
if options.bfile: IN = input.plink(options.bfile, type='b',impute=True)
elif options.tfile: IN = input.plink(options.tfile, type='t',impute=True)
elif options.emmaFile: 
	if not options.numSNPs: parser.error("You must provide the number of SNPs when specifying an emma formatted file.")
	IN = input.plink(options.emmaFile,type='emma',impute=True)
else: parser.error("You must provide at least one PLINK input file base (--tfile or --bfile) or an emma formatted file (--emaSNP).")

if options.emmaFile: IN.numSNPs = options.numSNPs 

IN.getSNPIterator()
print(IN.numSNPs)

numSNP = IN.numSNPs
m = 100
W = np.ones((numSNP, m))*np.nan
all_snps=[]
all_ids=[]
count=0
imputation = np.ones((numSNP,options.candidateNums))*np.nan
imputation_index=np.ones((numSNP,options.candidateNums))*np.nan


for snp,id  in IN:
	count +=1
	if options.verbose and count % 1000 ==0:
		sys.stderr.write("At SNP %d\n" % count)
	all_snps.append(snp)
	all_ids.append(id)
indiNum= len(all_snps[0])
if(options.verbose):
	sys.stderr.write("start to pick candidate...")


optimal_imputation=[[] for i in range(numSNP)]
optimal_correlation=[[] for i in range(numSNP)]
optimal_index=[[] for i in range(numSNP)]
for i in range(0,numSNP,options.neighborNums/2):
	temp = correlation(all_snps[i:min(numSNP,i+options.neighborNums)])
	for j in range (len(temp)):
		temp[j][j]=0
		if not len(optimal_imputation[i+j]):
			optimal_imputation[i+j] = list(temp[j])
		else :
			arr= list(temp[j][options.neighborNums/2:])
			optimal_imputation[i+j].extend(arr)
	if (i + len(temp)>= numSNP): break
counter=0
for i in range(numSNP):
	counter+=1
	if(options.verbose) and counter % 1000 ==0:
		sys.stderr.write("At SNP %d\n" % counter)	
	x = np.isnan(all_snps[i])
	if not len(all_snps[i][x]): continue
	if i < options.neighborNums:
		offset=0
	else: opffset= (i % options.neighborNums) * options.neighborNums
	count=0
	while len(all_snps[i][x]) and count < options.candidateNums:	
		max_cor = max(abs(np.array(optimal_imputation[i])))	
		a = list(abs(np.array(optimal_imputation[i])))
		index=a.index(max_cor)
		snp = index+offset
		all_snps[i][x] = optimal_imputation[i][index] * all_snps[snp][x]
		optimal_imputation[i][index]=0
		x = np.isnan(all_snps[i])
		count +=1
	x = np.isnan(all_snps[i])
	if  not len(all_snps[i]) and (options.verbose):
		sys.stderr.write("SNP %s is not replaced\n" % all_ids[i])
#for i in range(0,numSNP, options.neighborNums/2):
#	temp = correlation(all_snps[i:min(numSNP,i+options.neighborNums)])
#	for j in range (len(temp)):
#		temp[j][j]=0
#		counter=0
#		a=sort(abs(temp[j,]))[(len(temp)-options.candidateNums):]
#		a=list(a)
#		a.reverse()
#		counter=0
#		b=list(abs(temp[j,]))
#		for k in range(options.candidateNums):
#			index=b.index(a[counter])
#			if np.isnan(imputation[i+j][k]):
#				imputation[i+j][k]=temp[j][index]
#				imputation_index[i+j][k]=i+index
#				counter +=1
#				b[index]=0
#			elif abs(temp[j][index])> abs(imputation[i+j][k]):
#				imputation[i+j][k+1 :] = imputation[i+j][k:options.candidateNums-1]
#				imputation_index[i+j][k+1:] =imputation_index[i+j][k:options.candidateNums-1]
#
#				imputation[i+j]=temp[j][index]
#				imputation_index[i+j][k]=i+index
#				counter +=1
#				b[index]=0
#			else: continue
#		
				
#	if (i+len(temp) >= numSNP):
#		break
#count=0
#non_imputed=[]
#for i in range(numSNP):
#	count +=1
#	if options.verbose and count%1000 ==0:
#		sys.stderr.write("Imputation at SNP %d\n" % count)
#	x=np.isnan(all_snps[i])
#	counter=0
#	while len(all_snps[i][x]) and counter < options.candidateNums:
#		index=int(imputation_index[i][counter])
#		all_snps[i][x] = imputation[i][counter] * all_snps[index][x]
#		x=np.isnan(all_snps[i])
#		counter +=1
#	x=np.isnan(all_snps[i])
#	if len(all_snps[i][x]):
#		non_imputed.append(all_ids[i])
#		sys.stderr.write("snp %s is not replaced \n" % all_ids[i])
#all_snps=np.array(all_snps)
if(options.verbose):
	sys.stderr.write("Imputed file saved as %s type\n" % all_snps.dtype)
np.savetxt(outFile, all_snps)

	

