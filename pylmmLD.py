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
basicGroup.add_option("--candidateNum", dest="candidateNums", type="int",default=100,help="You can specify the number of candidate snp to choose from to pair with a specific snp")


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
imputation = np.ones((numSNP,2))*np.nan
for snp,id  in IN:
	count +=1
	if options.verbose and count % 1000 ==0:
		sys.stderr.write("At SNP %d\n" % count)
	all_snps.append(snp)
	all_ids.append(id)
indiNum= len(all_snps[0])
if(options.verbose):
	sys.stderr.write("start to pick candidate...")


for i in range(0,numSNP, options.candidateNums/2):
	temp = correlation(all_snps[i:min(numSNP,i+options.candidateNums)])
	for j in range (len(temp)):
		temp[j][j]=0
		max_cor= max(abs(temp[j,]))
		a=list(abs(temp[j,]))
		index=a.index(max_cor)
		if np.isnan(imputation[i+j][0]):
			imputation[i+j][0]=temp[j][index]
			imputation[i+j][1]=i+index
		elif max_cor > abs(imputation[i+j][0]):
			imputation[i+j][0]=temp[j][index]
			imputation[i+j][1]=i+index
		else: continue
	if (i+len(temp) >= numSNP):
		break
count=0
for i in range(numSNP):
	count +=1
	if options.verbose and count%1000 ==0:
		sys.stderr.write("Imputation at SNP %d\n" % count)
	for j in range(indiNum):
		if np.isnan(all_snps[i][j]):	
			snp_index=int(imputation[i][1])
			all_snps[i][j]=imputation[i][0]* all_snps[snp_index][j]
all_snps=np.array(all_snps)
if(options.verbose):
	sys.stderr.write("Imputed file saved as %s type\n" % all_snps.dtype)
np.savetxt(outFile, all_snps)

	

