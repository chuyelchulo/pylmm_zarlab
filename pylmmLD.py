#!/usr/bin/python

import sys
import pdb
from optparse import OptionParser,OptionGroup
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
if options.bfile: IN = input.plink(options.bfile, type='b')
elif options.tfile: IN = input.plink(options.tfile, type='t')
elif options.emmaFile: 
	if not options.numSNPs: parser.error("You must provide the number of SNPs when specifying an emma formatted file.")
	IN = input.plink(options.emmaFile,type='emma')
else: parser.error("You must provide at least one PLINK input file base (--tfile or --bfile) or an emma formatted file (--emaSNP).")

if options.emmaFile: IN.numSNPs = options.numSNPs 

IN.getSNPIterator()
print(IN.numSNPs)

numSNP = IN.numSNPs
m = 100
W = np.ones((numSNP, m))*np.nan
a = np.ones((numSNP, m))
all_snps=[]

for snp,id  in IN:
	all_snps.append(snp)
#for i in range(numSNP):


print(all_snps[0][9])
np.savetxt(outFile, a)

	

