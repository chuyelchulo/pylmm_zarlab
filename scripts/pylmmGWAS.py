#!/usr/bin/python

# pylmm is a python-based linear mixed-model solver with applications to GWAS

# Copyright (C) 2013  Nicholas A. Furlotte (nick.furlotte@gmail.com)

# The program is free for academic use. Please contact Nick Furlotte
# <nick.furlotte@gmail.com> if you are interested in using the software for
#commercial purposes.

#The software must not be modified and distributed without prior
#permission of the author.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
#CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import time
import sys
import os
import numpy as np
from pylmm import input, lmm

def printOutHead(): out.write("\t".join(["SNP_ID", "BETA", "BETA_SD", "F_STAT", "P_VALUE"]) + "\n")

def outputResult(id, beta, betaSD, ts, ps):
    out.write("\t".join([str(x) for x in [id, beta, betaSD, ts, ps]]) + "\n")

def printOutHeadAnnotated(): out.write(
    "\t".join(["SNP_ID", "A1", "NMISS", "BETA", "BETA_SD", "F_STAT", "P_VALUE"]) + "\n")

def outputResultAnnotated(id, beta, betaSD, ts, ps, nmiss, annotation_dict):
    try:
        a1 = annotation_dict[id]
    except KeyError:
        a1 = 'N/A'
    out.write("\t".join([str(x) for x in [id, a1, nmiss, beta, betaSD, ts, ps]]) + "\n")

parser = input.get_gwas_parser()
options = parser.parse_args()

outFilename = options.outfileBase

# Reading Annotation File
if options.afile:
    annotation_dict = {}
    # Remove weird excel formatting
    lines = '\n'.join([line for line in open(options.afile, 'r')])
    chars = list(lines)
    weird_returns = [i for i in range(len(chars)) if chars[i] == '\r']
    for i in range(len(weird_returns)):
        chars[weird_returns[i]] = '\n'
    lines = ''.join(chars)
    # print lines[:400]
    lines = lines.split('\n')
    for line in lines:
        # print line
        line = line.strip().split('\t')
        snp_id = line[0]
        allele = line[1]
        assert len(allele) == 1
        annotation_dict[snp_id] = allele
        # print len(annotation_dict.keys())

# READING PLINK input
if options.verbose:
    sys.stderr.write("Reading SNP input...\n")
if options.bfile:
    IN = input.plink(options.bfile, type='b', phenoFile=options.phenoFile, normGenotype=options.normalizeGenotype)
elif options.tfile:
    IN = input.plink(options.tfile, type='t', phenoFile=options.phenoFile, normGenotype=options.normalizeGenotype)
#elif options.pfile: IN = input.plink(options.pfile,type='p', phenoFile=options.phenoFile,normGenotype=options.normalizeGenotype)
elif options.emmaFile:
    IN = input.plink(options.emmaFile, type='emma', phenoFile=options.phenoFile, normGenotype=options.normalizeGenotype)

if not os.path.isfile(options.phenoFile or IN.fbase + '.phenos') and not os.path.isfile(options.emmaPheno):
    raise Exception(
        "No .pheno file exist for %s.  Please provide a phenotype file using the --phenofile or --emmaPHENO argument." % (
            options.phenoFile or IN.fbase + '.phenos'))

# Read the emma phenotype file if provided.
# Format should be rows are phenotypes and columns are individuals.
if options.emmaPheno:
    f = open(options.emmaPheno, 'r')
    P = []
    for line in f:
        v = line.strip().split()
        p = []
        for x in v:
            try:
                p.append(float(x))
            except:
                p.append(np.nan)
        P.append(p)
    f.close()
    IN.phenos = np.array(P).T

# READING Covariate File
if options.covfile:
    if options.verbose:
        sys.stderr.write("Reading covariate file...\n")
    P = IN.getCovariates(options.covfile)
    if options.noMean:
        X0 = P
    else:
        X0 = np.hstack([np.ones((IN.phenos.shape[0], 1)), P])
elif options.emmaCov:
    if options.verbose:
        sys.stderr.write("Reading covariate file...\n")
    P = IN.getCovariatesEMMA(options.emmaCov)
    if options.noMean:
        X0 = P
    else:
        X0 = np.hstack([np.ones((IN.phenos.shape[0], 1)), P])
else:
    X0 = np.ones((IN.phenos.shape[0], 1))

if np.isnan(X0).sum():
    raise Exception(
        "The covariate file %s contains missing values. At this time we are not dealing with this case.  "
        "Either remove those individuals with missing values or replace them in some way.")

# READING Kinship - major bottleneck for large datasets
if options.verbose: sys.stderr.write("Reading kinship...\n")
begin = time.time()
# This method seems to be the fastest and works if you already know the size of the matrix
if options.kfile[-3:] == '.gz':
    import gzip

    f = gzip.open(options.kfile, 'r')
    F = f.read()  # might exhaust mem if the file is huge
    K = np.fromstring(F, sep=' ')  # Assume that space separated
    f.close()
else:
    K = np.fromfile(open(options.kfile, 'r'), sep=" ")
K.resize((len(IN.indivs), len(IN.indivs)))
end = time.time()
# Other slower ways
#K = np.loadtxt(options.kfile)
#K = np.genfromtxt(options.kfile)
if options.verbose: sys.stderr.write(
    "Read the %d x %d kinship matrix in %0.3fs \n" % (K.shape[0], K.shape[1], end - begin))

if options.kfile2:
    if options.verbose: sys.stderr.write("Reading second kinship...\n")
    begin = time.time()
    # This method seems to be the fastest and works if you already know the size of the matrix
    if options.kfile2[-3:] == '.gz':
        import gzip

        f = gzip.open(options.kfile2, 'r')
        F = f.read()  # might exhaust mem if the file is huge
        K2 = np.fromstring(F, sep=' ')  # Assume that space separated
        f.close()
    else:
        K2 = np.fromfile(open(options.kfile2, 'r'), sep=" ")
    K2.resize((len(IN.indivs), len(IN.indivs)))
    end = time.time()
    if options.verbose:
        sys.stderr.write("Read the %d x %d second kinship matrix in %0.3fs \n" % (K2.shape[0], K2.shape[1], end - begin))

# PROCESS the phenotype data -- Remove missing phenotype values
# Remove all individuals without full phenotypes
phenoNum = IN.phenos.shape[1]
sys.stderr.write("%d number of phenotypes read\n" % phenoNum)
X0_origin = X0
K_origin = K
if options.kfile2:
    K2_origin = K2
else:
    K2_origin = None

for i in range(phenoNum):
    X0 = X0_origin
    K = K_origin
    K2 = K2_origin
    Y = IN.phenos[:, i]

    if phenoNum == 1:
        full_outFilename = outFilename
    else:
        start, end = os.path.splitext(outFilename)
        full_outFilename = start + '_{0}'.format(i) + end
    with open(full_outFilename, 'w') as out:
        if options.afile:
            printOutHeadAnnotated()
        else:
            printOutHead()

        TS, PS, output_list = lmm.GWAS2(Y, IN, K, K2, Kva=[], Kve=[], X0=X0, REML=True, refit=False)

        #write output list, include annotation dict if annotate is true
        for id, beta, sum_sqrt_betaVar, ts, ps, nmiss in output_list:
            if options.afile:
                outputResultAnnotated(id, beta, sum_sqrt_betaVar, ts, ps, nmiss, annotation_dict=annotation_dict)
            else:
                outputResult(id, beta, sum_sqrt_betaVar, ts, ps)
