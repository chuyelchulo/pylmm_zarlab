#!/usr/bin/python

# pylmm is a python-based linear mixed-model solver with applications to GWAS

# Copyright (C) 2013  Nicholas A. Furlotte (nick.furlotte@gmail.com)

# The program is free for academic use. Please contact Nick Furlotte
# <nick.furlotte@gmail.com> if you are interested in using the software for
# commercial purposes.

# The software must not be modified and distributed without prior
# permission of the author.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import sys
import pdb

from optparse import OptionParser, OptionGroup

usage = """usage: %prog [options] --[tfile | bfile] plinkFileBase outfile (--GxE --covfile covfile)
"""

parser = OptionParser(usage=usage)

basicGroup = OptionGroup(parser, "Basic Options")
GxEGroup = OptionGroup(parser, "GxE Options")
# advancedGroup = OptionGroup(parser, "Advanced Options")

# basicGroup.add_option("--pfile", dest="pfile",
#                  help="The base for a PLINK ped file")
basicGroup.add_option("--tfile", dest="tfile",
                      help="The base for a PLINK tped file")
basicGroup.add_option("--bfile", dest="bfile",
                      help="The base for a PLINK binary ped file")
basicGroup.add_option("--emmaSNP", dest="emmaFile", default=None,
                      help="For backwards compatibility with emma, we allow for \"EMMA\" file formats.  "
                           "This is just a text file with individuals on the rows and snps on the columns.")
basicGroup.add_option("--emmaNumSNPs", dest="numSNPs", type="int", default=0,
                      help="When providing the emmaSNP file you need to specify how many snps are in the file")

basicGroup.add_option("-e", "--efile", dest="saveEig", help="Save eigendecomposition to this file.")
basicGroup.add_option("-n", default=1000, dest="computeSize", type="int",
                      help="The maximum number of SNPs to read into memory at once (default 1000).  "
                           "This is important when there is a large number of SNPs, because memory could be an issue.")

basicGroup.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default=False,
                      help="Print extra info")

GxEGroup.add_option("--gxe", "--GxE",
                    action="store_true", dest="runGxE", default=False,
                    help="Run a gene-by-environment test instead of the gene test; "
                         "the environment variable should be binary and written as "
                         "the last column of the covariate file.")

GxEGroup.add_option("--covfile", dest="covfile", default=None,
                    help="The environment filename (no header)")

parser.add_option_group(basicGroup)
parser.add_optin_group(GxEGroup)
# parser.add_option_group(advancedGroup)

(options, args) = parser.parse_args()
if len(args) != 1:
    parser.print_help()
    sys.exit()

outFile = args[0]

import sys
import os
import numpy as np
from scipy import linalg
from pylmm import input, lmm
from estimate_variance_components import make_K_GxE

if not options.tfile and not options.bfile and not options.emmaFile:
    parser.error(
        "You must provide at least one PLINK input file base (--tfile or --bfile)"
        " or an emma formatted file (--emmaSNP).")

if options.verbose:
    sys.stderr.write("Reading PLINK input...\n")
if options.bfile:
    IN = input.plink(options.bfile, type='b')
elif options.tfile:
    IN = input.plink(options.tfile, type='t')
# elif options.pfile: IN = input.plink(options.pfile,type='p')
elif options.emmaFile:
    if not options.numSNPs:
        parser.error("You must provide the number of SNPs when specifying an emma formatted file.")
    IN = input.plink(options.emmaFile, type='emma')
else:
    parser.error(
        "You must provide at least one PLINK input file base (--tfile or --bfile)"
        " or an emma formatted file (--emmaSNP).")

K_G = lmm.calculateKinshipIncremental(IN, numSNPs=options.numSNPs,
                                      computeSize=options.computeSize, center=False, missing="MAF")

if options.runGxE:
    K_G_outfile = '{}_K_G{}'.format(*os.path.splitext(outFile))
else:
    K_G_outfile = outFile
if options.verbose:
    sys.stderr.write("Saving Genetic Kinship file to %s\n" % K_G_outfile)
np.savetxt(K_G_outfile, K_G)

if options.saveEig:
    if options.verbose:
        sys.stderr.write("Obtaining Eigendecomposition of K_G\n")
    K_Gva, K_Gve = lmm.calculateEigendecomposition(K_G)
    if options.verbose:
        sys.stderr.write("Saving eigendecomposition to %s.[kva | kve]\n" % K_G_outfile)
    np.savetxt("{}.kva".format(K_G_outfile), K_Gva)
    np.savetxt("{}.kve".format(K_G_outfile), K_Gve)

if options.runGxE:
    if options.verbose:
        sys.stderr.write("Reading covariate file...\n")
    X0 = IN.getCovariates(options.covfile)
    K_GxE = make_K_GxE(K_G=K_G, env=X0[:, -1])
    K_GxE_outfile = '{}_K_GxE'.format(outFile)
    if options.verbose:
        sys.stderr.write("Saving GxE Kinship file to %s\n" % K_G_outfile)
    np.savetxt(K_GxE_outfile, K_GxE)
    if options.saveEig:
        if options.verbose:
            sys.stderr.write("Obtaining Eigendecomposition of K_GxE\n")
        K_GxEva, K_GxEve = lmm.calculateEigendecomposition(K_GxE)
        if options.verbose:
            sys.stderr.write("Saving eigendecomposition to %s.[kva | kve]\n" % K_GxE_outfile)
        np.savetxt("{}.kva".format(K_GxE_outfile), K_GxEva)
        np.savetxt("{}.kve".format(K_GxE_outfile), K_GxEve)

    components = grid_search()