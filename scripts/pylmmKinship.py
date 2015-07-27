#!/usr/bin/python

# pylmm is a python-based linear mixed-model solver with applications to GWAS

# Copyright (C) 2013  Nicholas A. Furlotte (nick.furlotte@gmail.com)

#The program is free for academic use. Please contact Nick Furlotte
#<nick.furlotte@gmail.com> if you are interested in using the software for
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



import sys
import numpy as np
from pylmm import input, lmm

parser = input.get_kinship_parser()
options = parser.parse_args()

outFile = options.outfileBase

if options.verbose:
    sys.stderr.write("Reading PLINK input...\n")
if options.bfile:
    IN = input.plink(options.bfile,type='b')
elif options.tfile:
    IN = input.plink(options.tfile,type='t')
#elif options.pfile: IN = input.plink(options.pfile,type='p')
elif options.emmaFile: 
   IN = input.plink(options.emmaFile,type='emma')

K = lmm.calculateKinshipIncremental(IN, numSNPs=options.numSNPs, computeSize=options.computeSize, center=False, missing="MAF")

if options.verbose:
    sys.stderr.write("Saving Kinship file to %s\n" % outFile)
np.savetxt(outFile,K)

if options.saveEig:
   if options.verbose:
       sys.stderr.write("Obtaining Eigendecomposition\n")
   Kva,Kve = lmm.calculateEigendecomposition(K)
   if options.verbose:
       sys.stderr.write("Saving eigendecomposition to %s.[kva | kve]\n" % outFile)
   np.savetxt(outFile+".kva",Kva)
   np.savetxt(outFile+".kve",Kve)