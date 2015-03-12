import unittest
import sys
import numpy as np
from scipy import linalg
from pylmm import input
from pylmm import lmm

class test_lmm(unittest.TestCase):

    def setUp(self):
        self._liverSNPFile = 'data/hmdp.liver.snps.npdump'
        self._liver1000SNPFile = 'data/hmdp.liver.snps.1000.npdump'
        self._liverKinshipMatrix = 'data/hmdp.liver.K.npdump'
        self._liverPhenos = 'data/hmdp.liver.exprs.1'
        self._liverTestFile = 'data/hmdp.liver.gwas.test'

    def tearDown(self):
        pass
        #if os.path.isfile(self._outputFile):
        #    os.remove(self._outputFile)

    def test_calculateKinship(self):
        snps = np.load('data/hmdp.liver.snps.npdump').T
        K = lmm.calculateKinship(snps)

        ansK = np.load('data/hmdp.liver.K.npdump')

        self.assertTrue(np.allclose(K, ansK), "The kinship matrix generated during this test run " +
                        "did not match the kinship matrix stored in " + self._liverKinshipMatrix)

    def test_GWAS(self):
        Y = np.genfromtxt(self._liverPhenos)

        # Loading npdump and first 1000 snps for speed
        K = np.load(self._liverKinshipMatrix)

        snps = np.load(self._liver1000SNPFile).T

        # Genome-wide scan over the 1000 SNPs.
        # This call will use REML (REML = False means use ML).
        # It will also refit the variance components for each SNP.
        # Setting refit = False will cause the program to fit the model once
        # and hold those variance component estimates for each SNP.

        TS,PS = lmm.GWAS(Y,snps,K,REML=True,refit=True)
        results = np.array([TS,PS])
        ansKey = np.load(self._liverTestFile)

        #these results include np.nan values, so allclose cannot be used, also the results are similar with each
        #run, but do vary, so we can only check for similarity to a precision of about 1e-06
        for i in range(results.shape[0]):
            for j in range(results.shape[1]):
                a = results[i,j]
                b = ansKey[i,j]
                self.assertTrue( (np.isnan(a) and np.isnan(b)) or abs(a - b) < 1e-06 ,
                                 "Mismatch on values: " + str(a) + " and " + str(b))