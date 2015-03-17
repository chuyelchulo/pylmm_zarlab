import unittest
import math
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
        snps = np.load(self._liverSNPFile).T
        K = lmm.calculateKinship(snps)

        ansK = np.load(self._liverKinshipMatrix)

        self.assertTrue(np.allclose(K, ansK), "The kinship matrix generated during this test run " +
                        "did not match the kinship matrix stored in " + self._liverKinshipMatrix)

    def test_GWAS(self):
        Y = np.genfromtxt(self._liverPhenos)

        # Loading npdump and first 1000 snps for speed
        K = np.load(self._liverKinshipMatrix)

        snps = np.load(self._liver1000SNPFile).T
        vars = np.nanvar(snps, axis=0) #variances across the rows ignoring NaN, used to check which SNPs were not polymorphic across the given individuals

        TS,PS = lmm.GWAS(Y,snps,K,REML=True,refit=True)

        #SNPs that are not polymorphic (in the given individuals being tested) will have variance 0, this check ensures
        #that only these SNPs have a return value of NaN
        for i in range(len(PS)):
           self.assertTrue( not math.isnan(PS[i]) or vars[i] == 0, "NaN found in results corresponding to polymorphic SNP")

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