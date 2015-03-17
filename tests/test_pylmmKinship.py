import unittest
import os
import subprocess
import numpy as np

class test_pylmmKinship(unittest.TestCase):

    def setUp(self):
        self._scriptName = "scripts/pylmmKinship.py"
        self._ansKinshipFile = 'data/snps.132k.clean.noX.pylmm.kin'
        self._bfile = 'data/snps.132k.clean.noX'
        self._outputFile = 'data/pylmmKinshipTestOutput'
        self._computeSize1 = "1000"
        self._computeSize2 = "2000"
        self._numInd = 1219

    def tearDown(self):
        if os.path.isfile(self._outputFile):
            os.remove(self._outputFile)

    def test_pylmmKinshipScript1(self):
        subprocess.call(["python", self._scriptName,"-v", "-n", self._computeSize1, "--bfile", self._bfile, self._outputFile])

        K = np.fromfile(open(self._outputFile, 'r'), sep=" ")
        K.resize((self._numInd, self._numInd))

        ansK = np.fromfile(open(self._ansKinshipFile, 'r'), sep=" ")
        ansK.resize((self._numInd, self._numInd))

        self.assertTrue(np.allclose(K, ansK), "The kinship matrix generated during this test run " +
                        "did not match the kinship matrix stored in " + self._ansKinshipFile)

    #same test as above but using a different value for computeSize, this is done to ensure that changing the number
    #of SNPs processed as a single chunk does not affect the end results of the kinship matrix
    def test_pylmmKinshipScript2(self):
        subprocess.call(["python", self._scriptName,"-v", "-n", self._computeSize2, "--bfile", self._bfile, self._outputFile])

        K = np.fromfile(open(self._outputFile, 'r'), sep=" ")
        K.resize((self._numInd, self._numInd))

        ansK = np.fromfile(open(self._ansKinshipFile, 'r'), sep=" ")
        ansK.resize((self._numInd, self._numInd))

        self.assertTrue(np.allclose(K, ansK), "The kinship matrix generated during this test run " +
                        "did not match the kinship matrix stored in " + self._ansKinshipFile)