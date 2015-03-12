import unittest
import os
import subprocess

class test_pylmmGWAS(unittest.TestCase):

    def setUp(self):
        self._scriptName = "scripts/pylmmGWAS.py"
        self._ansGwasFile = 'data/snps.132k.clean.noX.pylmm.GWAS'
        self._bfile = 'data/snps.132k.clean.noX'
        self._kfile = 'data/snps.132k.clean.noX.pylmm.kin'
        self._pfile = 'data/snps.132k.clean.noX.fake.phenos'
        self._outputFile = 'data/pylmmKinshipTestOutput'
        self._computeSize = "1000"

    def tearDown(self):
        if os.path.isfile(self._outputFile):
            os.remove(self._outputFile)

    def test_pylmmGWASScript(self):
        subprocess.call(["python", self._scriptName, "-v", "--bfile", self._bfile, "--kfile", self._kfile, "--phenofile", self._pfile, self._outputFile])

        with (open(self._ansGwasFile, 'r')) as ansKey:
            with (open(self._outputFile, 'r')) as ansFile:
                line1 = ansKey.read().split()
                line2 = ansFile.read().split()
                for i in range(len(line1)):
                    self.assertTrue(line1[i].strip() == line2[i].strip(),
                        "An entry in the stored GWAS results differs from the entry in the results obtained from this test run.")