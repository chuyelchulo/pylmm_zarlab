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
        self._outputFile = 'data/pylmmGWASTestOutput'
        self._computeSize = "1000"

    def tearDown(self):
        if os.path.isfile(self._outputFile):
            os.remove(self._outputFile)

    def test_pylmmGWASScript(self):
        subprocess.call(["python", self._scriptName, "-v", "--bfile", self._bfile, "--kfile", self._kfile, "--phenofile", self._pfile, self._outputFile])

        with (open(self._ansGwasFile, 'r')) as ansKey:
            with (open(self._outputFile, 'r')) as ansFile:
                file1Contents = ansKey.read().split()
                file2Contents = ansFile.read().split()
                for i in range(len(file1Contents)):
                    item1 = file1Contents[i].strip()
                    item2 = file2Contents[i].strip()

                    #for numeric values, some small variations may occur in the data between runs of GWAS, so inexact
                    #matching is allowed
                    try:
                        item1 = float(item1)
                        item2 = float(item2)
                        self.assertTrue(abs(item1 - item2) < 1e-6,
                            "An entry in the stored GWAS results differs from the entry in the results obtained from this test run.")
                    #for strings, an exact match is made
                    except ValueError:
                        self.assertTrue(item1 == item2,
                            "An entry in the stored GWAS results differs from the entry in the results obtained from this test run.")