#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import subprocess
import os
import shutil
import sys
import filecmp

class MyTest(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.cwd = os.getcwd()
        self.testsPath = "{}/test_files".format(self.cwd)
        self.data1 = '{}/HEK.chr22.PET4.clust.txt'.format(self.testsPath)
        self.data2 = '{}/experimentalDataSet.txt'.format(self.testsPath)
        self.karytype = "{}/karyotype.hg19.txt".format(self.cwd)
        self.chr = 'chr22'

    def test_if_it_work(self):
        cmd = "{}/getSignal.py {} {} {} > /dev/null".format(self.cwd, self.karytype, self.data1, self.chr)
        return_code = subprocess.call(cmd, shell=True)
        outFile = "{}/chr22.signal.np.uint16".format(self.testsPath)
        refFile = "{}/example data/chr22.signal.np.uint16".format(self.cwd)
        assert return_code == 0, "getSignal.py returned {} exit code".format(return_code)
        assert filecmp.cmp(outFile, refFile), "Returned file is diffent than expected." 

    def test_ignore_option(self):
        cmd = "{}/getSignal.py -p notIgnore -r 100000 {} {} {} > /dev/null".format(self.cwd, self.karytype, self.data2, self.chr)
        return_code = subprocess.call(cmd, shell=True)
        outFile = "{}/notIgnore.chr22.signal.np.uint16".format(self.testsPath)
        refFile = "{}/ref1.chr22.signal.np.uint16".format(self.testsPath)
        assert return_code == 0, "getSignal.py returned {} exit code".format(return_code)
        assert filecmp.cmp(outFile, refFile), "Returned file is diffent than expected." 

        cmd = "{}/getSignal.py -i -p ignore -r 100000 {} {} {} > /dev/null".format(self.cwd, self.karytype, self.data2, self.chr)
        return_code = subprocess.call(cmd, shell=True)
        outFile = "{}/ignore.chr22.signal.np.uint16".format(self.testsPath)
        refFile = "{}/ref2.chr22.signal.np.uint16".format(self.testsPath)
        assert return_code == 0, "getSignal.py returned {} exit code".format(return_code)
        assert filecmp.cmp(outFile, refFile), "Returned file is diffent than expected." 
        
    @classmethod
    def tearDownClass(self):
        fnames = ['chr22.signal.np.uint16', 'ignore.chr22.signal.np.uint16', 'notIgnore.chr22.signal.np.uint16' ] 
        paths = ['{}/{}'.format(self.testsPath, i) for i in fnames ]
        for i in paths:
            if os.path.isfile(i): os.remove(i)

if __name__ == "__main__":
    unittest.main()
