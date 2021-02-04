## Standard library imports
import sys
import os
import pickle
import unittest
from unittest.mock import Mock
## Local applications import
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gffutils
import pyfaidx
from anvage.flankingSequences import set_window
from anvage.flankingSequences import window_sequence
from anvage.flankingSequences import vcf_flanking_sequences


class flankingSequencesTest(unittest.TestCase):
    """
    test case to check function from the module flankingSequences   
    """ 

    def test_set_window_basecase(self):
        """
        test case function set_window
        check 
        """
        varPos = 7
        sequence = "CAAAAATAAAAAC"
        windowsSize = 10        
        assert set_window(varPos, sequence, windowsSize) == [1,11]



if __name__ == '__main__':
    unittest.main()