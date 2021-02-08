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
        check basic case
        """
        varPos = 7
        sequence = "CAAAAATAAAAAC"
        windowsSize = 10        
        assert set_window(varPos, sequence, windowsSize) == [1,11]
    
    def test_set_window_far_left(self):
        """
        test case function set_window
        window's left coordinate is out of the sequence
        """
        varPos = 2
        sequence = "CAAAAATAAAAAC"
        windowsSize = 10
        assert set_window(varPos, sequence, windowsSize) == [0,6]

    def test_set_window_far_right(self):
        """
        test case function set_window
        window's right coordinate is out of the sequence
        """
        varPos = 10
        sequence = "CAAAAATAAAAAC"
        windowsSize = 10
        assert set_window(varPos, sequence, windowsSize) == [4,13]
    
    def test_set_window_far(self):
        """
        test case function set_window
        windowsSize bigger than the sequence
        """
        varPos = 5
        sequence = "CAAAAATAAAAAC"
        windowsSize = 100
        assert set_window(varPos, sequence, windowsSize) == [0,13]

    def test_set_window_var_out(self):
        """
        test case function set_window
        position of the variant is out of the sequence
        """
        varPos = 12
        sequence = "CAAAAATAAAAAC"
        windowsSize = 10
        assert set_window(varPos, sequence, windowsSize) == [6,13]
    
    def test_window_sequence_basecase(self):
        """
        test case function window_sequence
        check pyfaidx object is correctly processed
        """
        varPos = 7
        varChrom = "chrom1"
        fasta = pyfaidx.Fasta("tests/data/genome.fasta")
        #TATATATATATAATGCCCGATCCCTTT
        windowsSize = 10
        assert window_sequence(varPos, varChrom, fasta, windowsSize) == "ATATATATAT"
    
    def test_window_sequence_far_right(self):
        """
        test case function window_sequence
        check pyfaidx object is correctly processed when window is out of the sequence record
        """
        varPos = 25
        varChrom = "chrom1"
        fasta = pyfaidx.Fasta("tests/data/genome.fasta")
        #TATATATATATAATGCCCGATCCCTTT
        windowsSize = 10
        assert window_sequence(varPos, varChrom, fasta, windowsSize) == "ATCCCTTT"
    

if __name__ == '__main__':
    unittest.main()