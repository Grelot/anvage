## Standard library imports
import sys
import os
import pickle
import unittest
from unittest.mock import Mock
## Local applications import
from anvage.objets import CdsSeq
from anvage.synonymous import is_synonymous


class synonymousTest(unittest.TestCase):
    """
    test case to check function from the module synonymous
    """

    def test_is_synonymous_nonsynonymous_first_codon(self):
        """
        test case function is_synonymous
        check non synonymous mutation at first codon-position
        """       
        non_synonymous_codon_first_pos = Mock(CHROM = "chrom3", POS=174, REF='G', ALT=['A']) ## variant.ID = 6        
        cdsSeq_non_synonymous_codon_first_pos = Mock(seqid="chrom3", start=159, end=179, \
            codons = ['ATG', 'CCC', 'CCC', 'TTT', 'TGT', 'ATG', 'GCC']) ## cdsSeqList[3]
        assert is_synonymous(non_synonymous_codon_first_pos, \
            cdsSeq_non_synonymous_codon_first_pos ) == 0       
   
    def test_is_synonymous_nonsynonymous_second_codon(self):
        """
        test case function is_synonymous
        check non synonymous mutation at second codon-position
        """
        non_synonymous_codon_second_pos = Mock(CHROM = "chrom3", POS=17, REF='G', ALT=['C']) ## variant.ID = 4
        cdsSeq_non_synonymous_codon_second_pos = Mock(seqid="chrom3", start=10, end=18, \
            codons = ['ATG', 'CAG', 'TCT']) ## cdsSeqList[2]
        assert is_synonymous(non_synonymous_codon_second_pos, \
            cdsSeq_non_synonymous_codon_second_pos ) == 0 
    
    def test_is_synonymous_nonsynonymous_third_codon(self):
        """
        test case function is_synonymous
        check non synonymous mutation at third cordon-position
        """
        non_synonymous_codon_third_pos = Mock(CHROM = "chrom1", POS=21, REF='T', ALT=['G']) ## variant.ID = 1
        cdsSeq_non_synonymous_codon_third_pos = Mock(seqid="chrom1", start=13, end=27, \
            codons =['ATG', 'CCC', 'GAT', 'CCC', 'TTT'] ) ## cdsSeqList[0]
        assert is_synonymous(non_synonymous_codon_third_pos, \
            cdsSeq_non_synonymous_codon_third_pos ) == 0

    def test_is_synonymous_synonymous(self):
        """
        test case function is_synonymous
        check synonymous mutation
        """        
        synonymous_var = Mock(CHROM = "chrom1", POS=27, REF='T', ALT=['C']) ## variant.ID = 2
        cdsSeq_synonymous_var = Mock(seqid="chrom1", start=13, end=27, \
            codons =['ATG', 'CCC', 'GAT', 'CCC', 'TTT'] ) ## cdsSeqList[0]
        assert is_synonymous(synonymous_var, \
            cdsSeq_synonymous_var ) == 1


if __name__ == '__main__':
    unittest.main()