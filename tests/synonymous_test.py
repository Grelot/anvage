## Standard library imports
import sys
import os
import pickle
import unittest
from unittest.mock import Mock
## Local applications import
from vcfsynonymous.objets import CdsSeq
from vcfsynonymous.synonymous import variant_position_within
from vcfsynonymous.synonymous import is_synonymous


class VariantCodonSynonymousTest(unittest.TestCase):
    """
    test case to check function from the submodule synonymous from the module VariantCodon    
    """    

    def test_variant_position_within_var_outside(self):
        """
        test case function variant_position_within
        check variant outside intervals are not detected
        """        
        variant_within_no_cds = Mock(CHROM="chrom1", POS=2) ## variant.ID = 8        
        cdsSeq_cds = Mock(seqid = "chrom1", start=13, end=27) ## cdsSeqList[0]
        assert variant_position_within(variant_within_no_cds, cdsSeq_cds) == 0

    def test_variant_position_within_var_within(self):
        """
        test case function variant_position_within
        check variant inside intervals are detected
        """        
        variant_within_cds = Mock(CHROM="chrom1", POS=21) ## variant.ID = 1       
        cdsSeq_cds = Mock(seqid = "chrom1", start=13, end=27) ## cdsSeqList[0]
        assert variant_position_within(variant_within_cds, cdsSeq_cds) == 1

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