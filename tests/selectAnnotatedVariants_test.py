## Standard library imports
import sys
import os
import unittest
from unittest.mock import Mock
import gffutils
import pyfaidx
## Local applications import
from anvage.selectAnnotatedVariants import select_annotation_type
from anvage.selectAnnotatedVariants import variant_position_within


class selectAnnotatedVariantsTest(unittest.TestCase):
    """
    test case to check function from the module selectAnnotatedVariants
    """
    def setUp(self):
        dirname = os.path.dirname(__file__)
        dbGffFile = os.path.join(dirname, 'data/gff.db')
        genomeFastaFile = os.path.join(dirname, 'data/genome.fasta')
        self.db = gffutils.FeatureDB(dbGffFile)
        self.fasta = pyfaidx.Fasta(genomeFastaFile)

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
    
    def test_select_annotation_type(self):
        """
        test case function select_annotation_type        
        """
        selectionAnnotationType = "CDS"
        assert len(select_annotation_type(self.db, self.fasta, selectionAnnotationType)) == 5
