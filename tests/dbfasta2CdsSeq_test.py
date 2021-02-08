## Standard library imports
import sys
import os
#import pickle
import unittest
## Third party imports (bioinformatics)
import gffutils
import pyfaidx
## Local applications import
from anvage.objets import CdsSeq
from anvage.dbfasta2CdsSeq import dbfasta2CdsSeq


class VariantCodonTest(unittest.TestCase):
    """
    test case to check function from the module VariantCodon
    """
    def setUp(self):
        dirname = os.path.dirname(__file__) 
        dbGffFile = os.path.join(dirname, 'data/gff.db')
        genomeFastaFile = os.path.join(dirname, 'data/genome.fasta')
        self.db = gffutils.FeatureDB(dbGffFile)
        self.fasta = pyfaidx.Fasta(genomeFastaFile)
        #self.cdsSeqList = pickle.load( open( cdsSeqListFile, "rb" ) )
        self.cdsSeqList = [
            CdsSeq(seqid='chrom1', start=13, end=27, sequence='ATGCCCGATCCCTTT'),
            CdsSeq(seqid= 'chrom2', start= 6, end= 14, sequence= 'ATGATGAGT'),   
            CdsSeq(seqid= 'chrom3', start= 10, end= 18, sequence= 'ATGCAGTCT'),
            CdsSeq(seqid= 'chrom3', start= 159, end= 179, sequence= 'ATGCCCCCCTTTTGTATGGCT'),
            CdsSeq(seqid= 'chrom4', start= 34, end= 45, sequence= 'ATGCCCGATCCC')
        ]

    def test_dbfasta2CdsSeq(self):
        """
        test case function dbfasta2CdsSeq
        """
        assert dbfasta2CdsSeq(self.db, self.fasta) == self.cdsSeqList

if __name__ == '__main__':
    unittest.main()