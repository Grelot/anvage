
#!/usr/bin/env python3
#===============================================================================
#INFORMATION
#===============================================================================
# Codes for vcfsynonymous 
# 
# CEFE - EPHE - RESERVEBENEFIT 2021
# Guerin Pierre-Edouard
#
# ANVAGE is a toolkit to process GFF3 VCF and FASTA files together
# ANVAGE is a python3 program.
#
# git repository : https://github.com/Grelot/anvage
#
#==============================================================================
#MODULES
#==============================================================================

## Standard library imports
import os
import sys
import argparse


## Third party imports (bioinformatics)
import gffutils
import pyfaidx
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import vcf

## Local applications import
from anvage.arguments import parse_args
from anvage.objets import CdsSeq
from anvage.dbfasta2CdsSeq import dbfasta2CdsSeq
from anvage.synonymous import variant_position_within
from anvage.synonymous import is_synonymous
from anvage.flankingSequences import vcf_flanking_sequences

def main():
    args = parse_args()    
    ## read FASTA genome file
    fasta = pyfaidx.Fasta(args.genome)
    ## read VCF file
    vcf_reader = list(vcf.Reader(open(args.vcf, 'r')))
    ## commands
    #### SYNONYMOUS COMMAND
    if args.command == 'synonymous':
        print("read VCF and detect synonymous and non-synonymous coding variants...", end="")
        dbfnFile = 'currentgff.db'
        ## read GFF3 file
        if os.path.exists(dbfnFile):
            os.remove(dbfnFile)
        db = gffutils.create_db(args.annotation, dbfn=dbfnFile)
        ## From the genome(GFF3, FASTA),
        ## extract a list of CDS (coding sequences) objects
        cdsSeqList = dbfasta2CdsSeq(db, fasta)
        ## check wether variant is within a CDS
        vcf_writer_synonymous = vcf.Writer(open(args.output_prefix+'_synonymous.vcf', 'w'), vcf.Reader(open(args.vcf, 'r')))
        vcf_writer_non_synonymous = vcf.Writer(open(args.output_prefix+'_nonsynonymous.vcf', 'w'), vcf.Reader(open(args.vcf, 'r')))
        for variant in vcf_reader:
        #print(variant.CHROM, variant.POS, variant.REF, variant.ALT[0])
            for cdsSeq in cdsSeqList:        
                if variant_position_within(variant, cdsSeq):
                    #print("cds #", i)
                    #print(variant.CHROM,variant.POS, "|", cdsSeq.seqid, cdsSeq.start, cdsSeq.end)
                    if is_synonymous(variant, cdsSeq):                    
                        vcf_writer_synonymous.write_record(variant)
                    else:
                        vcf_writer_non_synonymous.write_record(variant)               
                    break 
    #### FLANK COMMAND      
    elif args.command == 'flank':
        print("read VCF and extract flanking sequences of variants from the genome...", end= "")
        windowsSize = args.windowsSize
        sequences = vcf_flanking_sequences(vcf_reader, fasta, windowsSize)
        with open(args.output_prefix+"_flanking.fasta", "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")
    else:
        print("Au revoir !")
        sys.exit(0)    
    print("done")

if __name__ == '__main__':
    main()




    