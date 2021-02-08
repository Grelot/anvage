import argparse
import sys



HELPER_TEXT ="""
      _       ____  _____  ____   ____  _        ______  ________  
     / \     |_   \|_   _||_  _| |_  _|/ \     .' ___  ||_   __  | 
    / _ \      |   \ | |    \ \   / / / _ \   / .'   \_|  | |_ \_| 
   / ___ \     | |\ \| |     \ \ / / / ___ \  | |   ____  |  _| _  
 _/ /   \ \_  _| |_\   |_     \ ' /_/ /   \ \_\ `.___]  |_| |__/ | 
|____| |____||_____|\____|     \_/|____| |____|`._____.'|________| 
                                                                    
Pierre-Edouard GUERIN, Emilie Boulanger, Laura Benestan, Stephanie MANEL
CNRS, EPHE, Sorbonne University, Montpellier, France
Founded by biodiversa RESERVEBENEFIT 2017-2021
version 0.1 "Empty Coffee" february 2021
Usage:
> python3 anvage [options]
For help:
> python3 anvage --help
"""

def parse_args(usage=HELPER_TEXT):
    parser = argparse.ArgumentParser(description='anvage - toolkit to process gff3, VCF and fasta genome files.')

    subprasers = parser.add_subparsers(dest='command')

    select = subprasers.add_parser('select', help='extract genetic variants located on the selected type of annotated regions')
    select.add_argument("-f","--vcf", type=str, help='path of the variant Calling File (VCF) with variants you want to select', required=True)
    select.add_argument("-g","--genome", type=str, help='path of the genome sequences FASTA file', required=True)
    select.add_argument("-a","--annotation",type=str, help='path of the genome annotation GFF3 file', required=True)
    select.add_argument("-s","--selectionAnnotationType", type=str, help='GFF3 feature type column used to select a variant (e.g. "exon")', required=True)
    select.add_argument("-o","--output_prefix", type=str, help='prefix of the output VCF such as [PREFIX]_[selectionAnnotationType].vcf')

    synonymous = subprasers.add_parser('synonymous', help='detect synonymous genetic variants')
    synonymous.add_argument("-o","--output_prefix", type=str, help='prefix of the two output VCF files such as [PREFIX]_synonymous.vcf and [PREFIX]_nonsynonymous.vcf')
    synonymous.add_argument("-f","--vcf", type=str, help='path of the variant Calling File (VCF) with variants you want to determine synonymous or non-synonymous', required=True)
    synonymous.add_argument("-g","--genome", type=str, help='path of the genome sequences FASTA file', required=True)
    synonymous.add_argument("-a","--annotation",type=str, help='path of the genome annotation GFF3 file', required=True)

    flank = subprasers.add_parser('flank', help='extract genetic variants flanking sequences')
    flank.add_argument("-f","--vcf", type=str, help='path of the variant Calling File (VCF) with variants from which you want to extract genome sequences')
    flank.add_argument("-g","--genome", type=str, help='path of the genome sequences FASTA file')
    flank.add_argument("-w","--windowsSize", type=int, help='size of window centered on position of variant to extract sequences from the genome')
    flank.add_argument("-o","--output_prefix", type=str, help='prefix of the output FASTA file such as [PREFIX]_.flanking.fasta')

   
    args = parser.parse_args()

    if args.command not in ['select', 'synonymous', 'flank']:
        print(usage)
        sys.exit(0)

    return args