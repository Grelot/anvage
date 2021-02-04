from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import gffutils
import pyfaidx

def set_window(varPos, sequence, windowsSize):
    """
    From variant position and windowsSize set the windows left and right position
    Correct the windows position when the windows is outside the actual genome scaffold sequence
    Return a list with left position and right position of the window
    """
    ## VCF position is 1-based while python list is 0-based 
    varPos = varPos -1
    seqLeft = 0
    seqRight = len(sequence)
    varLeft = varPos - int(windowsSize/2)
    varRight = varPos + int(windowsSize/2)
    if varLeft < seqLeft:
        varLeft = seqLeft
    if varRight > seqRight:
        varRight = seqRight
    window = [varLeft, varRight]
    return window
    
def window_sequence(varPos, varChrom, fasta, windowsSize):
    """
    From variant position onto the genome and a given window size,
    extract the genome sequences into the set windows
    return a DNA sequence as a string object
    """
    sequence = fasta[varChrom]
    window = set_window(varPos, sequence, windowsSize)
    window_sequence = str(sequence[window[0]:window[1]].seq)
    return window_sequence

def vcf_flanking_sequences(vcf_reader, fasta, windowsSize):
    """
    Browse all variants in vcf_reader object,
    return for each variant the flanking sequences onto genome
    """
    sequences=[None]*len(vcf_reader)
    i = 0
    for variant in vcf_reader:
        sequence = window_sequence(variant.POS, variant.CHROM, fasta, windowsSize)
        sequenceID = str(variant.CHROM)+":"+str(variant.POS)
        sequenceDescription = ""
        sequences[i] =SeqRecord(Seq(sequence),id=sequenceID,description=sequenceDescription)
        i = i+1
    return sequences
