from Bio.Seq import Seq
from anvage.objets import CdsSeq


def is_synonymous(variant, cds):
    """
    check if a variant REF/ALT is synonymous or not from CDS codon in within it is located
    """
    posCodon = int((variant.POS - cds.start)/3) 
    codonVariant = cds.codons[posCodon]
    posVariantCodon = (variant.POS - cds.start)%3
    codonRef = list(codonVariant)
    codonAlt = list(codonVariant)
    codonRef[posVariantCodon] = str(variant.REF)
    codonAlt[posVariantCodon] = str(variant.ALT[0])
    codonRefSeq = Seq(''.join(codonRef))
    codonAltSeq = Seq(''.join(codonAlt))
    aaRef = str(codonRefSeq.translate(table=1))
    aaAlt = str(codonAltSeq.translate(table=1))
    #print(codonVariant, "|||", codonRef, aaRef,"|", codonAlt, aaAlt)
    if aaRef != aaAlt:
        #print("NONSYNONYME")
        return(0)
    else:
        #print("SYNONYME")
        return(1)