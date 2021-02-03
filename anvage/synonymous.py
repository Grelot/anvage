from Bio.Seq import Seq
from anvage.objets import CdsSeq


def variant_position_within(coordsVar, coordsInterval):
    """
    check if coordsVars is within coordsInterval. Return 0
    """
    if coordsVar.CHROM == coordsInterval.seqid:
        if coordsVar.POS >= coordsInterval.start:            
            if coordsVar.POS <= coordsInterval.end:
                return(1)
            else:
                return(0)
        else:
            return(0)
    return(0)


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