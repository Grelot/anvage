import gffutils
import pyfaidx

def select_annotation_type(db, fasta, selectionAnnotationType):
    """
    list of gff3 features as fasta record of selected gff3 type (e.g. mRNA)
    """
    countFeature = db.count_features_of_type(selectionAnnotationType)
    featureList = [None] * countFeature
    i = 0
    for feature in db.features_of_type(selectionAnnotationType):
        featureList[i] = feature
        i=i+1
    featureList = list(filter(None.__ne__, featureList))
    return(featureList)


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
