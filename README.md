# ANVAGE

[![Build Status](https://travis-ci.com/Grelot/anvage.svg?branch=main)](https://travis-ci.com/Grelot/anvage)

ANnotation Variants GEnome is a toolkit software to perform frequent processing on VCF, GFF3 and FASTA genome


# usage


```
anvage synonymous --vcf tests/data/sample.vcf \
--genome tests/data/genome.fasta \
--annotation tests/data/genome.gff3 \
--output_prefix res
```


```
anvage flank --vcf tests/data/sample.vcf \
--genome tests/data/genome.fasta \
--windowsSize 25 \
--output_prefix res
```


https://github.com/laurabenestan/BLAST/blob/master/flanking_sequence.py