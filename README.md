# ANVAGE

[![Build Status](https://travis-ci.com/Grelot/anvage.svg?branch=main)](https://travis-ci.com/Grelot/anvage)
[![PyPI version](https://badge.fury.io/py/anvage.svg)](https://badge.fury.io/py/anvage)
[![PyPi Python Versions](https://img.shields.io/pypi/pyversions/anvage.svg)](https://pypi.org/project/anvage)




**ANnotation Variants GEnome** is a toolkit software to perform routine operations that involve VCF, GFF3 and FASTA genome files.

# Installation

```
pip3 install anvage
```


# Quick start


```
anvage select --vcf tests/data/sample.vcf \
--genome tests/data/genome.fasta \
--annotation tests/data/genome.gff3 \
--selectionAnnotationType CDS \
--output_prefix res
```


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

# Credits

**ANVAGE** was coded and written by Emilie Boulanger and Pierre-Edouard Guerin.

We thank the following people for their help in the development of this software: Stephanie Manel, Laura Benestan.

The [**CNRS**](https://www.cnrs.fr/) research unit [**CEFE**](https://www.cefe.cnrs.fr/fr/), team [Biogeography of Vertebrates](https://twitter.com/ephe_bev) supports and contributes to the **ANVAGE** project.

# Contributions and Support

🐛 If you are sure you have found a bug, please submit a bug report. You can submit your bug reports on 


# Citations

You can cite the **ANVAGE** publication as follows:

> **Climate differently influences the genomic patterns of two sympatric marine fish species**
>
> *Emilie Boulanger, Laura Benestan, Pierre-Edouard Guerin, Alicia Dalongeville, David Mouillot, Stéphanie Manel*
>
> Journal of Animal Ecology. 2021 Oct 31. DOI: [10.1111/1365-2656.13623](https://doi.org/10.1111/1365-2656.13623)


# License

**ANVAGE** is licensed under MIT License. Please see the [LICENSE file](https://github.com/Grelot/anvage/blob/main/LICENSE).

