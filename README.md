dacpet
======

An R package for ChIA-PET data analysis

This package provides methods for linker splitting, tag pair counting, normalization and detection of specific interactions for ChIA-PET data. It is based primarily on the statistical methods in the edgeR package. Potential users are directed to peruse the user's guide (found in `package/inst/doc`) to implement their own analyses, and to understand some of the theory behind the pipeline.

To install this package, make sure that the latest version of R is installed. Then, at the R prompt, type:
```R
source("http://bioconductor.org/biocLite.R")
useDevel()
biocLite(c('edgeR', 'Rsamtools', 'GenomicRanges', 'rhdf5', 'csaw'))
devtools::install_github("LTLA/dacpet/package")
```
