\name{splitLinkers}
\alias{splitLinkers}

\title{Split ChIA-PET linkers}
\description{Separate ChIA-PET read pairs into FastQ files for each linker combination.}

\usage{
splitLinkers(fastq1, fastq2, linkerA, linkerB,
	min.score=10, start=NA, prefix="out", read.prefix=NULL)
}

\arguments{
\item{fastq1}{character string specifying the FastQ file containing the first read in each pair}
\item{fastq2}{character string specifying the FastQ file containing the second read in each pair}
\item{linkerA}{character string containing the linker sequence for A}
\item{linkerB}{character string containing the linker sequence for B}
\item{min.score}{integer scalar, specifying the minimum Smith-Waterman alignment score for linker identification}
\item{start}{integer scalar, specifying the position on the read where the linker sequence starts}
\item{prefix}{charater string containing the prefix for the output files}
\item{read.prefix}{character string describing the read name to use in the output files}
}

\value{
FastQ files are constructed containing pairs of tag sequences for each linker combination.
A character vector is returned containing the paths to the output FastQ files.
}

\details{
ChIA-PET generates data in terms of read pairs.
Each read is typically composed of a tag sequence, followed by a linker sequence.
Two linkers (A or B) are used to provide a control for non-specific ligation.
The linker sequence must be identified as A or B, and then removed prior to alignment of the tag to the reference genome.

This function identifies the linker sequence from each read in a pair.
It then saves the corresponding tag sequences into the appropriate FastQ files.
In this manner, four pairs of FastQ files are produced with the specified \code{prefix}, i.e., for AA, AB, BB or other (unassigned).
These files can be used for mapping with an appropriate aligner.

Linker identification is performed by k-mer matching prior to Smith-Waterman alignment.
K-mer matching avoids the need for a full local alignment in the majority of cases with near-exact matches to the linker sequence.
Alignment score must be at least \code{min.score} for proper linker identification.
Linkers will also be unassigned if the alignment scores for \code{linkerA} and \code{linkerB} are identical.

The linker sequence is identified as that starting from \code{start} and running to the length of \code{linkerA} or \code{linkerB}.
If \code{start=NA}, the corresponding length from the 3' tail of each read is used.
Specifying an integer value for \code{start} can be useful when sequencing is performed past the ligation junction.
In all cases, both linkers should have the same length.

The read name in the output files can also be altered by setting \code{read.prefix}.
This can guarantee that the names in each pair are unique, if concatenated FastQ files are provided.
}

\author{
Aaron Lun
}

\references{
Li G et al. (2010). ChIA-PET tool for comprehensive chromatin interaction analysis with paired-end tag sequencing. \emph{Genome Biology}, 11, R22.

Fullwood MJ et al. (2009). An oestrogen-receptor-alpha-bound human chromatin interactome. \emph{Nature}, 462, 58-64.
}

