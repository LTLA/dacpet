\name{compressMatrix}
\alias{compressMatrix}

\title{Compress the count matrix}
\description{Compress the matrix by adding together homo-linker counts for the same library}

\usage{
compressMatrix(data, hetero, libname)
}

\arguments{
\item{data}{a list of interaction data, like that produced by \code{\link{countPET}} or \code{\link{recountPET}}}
\item{hetero}{a logical vector indicating whether the directory represents heterolinker counts}
\item{libname}{a character vector containing the library name for each directory}
}

\details{
After processing with \code{\link{preparePET}}, each library will yield
multiple count directories corresponding to different linker combinations. The
library name in \code{libname} should be constant for all directories generated
from the same library. Directories corresponding to hetero-linkers should be
denoted as \code{TRUE} in \code{hetero}. Homo-linker counts for each library will 
then be added together and reported in \code{counts}. 

Similar addition is performed for \code{totals}. However, it is worth pointing
out that the all entries for the same library will have the same value in
\code{totals}. This is because homo- and hetero-linker counts are taken from a
single library. Computing the library size for only, say, the hetero-linker
counts would be misleading as the difference in the total number of reads is
not attributable to sequencing depth.
}

\value{
A list containing all the elements of \code{data}, with changes to the following elements:
\item{counts}{an integer count matrix where each column corresponds to a hetero-/homo-linker
count in a library, and each row corresponds to a pair of windows}
\item{totals}{an integer vector specifying the total number of tag pairs in each library, rather than in each directory}
}

\seealso{
\code{\link{countPET}},
\code{\link{recountPET}}
}

\author{Aaron Lun}

\examples{
out1 <- system.file("exdata", "li_AA.bam", package="dacpet")
stats <- preparePET(out1, "gunkAA")
out2 <- system.file("exdata", "li_BB.bam", package="dacpet")
stats <- preparePET(out2, "gunkAB")
out3 <- system.file("exdata", "li_BB.bam", package="dacpet")
stats <- preparePET(out2, "gunkBB")

status <- c("gunkAA", "gunkAB", "gunkBB")
is.hetero <- grepl("AB$", status)
lib.name <- sub("[AB][AB]$", "", status)
data <- countPET(status, width=1e6, filter=1L)

freqs <- compressMatrix(data, is.hetero, lib.name)
nrow(data$counts)
nrow(freqs$counts)
}
