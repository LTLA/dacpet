\name{recountPET}
\alias{recountPET}

\title{Count PETs into regions}
\description{Count PETs into pairs of user-specified regions across the interaction space}

\usage{
recountPET(files, regions, ext=1L, filter=20L)
}

\arguments{
\item{files}{a character vector containing paths to the processed tag files}
\item{regions}{a \code{GRanges} object containing the user's regions of interest}
\item{ext}{an integer scalar specifying the fragment length for read extension }
\item{filter}{an integer scalar specifying the minimum count sum for an interaction}
}
	
\details{
This function will extend each tag by \code{ext} in the direction of the tag.
For any given pair of regions, it will count the number of extended tags
that overlap each region. This is performed for each processed tag
file. Region pairs that have a sum of counts less than \code{filter}
across all files will be discarded. Again, \code{ext} is provided here as
a formality but isn't particularly useful for actual analyses, see \code{\link{countPET}}.
}

\value{
A list containing:
\item{counts}{an integer matrix with one column per tag file, specifying the counts for each region pair}
\item{pairs}{a dataframe with integer fields, containing the indices of the interacting regions in each pair}
\item{totals}{an integer vector specifying the total number of PETs in each tag file}
\item{region}{a \code{GRanges} object containing the genomic coordinates of each region}
\item{files}{a character vector containing the paths to the processed tag files}
}

\seealso{
\code{\link{compressMatrix}},
\code{\link{countPET}}
}

\author{Aaron Lun}

\references{
Li G et al. (2010). ChIA-PET tool for comprehensive chromatin interaction analysis with paired-end tag sequencing.
	\emph{Genome Biology}, 11, R22.
}

\examples{
out1 <- system.file("exdata", "li_AA.bam", package="dacpet")
stats <- preparePET(out1, "gunkAA")
out2 <- system.file("exdata", "li_BB.bam", package="dacpet")
stats <- preparePET(out2, "gunkAB")
out3 <- system.file("exdata", "li_BB.bam", package="dacpet")
stats <- preparePET(out2, "gunkBB")

# A bit of a silly example, just for demonstration.
status <- c("gunkAA", "gunkAB", "gunkBB")
myregs <- GRanges(c("chr1", "chr2", "chr3"),
	IRanges(c(1e7, 2e7, 3e7), c(12e7, 13e7, 14e7)))
data <-recountPET(status, myregs, filter=0L)
head(data$counts)

# compressMatrix works here.
is.hetero <- grepl("AB$", status)
lib.name <- sub("[AB][AB]$", "", status)
freqs <- compressMatrix(data, is.hetero, lib.name)
head(freqs$counts)
}
