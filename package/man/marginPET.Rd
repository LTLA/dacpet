\name{marginPET}
\alias{marginPET}

\title{Count PETs into regions}
\description{Get the marginal count of PETs in each region}

\usage{
marginPET(files, regions, ext=1L, restrict=NULL)
}

\arguments{
\item{files}{a character vector containing paths to the processed tag files}
\item{regions}{a \code{GRanges} object containing the user's regions of interest}
\item{ext}{an integer scalar specifying the fragment length for read extension }
\item{restrict}{a character vector containing the chromosomes to be analyzed}
}
	
\details{
This function will extend each tag by \code{ext} in the direction of the tag.
For each region, it will count the number of pairs with at least one extended tag overlapping the region.
This is defined as the marginal count and is a useful measure of overall coverage in a ChIA-PET experiment.
Counting is performed in this manner for each processed tag file. 
Again, \code{ext} is provided here as a formality but isn't particularly useful for actual analyses, see \code{\link{countPET}}.
}

\value{
A list containing:
\item{counts}{an integer matrix with one column per tag file, specifying the marginal counts for each region}
\item{totals}{an integer vector specifying the total number of PETs in each tag file}
\item{files}{a character vector containing the paths to the processed tag files}
}

\seealso{
\code{\link{compressMatrix}},
\code{\link{recountPET}}
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
data <- marginPET(status, myregs, filter=0L)
head(data$counts)

# compressMatrix works here.
is.hetero <- grepl("AB$", status)
lib.name <- sub("[AB][AB]$", "", status)
freqs <- compressMatrix(data, is.hetero, lib.name)
head(freqs$counts)
}
