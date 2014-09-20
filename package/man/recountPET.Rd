\name{recountPET}
\alias{recountPET}

\title{Count PETs into region pairs}
\description{Count PETs into pairs of user-specified regions across the interaction space}

\usage{
recountPET(files, regions, ext=1L, filter=20L, restrict=NULL)
}

\arguments{
\item{files}{a character vector containing paths to the processed tag files}
\item{regions}{a \code{GRanges} object containing the user's regions of interest}
\item{ext}{an integer scalar specifying the fragment length for read extension }
\item{filter}{an integer scalar specifying the minimum count sum for an interaction}
\item{restrict}{a character vector containing the chromosomes to be analyzed}
}
	
\details{
This function will extend each tag by \code{ext} in the direction of the tag. 
For any given pair of regions, it will count the number of pairs with one extended tag overlapping each region. 
This is performed for each processed tag file. 
Region pairs that have a sum of counts less than \code{filter} across all files will be discarded. 
Again, \code{ext} is provided here as a formality but isn't particularly useful for actual analyses; see \code{\link{countPET}} for more details.

Note that the input \code{regions} is re-ordered in the final output.
This is necessary to preserve the consistency of the anchor/target definitions.
In particular, all anchor indices are guaranteed to be greater than or equal to the target indices.
If the original sorting order is required, it can be accessed in the \code{original} field of the metadata of \code{region(x)}, where \code{x} is the \code{IList} returned by this function.
}

\value{
An \code{IList} object is returned containing counts for interactions between pairs of elements in \code{regions}.
For each pair of regions, one count is returned for each tag file.
The \code{info} slot contains a dataframe with \code{totals} and \code{files} fields, specifying the total number of PETs in and the path to each tag file, respectively.
The \code{region} slot contains a re-ordered version of \code{regions}.
}

\seealso{
\code{\link{compressMatrix}},
\code{\link{countPET}},
\code{\link{IList-class}}
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
head(counts(data))

# compressMatrix works here.
is.hetero <- grepl("AB$", status)
lib.name <- sub("[AB][AB]$", "", status)
freqs <- compressMatrix(data, is.hetero, lib.name)
head(counts(freqs))
}
