\name{countPET}
\alias{countPET}

\title{Count PETs into windows}
\description{Count PETs into pairs of windows across the interaction space}

\usage{
countPET(files, ext=1L, shift=0L, width=5000, spacing=width, filter=20L, restrict=NULL)
}

\arguments{
\item{files}{a character vector containing paths to the processed tag files}
\item{ext}{an integer scalar specifying the fragment length for read extension }
\item{shift}{an integer scalar specifying the starting location of the first window}
\item{width}{an integer scalar denoting the width of each window}
\item{spacing}{an integer scalar specifying the distance between adjacent windows}
\item{filter}{an integer scalar specifying the minimum count sum for an interaction}
\item{restrict}{a character vector containing the chromosomes to be analyzed}
}
	
\details{
Consider a set of sliding windows of width equal to \code{width} across the
genome, starting from \code{1-shift} and hopping along by \code{spacing}.  For
each pair of windows, a PET is counted if each extended tag overlaps one of the
windows. This is performed for all pairs of windows in the genome, and for each
supplied tag file. Only those window pairs with a sum of counts above
\code{filter} are retained, to reduce memory usage.

Setting a large \code{width} will increase the size of the collection area in
the interaction space. This is useful for increasing the count size.  However,
spatial resolution will be lost as adjacent interactions cannot be distinguished.
A size of 2 - 10 kbp is usually recommended. Setting \code{spacing} as a
proportion of \code{width} is also wise as it will improve computational
efficiency without further loss of resolution.

The value of \code{ext} is provided here for purposes of flexibility. Each read
will be extended to the length of \code{ext} prior to computing overlaps.
Practical use is limited as ChIP-seq-like bimodality is rarely observed.
Instead, binding sites in ChIA-PET will be marked by outward-facing read pairs
that will be captured with sufficiently large windows.

By default, \code{spacing=width} which means that the windows represent
contiguous non-overlapping bins of size equal to \code{width}. Bins are often
easier to interpret than windows as the absence of overlaps ensures that
results are non-redundant. Users should set \code{filter=1L} for a true
binning strategy whereby all non-empty bins are reported.
}

\section{Explaining the output}{
The output in \code{pairs} warrants some further explanation. The underlying
algorithm will report groups of adjacent window pairs with the same set of
counts across all libraries. So, the \code{pairs} dataframe will contain the
indices of the first and last windows of the group on each chromosome, i.e.,
\code{anchor.start} to \code{anchor.end} for the anchor chromosome, and
\code{target.start} to \code{target.end} for the target chromosome.

The relevant window pairs are represented by each combination of anchor and
target indices within the specified ranges. For example, a group with anchor
indices ranging from 10 to 12 and targets ranging from 1:5 will contain window
pairs corresponding to in \code{expand.grid(10:12, 1:5)}. Each index refers to
a range in \code{region}.  Similarly, each row in \code{pairs} will also have
an index for the row in \code{counts}, to indicate which combination of counts
corresponds to the current (set of) window pairs.

This structure can be unpacked into a format where each row of a matrix
corresponds to a single pair of windows, rather than the set of window pairs.
This is easier to manipulate during the statistical analysis. Unpacking can
be performed using the \code{\link{compressMatrix}} function.
}

\value{
A list containing:
\item{counts}{an integer matrix with one column per input file, specifying the count for each window pair}
\item{totals}{an integer vector specifying the total number of PETs in each directory}
\item{pairs}{a dataframe with integer fields, containing the start/end indices and the index of the count matrix for each group of window pairs}
\item{region}{a \code{GRanges} object containing the genomic coordinates of each window}
}

\seealso{
\code{\link{compressMatrix}},
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

# Filter is set too high for all of these, and they're shown here only for fun.
data <- countPET(c("gunkAA", "gunkAB", "gunkBB"))
data <- countPET(c("gunkAA", "gunkAB", "gunkBB"), width=10000)
data <- countPET(c("gunkAA", "gunkAB", "gunkBB"), width=10000, spacing=5000)
data <- countPET(c("gunkAA", "gunkAB", "gunkBB"), width=10000, spacing=5000, restrict=c("chr1", "chr2"))

# No filter, so things show up.
data <- countPET(c("gunkAA", "gunkAB", "gunkBB"), width=10000) 
}
