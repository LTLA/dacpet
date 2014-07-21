\name{plotChIA}
\alias{plotChIA}

\title{Plot tag pairs}
\description{Make a plaid plot of tag pairs from a processed tag file}

\usage{
plotChIA(file, anchor, anchor.lim=NULL, target=anchor, target.lim=anchor.lim,
	width=1000, cap=50, xlab=NULL, ylab=NULL, col="red", diag=TRUE, ...)
}

\arguments{
\item{file}{a character string containing a path to a processed tag file}
\item{anchor}{a \code{GRanges} object specifying the anchor region, or a character string specifying the anchor chromosome}
\item{anchor.lim}{a numeric vector of length 2, specifying the coordinates of the anchor region}
\item{target}{a \code{GRanges} object specifying the target region, or a character string specifying the target chromosome}
\item{target.lim}{a numeric vector of length 2, specifying the coordinates of the target region}
\item{width}{an integer scalar denoting the width of each box in the plot}
\item{cap}{an integer scalar specifying the count in each box at which the colour is saturated}
\item{xlab}{a character string specifying the x-axis label, set to the anchor chromosome name by default}
\item{ylab}{a character string specifying the y-axis label, set to the target chromosome name by default}
\item{col}{a R colour indicating the transparency series to use in the plot}
\item{diag}{a logical scalar specifying whether read pairs past the diagonal should be shown in intra-chromosomal plots}
\item{...}{other arguments to be passed to \code{\link{plot}}}
}

\value{
A \code{NULL} object will be invisibly returned.
}

\details{
A plaid plot is generated similar to that originally described by
Lieberman-Aiden \emph{et al.} in their 2009 paper on Hi-C. Specifically, each
box on a plaid plot represents a pair of bins of size equal to \code{width}.
The number of read pairs with one read in each bin is counted for each bin
pair. The transparency/intensity of the colour of the box in the plot is
proportional to the number of counts, at least up to the specified \code{cap}.

Lowering \code{width} can improve resolution but will reduce the intensity
and spread of colours. A concomitant decrease in \code{cap} is often necessary
to ensure that the appropriate range of colours is shown. A moderate \code{cap}
is also recommended in the presence of very intense interactions, otherwise 
it will overshadow its neighbours by scaling down their intensities.

If \code{anchor} or \code{target} are \code{GRanges} objects, then
\code{anchor.lim} and \code{target.lim} do not need to be specified. Neither
needs to be consistent with the definition of the anchor and target chromosome
in \code{\link{preparePET}}; this function will automatically attempt to find
the corresponding read pairs if it suspects misspecification. 
}

\author{Aaron Lun}

\references{
Lieberman-Aiden E et al. (2009). Comprehensive mapping of long-range
interactions reveals folding principles of the human genome. \emph{Science}
326, 289-293.
}

\examples{
out1 <- system.file("exdata", "li_AA.bam", package="dacpet")
stats <- preparePET(out1, "gunkAA")

plotChIA("gunkAA", anchor=GRanges('chr1', IRanges(1e7, 12e7)), width=5e6, cap=1)
plotChIA("gunkAA", anchor=GRanges('chr1', IRanges(1e7, 12e7)), width=1e6, cap=1)
plotChIA("gunkAA", anchor=GRanges('chr1', IRanges(1e7, 12e7)), width=2e6, cap=1, diag=FALSE)

# Different ways of saying the same thing.
plotChIA("gunkAA", anchor=GRanges('chr1', IRanges(1e7, 12e7)), width=2e6, cap=5)
plotChIA("gunkAA", anchor='chr1', anchor.lim=c(1e7, 12e7), width=2e6, cap=5)

# Using a different anchor/target range.
plotChIA("gunkAA", anchor=GRanges('chr1', IRanges(1e7, 12e7)), 
	target=GRanges('chr1', IRanges(5e7, 20e7)), width=2e6, cap=5)
}
