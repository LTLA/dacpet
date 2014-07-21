\name{stripOutwardPET}
\alias{stripOutwardPET}

\title{Strip outward-facing PETs}
\description{Remove outward-facing tag pairs from a tag file}

\usage{
stripOutwardPET(file.in, file.out=file.in, min.gap=1e4, discard.to=NULL)
}

\arguments{
\item{file.in}{a character string specifying the tag file to process}
\item{file.out}{a character string specifying the path for the output tag file}
\item{min.gap}{a numeric scalar indicating the minimum gap for outward-facing read pairs}
\item{discard.to}{a character string specifying a file ppath to deposit discarded read pairs}
}

\details{
Intra-chromosomal outward-facing read pairs are often overrepresented as they
correspond to self-ligation events. These provide no information with respect
to interactions, as self-ligation events occur within a single fragment. They
must be removed to avoid any confounding effects on the subsequent analysis.

This is done by removing those pairs where the gap between reads is less
than \code{min.gap}. The choice of threshold distance can be selected by 
extracting the gap for each strand orientation with \code{\link{diagnosePET}}, 
and plotting the distribution of gaps. The chosen threshold should be 
the smallest distance at which the strand orientations are balanced.

If \code{file.out} is not specified, the operation will overwrite the tag file specified
in \code{file.in} to remove outward-facing read pairs. If \code{discard.to} is 
specified, the 5' end locations of each discarded outward-facing read pair
will be recorded in the specified file. This is useful for further analyses,
e.g., with peak calling to identify binding sites.
}

\value{
An integer scalar indicating the number of read pairs discarded is invisibly returned.
}

\seealso{
\code{\link{diagnosePET}}
}

\author{Aaron Lun}

\references{
Jin F et al. (2013). A high-resolution map of the three-dimensional chromatin interactome in human cells. \emph{Nature}, 503, 290-294
}

\examples{
out <- system.file("exdata", "li_AA.bam", package="dacpet")
stats <- preparePET(out, "gunk")

# Saves to a separate directory.
x <- stripOutwardPET("gunk", "temp")
print(x)

# Generates an ouput file.
x <- stripOutwardPET("gunk", "temp", discard.to="whee.tsv")
print(x)

# Milder removal.
x <- stripOutwardPET("gunk", "temp", min.gap=1000)
print(x) 

# Overwrites current directory.
x <- stripOutwardPET("gunk")
print(x)
}
