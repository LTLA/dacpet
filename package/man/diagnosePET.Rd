\name{diagnosePET}
\alias{diagnosePET}

\title{Diagnose problems with ChIA-PET data}
\description{Extract strand and gap diagnostic information from a processed directory}

\usage{
diagnosePET(dir, restrict=NULL)
}

\arguments{
\item{dir}{a character string specifying a file path to a tag directory}
\item{restrict}{a character vector containing the chromosomes to be analyzed}
}

\details{
Each intra-chromosomal read/tag pair can be parameterised in terms of the gap
distance between the 333' ends of each tag. Smaller gaps correspond to local
interactions whereas larger gaps correspond to distal interactions. Of course,
inter-chromosomal read pairs are not reported as the gap has no meaning.

Each intra-chromosomal read pair can also be reported in terms of the relative
orientation of the strands. A code of 1 indicates inward-facing read pairs, 
whereas a code of 2 indicates outward-facing read pairs. Codes of 0 or 3 
indicate read pairs that are both on the forward and reverse strand, respectively.

The distribution of gap sizes can be plotted for each strand orientation, to 
identify the minimum gap distance at which outward-facing read pairs should be
retained. Any closer read pairs correspond to self-ligation events and should
be removed with \code{\link{stripOutwardPET}}. The minimum size threshold is
that where the distributions for all strand orientations are the same.

The total number of inter-chromosomal read pairs is also reported. For the
purposes of diagnostics, one inter-chromosomal read pair is indistinguishable
from the next given that the relative definition of strandedness is arbitrary on 
each chromosome.
}

\value{
A list containing:
\item{flag}{a vector of integer codes, specifying the relative strand orientation for intra-chromosomal read pairs}
\item{gap}{an integer vector containing the gap distances for intra-chromosomal read pairs of each orientation}
}

\seealso{
\code{\link{stripOutwardPET}}
}

\author{Aaron Lun}

\references{
Jin F et al. (2013). A high-resolution map of the three-dimensional chromatin interactome in human cells. \emph{Nature}, 503, 290-294
}

\examples{
out <- system.file("exdata", "li_AA.bam", package="dacpet")
stats <- preparePET(out, dir="gunk")
blah <- diagnosePET("gunk")
print(blah$inter)

lgap <- log1p(blah$gap)
all.x <- seq(0, max(lgap)+1, length.out=20)
inwards <- hist(lgap[blah$flag==1L], breaks=all.x, plot=FALSE)
outwards <- hist(lgap[blah$flag==2L], breaks=all.x, plot=FALSE)
samestr <- hist(lgap[blah$flag==0L | blah$flag==3L], breaks=all.x, plot=FALSE)

plot(inwards$mids, outwards$counts, type="l", col="red", xlab="Log Distance", ylab="Frequency")
lines(outwards$mids, inwards$counts, col="blue")
lines(samestr$mids, samestr$counts/2, col="darkgreen")
legend("topright", c("Outward", "Inward", "Same str"), col=c("red", "blue", "darkgreen"), lwd=2)
}
