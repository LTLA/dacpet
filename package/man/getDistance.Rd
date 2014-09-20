\name{getDistance}
\alias{getDistance}
\alias{getAnchor}
\alias{getTarget}

\title{Get distance between pairs}
\description{Compute the gap distance between the regions involved in an interaction}

\usage{
getDistance(data)
}

\arguments{
\item{data}{an \code{IList} of interaction data like that produced by \code{\link{countPET}} or \code{\link{recountPET}}}
}

\value{
An integer vector specifying the gap distance for each pair of regions in \code{data}.
}

\details{
The \code{getDistance} function computes the gap distance between the interacting regions for each pair of regions. 
The gap refers to the number of base pairs between the end of one region and the start of the next. 
Negative values indicate the span of the overlap whereas \code{NA} values mark inter-chromosomal pairs of regions.
}

\seealso{
	\code{\link{countPET}},
	\code{\link{recountPET}}
}

\author{Aaron Lun}

\examples{
out1 <- system.file("exdata", "li_AA.bam", package="dacpet")
stats <- preparePET(out1, "gunkAA")
data <-countPET("gunkAA", width=1000000, filter=0L)
summary(getDistance(data))
}

