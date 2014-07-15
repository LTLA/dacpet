\name{mergePET}
\alias{mergePET}

\title{Merge tag directories}
\description{Merge tag directories from multiple libraries}

\usage{
mergePET(dirs, out)
}

\arguments{
\item{dirs}{a character vector containing paths to tag directories to be merged}
\item{out}{a character string specifying the path to store the merged directory}
}

\details{
This simply concatenates the files in each tag directory corresponding to the
same pair of chromosomes, and produces a new file in \code{out}. It does so
for each pair of chromosomes. This is useful for merging together technical
replicates prior to downstream processing. Note that the user is responsible for 
ensuring that all merged directories have the same linker combination.
}

\value{An invisible \code{NULL} is returned}

\seealso{
\code{\link{preparePET}}
}

\author{Aaron Lun}

\examples{
out1 <- system.file("exdata", "li_AA.bam", package="dacpet")
stats <- preparePET(out1, dir="gunkAA")
out2 <- system.file("exdata", "li_BB.bam", package="dacpet")
stats <- preparePET(out2, dir="gunkBB")
mergePET(c("gunkAA", "gunkAB"), "gunk_homo")
}
