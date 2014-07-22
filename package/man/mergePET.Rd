\name{mergePET}
\alias{mergePET}

\title{Merge tag files}
\description{Merge tag files generated from multiple libraries}

\usage{
mergePET(files, file.out)
}

\arguments{
\item{files}{a character vector containing paths to tag files to be merged}
\item{file.out}{a character string specifying the path to store the merged file}
}

\details{
This simply concatenates the objects in each HDF5 tag file corresponding to the
same pair of chromosomes, and produces a new file in \code{file.out}. It does so
for each pair of chromosomes. This is useful for merging together technical
replicates prior to downstream processing. Note that the user is responsible for 
ensuring that all merged files have the same linker combination.
}

\value{An invisible \code{NULL} is returned}

\seealso{
\code{\link{preparePET}}
}

\author{Aaron Lun}

\examples{
out1 <- system.file("exdata", "li_AA.bam", package="dacpet")
stats <- preparePET(out1, "gunkAA")
out2 <- system.file("exdata", "li_BB.bam", package="dacpet")
stats <- preparePET(out2, "gunkBB")
mergePET(c("gunkAA", "gunkAB"), "gunk_homo")
}
