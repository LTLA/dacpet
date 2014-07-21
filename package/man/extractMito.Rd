\name{extractMito}
\alias{extractMito}

\title{Extract mitochondrial counts}
\description{Extract tag pair counts for interactions between the mitochondrial and nuclear chromosomes}

\usage{
extractMito(file, self=FALSE, restrict=NULL, name="chrM")  
}

\arguments{
\item{file}{a character string specifying a file path to a tag file}
\item{self}{a logical scalar indicating whether read pairs mapped wholly to the mitochondrial genome should be counted}
\item{restrict}{a character vector containing the nuclear chromosomes to be analyzed}
\item{name}{character string specifying the name of the mitochondrial genome}
}

\details{
In eukaryotic organisms, the mitochondria and the nucleus are two separate
organelles.  This means that there should be no ``interactions'' between the
DNA in each organelle. Any read pairs mapped between organelles should
correspond to non-specific ligation events. Counts can be extracted to see
whether this is the case.

Specifically, counts can be extracted for each interaction between the
mitochondrial genome as a whole and any nuclear chromosome in \code{restrict}
(all are used by default). If \code{self=TRUE}, intra-chromosomal read pairs on
the mitochondrial genome will also be counted. This should be performed for the
tag file corresponding to each linker combination for a single library.

Analysis of the counts usually proceeds with Pearson's chi-squared test,
to check whether or not a random ligation model is appropriate. More general
analyses will examine the fold change between the homo- and hetero-linker
counts and see how well they correlate to the rest of the dataset.
}

\seealso{
\code{\link{chisq.test}}
}

\author{
Aaron Lun
}

\examples{
out <- system.file("exdata", "li_AA.bam", package="dacpet")
stats <- preparePET(out, "gunk")
sum(extractMito("gunk"))
sum(extractMito("gunk", self=TRUE))

# Can also use it for other chromosomes.
sum(extractMito("gunk", name="chrY"))
sum(extractMito("gunk", self=TRUE, name="chr1"))

# Pretend to have a sum of mito:nuclear counts for AA, AB and BB files.
all.aa <- 22
all.ab <- 57
all.bb <- 53
all.counts <- c(all.aa, all.ab, all.bb)
prop.a <- (all.aa*2+all.ab)/sum(all.counts)/2
prop.b <- 1-prop.a
chisq.test(all.counts, p=c(prop.a^2, 2*prop.a*prop.b, prop.b^2))
}
