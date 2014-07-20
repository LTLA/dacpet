\name{preparePET}
\alias{preparePET}

\title{Prepare ChIA-PET data}
\description{Generate a tag directory from a BAM file for ChIA-PET data}

\usage{
preparePET(bam, dir, dedup=TRUE, yield=1e7, minq=NA)
}

\arguments{
\item{bam}{a character string containing a path to a name-sorted BAM file}
\item{dir}{a character string containing a path to an output directory}
\item{dedup}{a logical scalar specifying whether or not to remove duplicates}
\item{yield}{a numeric scalar indicating the number of reads to process per iteration}
\item{minq}{a numeric scalar specifying the minimum mapping quality for each read}
}

\details{
This function takes a name-sorted BAM file and generates a set of files in
\code{dir} where each file contains the 3' end positions for each tag pair.
Tags aligned to the reverse strand are denoted with a \code{-} symbol.  Each
file containts all tag pairs for a specific pair of chromosomes. This facilitates
rapid downstream processing as the BAM file does not need to be handled again. 

As previously mentioned, each file contains all tag pairs for a pair of
chromosomes.  One chromosome is referred to as the anchor whereas the other is
referred to as the target. The anchor is after the target in alphanumerical
order. The position of each tag on either chromosome is noted for each pair.
When the anchor and target refer to the same chromosome, the anchor tag is
designated as that with the larger 3' position.

Both tags in the pair must be mapped, have MAPQ above \code{minq} and not be marked as duplicates if \code{dedup=TRUE}. 
A high value for \code{minq} is generally recommended as it avoids unstable/incorrect alignments for short tag sequences. 
The default \code{NA} values means that no filtering is performed. 

Setting \code{dedup} is also recommended as the chance of random exact overlaps across the interaction space is generally low.
Any overlaps are probably technical artifacts caused by PCR duplication.

To reduce memory usage, read pairs are read in by \code{yield} at a time.
Reducing \code{yield} will reduce peak memory usage, albeit at the cost of computational speed. 
}

\value{
A tag directory is produced at \code{dir} containing files for all pairs of
chromosomes with tag pairs mapped to them. A list is also returned containing:
\item{pairs}{an integer vector containing the total number of tag pairs; the
number that is marked as duplicates; the number with scores below \code{minq};
and the number retained in the final directory}
\item{other}{an integer vector containing the number of unpaired tags, and the
number of tag groups with more than two alignments}
}

\author{Aaron Lun}

\examples{
out <- system.file("exdata", "li_AA.bam", package="dacpet")
stats <- preparePET(out, dir="gunk")
print(stats$pairs)

stats <- preparePET(out, dir="gunk", dedup=FALSE)
print(stats$pairs)

stats <- preparePET(out, dir="gunk", dedup=TRUE, minq=20)
print(stats$pairs)

stats <- preparePET(out, dir="gunk", dedup=TRUE, yield=100, minq=20)
print(stats$pairs)
}
