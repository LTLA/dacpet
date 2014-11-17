\name{IList-class}
\docType{class}

\alias{IList-class}
\alias{IList}

\alias{anchors}
\alias{anchors,IList-method}
\alias{targets}
\alias{targets,IList-method}
\alias{regions}
\alias{regions,IList-method}
\alias{counts}
\alias{counts,IList-method}
\alias{info}
\alias{info,IList-method}

\alias{dim,IList-method}
\alias{dimnames,IList-method}

\alias{[,IList,ANY,ANY-method}
\alias{show,IList-method}

\alias{asDGEList}
\alias{asDGEList,IList-method}
\alias{normalize}
\alias{normalize,IList-method}

\title{IList class and methods}

\description{
Using the IList class to store pairwise genomic interactions.
}

\details{
Each IList stores counts for pairwise genomic interactions.
Slots are defined as:
\describe{
\item{\code{counts}}{an integer matrix containing the number of read pairs for an interaction in each tag file}
\item{\code{info}}{a dataframe containing assorted information about each tag file, with at least the total number of PETs in a numeric \code{totals} field}
\item{\code{anchor.id}}{an integer vector specifying the index of the first interacting region}
\item{\code{target.id}}{an integer vector specifying the index of the second interacting region}
\item{\code{region}}{a \code{GRanges} object containing the coordinates of all interacting regions}
}
Each row of \code{counts} corresponds to an interaction, while each column corresponds to a tag file.
Each interaction is parameterized by an \code{anchor.id} and \code{target.id} index, pointing to the anchor and target intervals in \code{region}, respectively.
The anchor interval is defined as that with the higher genomic coordinate compared to the target interval, to avoid redundancy from permutations.

The interpretation of each column of \code{counts} assumes that the data has been loaded from \code{\link{countPET}} or its relatives.
This changes slightly if the IList is produced by \code{\link{compressMatrix}}.
Instead of tag files, each column now represents the total count for all homo- or hetero-linkers in a single library.
}

\section{Constructor}{
\describe{
	\item{}{\code{IList(counts, info=NULL, anchors, targets, regions)}:
	Returns a IList object where each entry is used to fill the corresponding slots.
	Type coercion is performed as required.
	If \code{info} is left as \code{NULL}, the \code{totals} field is constructed from the column sums of \code{counts}.
	Values of \code{anchors}, \code{targets} and \code{regions} are used for slots \code{anchor.id}, \code{target.id} and \code{region}.
	}
}
}

\section{Accessors}{
In the code snippets below, \code{x} is a IList object.

\describe{
	\item{}{\code{anchors(x, id=FALSE)}:
	Get the \code{GRanges} corresponding to the anchor intervals for each interaction.
	If \code{id=TRUE}, indices to the corresponding intervals in \code{region} are returned instead.
	}
	
	\item{}{\code{targets(x, id=FALSE)}:
	Get the \code{GRanges} corresponding to the target intervals for each interaction.
	If \code{id=TRUE}, indices to the corresponding intervals in \code{region} are returned instead.
	}

	\item{}{\code{regions(x)}:
	Get the \code{GRanges} for all interacting regions.
	}

	\item{}{\code{counts(x)}:
	Get the matrix of counts for all interactions in all tag files.
	}

	\item{}{\code{info(x)}:
	Get the dataframe of tag file information.
	}

	\item{}{\code{dim(x)}:
	Get the dimensions, i.e., number of interactions (rows) against number of tag files (columns).
	}

	\item{}{\code{dimnames(x)}:
	Get the dimension names. This returns a list of length 2, where each element is \code{NULL} or a character vector.
	}
}
}

\section{Subsetting}{
In the code snippets below, \code{x} is a IList object.

\describe{
    \item{}{\code{x[i, j]}:
  	Get data for all interactions \code{i} in tag files \code{j}.
	Either \code{i} or \code{j} can be missing, in which case all interactions or tag files are returned, respectively.
	}
}
}

\section{Other methods}{
\describe{
	\item{}{\code{show(x)}:
		By default the \code{show} method displays count and interval data for the top 5 interactions.
		Data for all interactions will be shown if there are no more than 10 interactions.
	}
	\item{}{\code{asDGEList(x, ...)}:
		Extracts count data and totals, and returns it in the form of a \code{DGEList}.
		Arguments in \code{...} are passed to the \code{\link{DGEList}} constructor.
	}
	\item{}{\code{normalize(x, ...)}:
		Passes count data, totals and other arguments to \code{\link{normalizeCounts}} to compute normalization factors or offsets.
	}
}
}

\author{Aaron Lun}

\examples{
blah <- IList(counts=matrix(c(1,1,2,2,3,3,4,4), ncol=2), 
    info=data.frame(totals=c(10L, 10L), files=c("A", "B")), 
	anchors=c(1,2,3,4), targets=c(1,1,2,2), 
	regions=GRanges("chrA", IRanges(10+1:20, 2+21:40)))
nrow(blah)
ncol(blah)

blah
blah[1,]
blah[,2]
blah[1:2,2]

anchors(blah)
anchors(blah, id=TRUE)
targets(blah)
targets(blah, id=TRUE)
counts(blah)
info(blah)
regions(blah)

asDGEList(blah)
normalize(blah)
head(normalize(blah, type="loess"))
}

\seealso{
\code{\link{DGEList}},
\code{\link{countPET}},
\code{\link{compressMatrix}},
\code{\link{normalizeCounts}}
}
