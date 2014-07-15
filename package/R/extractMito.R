extractMito <- function(dir, self=FALSE, restrict=NULL, name="chrM")
# Pretty much does what it says. Extracts all counts for interactions
# involving the mitochondrial genome. By default, it doesn't extract
# counts involving interactions between the mitochondrial genome and 
# itself. It can also be regulated further by specifying a non-trivial
# restrict argument.
# 
# written by Aaron Lun
# 21 January, 2014
{
	indices <- read.table(.getIndex(dir), header=FALSE, stringsAsFactors=FALSE, comment.char="")
	if (is.null(restrict)) { 
		restrict <- read.table(.getLengths(dir), header=TRUE, stringsAsFactors=FALSE)$chr 
		mito <- indices[,1]==name | indices[,2]==name 
	} else {
		mito <- (indices[,1]==name & indices[,2] %in% restrict) | (indices[,2]==name & indices[,1] %in% restrict)
	}
	if (!self) { mito <- mito & indices[,1]!=indices[,2] }
	fnames <- indices[mito,3]
	relevants <- ifelse(indices[mito,1]==name, indices[mito,2], indices[mito,1])
	if (anyDuplicated(relevants)) { stop("multiple entries in the index file for a combination") }
		
	# Pulling out everything for each file.
	outed <- list()
	is.present <- match(restrict, relevants)
	for (x in 1:length(restrict)) {
		y <- is.present[x]
		if (is.na(y)) { 
			outed[[x]] <- 0L
		} else {
			outed[[x]] <- nrow(read.table(file.path(dir, fnames[y]), header=TRUE, colClasses="integer", stringsAsFactors=FALSE))
		}
	}
	
	# Deciciding what to return.
	outed <- unlist(outed)
	names(outed) <- restrict
	return(outed)
}
