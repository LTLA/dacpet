compressMatrix <- function(data, hetero=NULL, libname=NULL) 
# This converts the raw output from countPET into something that's 
# actually usable in edgeR. Briefly, it gets rid of weighting by 
# expanding out all relevant elements. It also gets rid of the majority
# of uninteresting and highly local interactions, based on the
# 'min.dist' parameter. The latter ensures that the expansion won't
# use an outrageous amount of memory. All in all, this is easier to
# manipulate as weights no longer need consideration.
#
# written by Aaron Lun
# 21 January, 2014
{
    if (is.null(hetero)) { hetero <- grepl("AB", data$files) }
	if (is.null(libname)) { 
		libname <- sub("_(AA|AB|BB|other)(\\.[^\\.]+)?$", "", data$files) 
		if (any(libname==data$files)) { warning("possible error in automatic library name identification") }
	}

	ndir <- length(hetero)
	if (ndir!=length(libname) || ndir!=ncol(data$counts)) { stop("inconsistent number of directories in libname, hetero and data") }
	o <- order(libname, hetero)
    hetero <- hetero[o]
	libname <- libname[o]

	# Collapsing columns of the count matrix according to the homo/hetero specification.
	new.counts <- list(data$counts[,o[1]])
	lib.out <- list(libname[1])
	status <- list(hetero[1])
	if (ndir>=2L) {
		index <- 1L
	    for (x in 2:ndir) {
			cur.counts <- data$counts[,o[x]]
		    if (hetero[x]==hetero[x-1L] && libname[x]==libname[x-1L]) {
			    new.counts[[index]] <- new.counts[[index]] + cur.counts
			} else {
				index <- index + 1L
				new.counts[[index]] <- cur.counts
				lib.out[[index]] <- libname[x]
				status[[index]] <- hetero[x]
			}
		}
	}
		
	# Computing the total sum for each library.
	sum.totals <- lapply(split(data$totals[o], libname), FUN=sum)
	new.totals <- integer(length(lib.out))
	for (x in 1:length(lib.out)) { new.totals[x] <- sum.totals[[lib.out[[x]]]] }

	# Cleaning up the output.
 	new.counts <- do.call(cbind, new.counts)
	colnames(new.counts) <- paste0(lib.out, ifelse(status, "het", "hom"))
	rownames(new.counts) <- rownames(data$counts)
	
	output <- list()
	output$hetero <- unlist(status)
	output$libname <- unlist(lib.out)
	output$totals <- new.totals
	
	# Expanding rows of the count matrix, if necessary.
	output$counts <- new.counts		
	output$pairs <- data$pairs
   	output$region <- data$region
	return(output)
}
