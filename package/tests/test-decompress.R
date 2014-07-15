# This script tests the decompression capabilities of 'decompressMatrix'.

suppressPackageStartupMessages(require(dacpet))

comp <- function(nlibs=3, ndirs=10) {
	# Generates something like data out of countPET.
	data <- list()
	cur.lib <- as.character(sample(nlibs, ndirs, replace=TRUE))
	cur.het <- rbinom(ndirs, 1, 0.5)==1L 
	npairs <- 1000
	data$counts <- matrix(as.integer(runif(npairs*ndirs,0, 10)), ncol=ndirs, nrow=npairs)
	data$totals <- as.integer(runif(ndirs, 1000, 1500))
	
	my.chrA <- sample(2, npairs, replace=TRUE)
	my.chrT <- sample(2, npairs, replace=TRUE)
	clen <- 50L
	as <- as.integer(runif(npairs, 1L, clen+0.99)) + ifelse(my.chrA==1L, 0L, clen)
	ts <- as.integer(runif(npairs, 1L, clen+0.99)) + ifelse(my.chrT==1L, 0L, clen)
	data$pairs <- data.frame(anchor=as, target=ts)
	range.start <- c(1:clen, 1:clen)
	data$region <- GRanges(rep(c("chrA", "chrB"), c(clen, clen)), IRanges(range.start, range.start+10L))

	# Basic decompression. 
   	proposed <- compressMatrix(data, libname=cur.lib, hetero=cur.het)
	temp <- proposed$counts
	totes <- proposed$totals
	for (d in 1:ndirs) {
		chosen <- cur.het[d]==proposed$hetero & cur.lib[d]==proposed$libname
		if (sum(chosen)!=1L) { stop("multiple matches not valid") }
		temp[,chosen] <- temp[,chosen] - data$counts[,d]
		chosen <- cur.lib[d]==proposed$libname
		totes[chosen] <- totes[chosen] - data$totals[d]
	}
	if (any(temp!=0L)) { stop("mismatch in accumulated counts") }
	if (any(totes!=0L)) { stop("mismatch in accumulated totals") }

	return(head(proposed$counts))
}

####################################################################################
# Actually running it.

set.seed(489176)

comp(3, 5)
comp(3, 10)
comp(3, 10)
comp(3, 10)
comp(4, 8)
comp(4, 5)
comp(2, 2)
comp(1, 2)
comp(1, 2)
comp(1, 1)

comp(5, 5)
comp(5, 10)
comp(5, 10)
comp(5, 10)
comp(5, 8)
comp(5, 5)
comp(5, 2)
comp(5, 2)
comp(5, 1)

####################################################################################
# End.

