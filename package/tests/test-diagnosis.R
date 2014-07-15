# Tests the diagnostic components: diagnosePET, extractMito, stripOutwardPET and mergePET.

suppressPackageStartupMessages(require(dacpet))
	
all.chrs <- c(chrA=1000, chrB=2000, chrC=3000, chrM=500)
simdir <- function(dirpath, ntags) {
	unlink(list.files(dirpath, full=TRUE))
	offsets <- c(0L, cumsum(all.chrs))
	chosen.alpha <- sample(length(all.chrs), ntags, replace=TRUE)
	chosen.bravo <- sample(length(all.chrs), ntags, replace=TRUE)
	pos.alpha <- as.integer(runif(ntags, 1, all.chrs[chosen.alpha]))
	pos.bravo <- as.integer(runif(ntags, 1, all.chrs[chosen.bravo]))

	is.anchor <-chosen.alpha > chosen.bravo | (chosen.alpha==chosen.bravo & pos.alpha > pos.bravo)
	anchor.chr <- ifelse(is.anchor, chosen.alpha, chosen.bravo)
	anchor.pos <- sample(c(-1L, 1L), ntags, replace=TRUE)*ifelse(is.anchor, pos.alpha, pos.bravo)
	target.chr <- ifelse(!is.anchor, chosen.alpha, chosen.bravo)
	target.pos <- sample(c(-1L, 1L), ntags, replace=TRUE)*ifelse(!is.anchor, pos.alpha, pos.bravo)

	out <- split(data.frame(anchor.pos=anchor.pos, target.pos=target.pos), anchor.chr)
	index <- 1L
	for (x in names(out)) { 
		reout <- split(out[[x]], target.chr[as.integer(rownames(out[[x]]))])
		for (y in names(reout)) {
			fname <- paste0(index, ".gz")
			dacpet:::.saveExt(reout[[y]], file.path(dirpath, fname))
			write.table(data.frame(names(all.chrs)[as.integer(x)], names(all.chrs)[as.integer(y)], fname), file=dacpet:::.getIndex(dirpath),
				sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=index!=1L)
			index <- index + 1L
		}	
	}
	write.table(data.frame(chr=names(all.chrs), length=all.chrs), sep="\t", quote=FALSE, row.names=FALSE,
		file=dacpet:::.getLengths(dirpath))
	return(data.frame(anchor.chr=names(all.chrs)[anchor.chr], anchor.pos, target.chr=names(all.chrs)[target.chr], target.pos,
			stringsAsFactors=FALSE))	
}

###########################################################################################################
# Testing diagnosePET.

dirpath <- "stuff"
dir.create(dirpath)
diagcomp <- function(ntags, restrict=NULL) {
	truth <- simdir(dirpath, ntags)
	if (!is.null(restrict)) { truth <- truth[truth$anchor.chr %in% restrict & truth$target.chr %in% restrict,] }
	out <- diagnosePET(dirpath, restrict=restrict)

	# Checking inters.
	if (out$inter[["0"]]!=sum(truth$anchor.chr!=truth$target.chr & truth$anchor.pos > 0 & truth$target.pos > 0)) { stop("mismatch in both forward interchromosomals") }
	if (out$inter[["1"]]!=sum(truth$anchor.chr!=truth$target.chr & truth$anchor.pos < 0 & truth$target.pos > 0)) { stop("mismatch in inward interchromosomals") }
	if (out$inter[["2"]]!=sum(truth$anchor.chr!=truth$target.chr & truth$anchor.pos > 0 & truth$target.pos < 0)) { stop("mismatch in outward interchromosomals") }
	if (out$inter[["3"]]!=sum(truth$anchor.chr!=truth$target.chr & truth$anchor.pos < 0 & truth$target.pos < 0)) { stop("mismatch in both reverse interchromosomals") }
	
	# Checking gappiness.
	all.gap <- abs(truth$anchor.pos) - abs(truth$target.pos)
	if (!identical(sort(all.gap[truth$anchor.chr==truth$target.chr & truth$anchor.pos > 0 & truth$target.pos > 0]), sort(out$gap[out$flag==0L]))) { stop("mismatch in both forward intrachromosomals") }
	if (!identical(sort(all.gap[truth$anchor.chr==truth$target.chr & truth$anchor.pos < 0 & truth$target.pos > 0]), sort(out$gap[out$flag==1L]))) { stop("mismatch in both forward intrachromosomals") }
	if (!identical(sort(all.gap[truth$anchor.chr==truth$target.chr & truth$anchor.pos > 0 & truth$target.pos < 0]), sort(out$gap[out$flag==2L]))) { stop("mismatch in both forward intrachromosomals") }
	if (!identical(sort(all.gap[truth$anchor.chr==truth$target.chr & truth$anchor.pos < 0 & truth$target.pos < 0]), sort(out$gap[out$flag==3L]))) { stop("mismatch in both forward intrachromosomals") }

	return(cbind(head(out$gap), head(out$flag)))
}

set.seed(234096)
diagcomp(1000)
diagcomp(1000, restrict=c("chrA", "chrB", "chrC"))
diagcomp(1000, restrict=c("chrA", "chrB"))
diagcomp(1000, restrict=c("chrA"))
diagcomp(100)
diagcomp(100, restrict=c("chrA", "chrB", "chrC"))
diagcomp(100, restrict=c("chrA", "chrB"))
diagcomp(100, restrict=c("chrA"))
diagcomp(10)
diagcomp(10, restrict=c("chrA", "chrB", "chrC"))
diagcomp(10, restrict=c("chrA", "chrB"))
diagcomp(10, restrict=c("chrA"))
diagcomp(1)
diagcomp(1, restrict=c("chrA", "chrB", "chrC"))
diagcomp(1, restrict=c("chrA", "chrB"))
diagcomp(1, restrict=c("chrA"))

###########################################################################################################
# Testing stripOutwardPET.

stripcomp <- function(ntags, out=NULL, min.gap=100, locale=NULL) {
	truth <- simdir(dirpath, ntags)
	stripped <- stripOutwardPET(dirpath, out=out, min.gap=min.gap, discard.to=locale)
	retained <- list()
	
	# Checking output directory.
	if (is.null(out)) { out <- dirpath }
	reported <- read.table(dacpet:::.getIndex(out), stringsAsFactors=FALSE)
	for (x in 1:nrow(reported)) { 
		keep <- truth$anchor.chr==reported[x,1] & truth$target.chr==reported[x,2] 
		current <- truth[keep,]
		lost <- current$anchor.chr==current$target.chr & current$anchor.pos > 0 & current$target.pos < 0 & abs(current$anchor.pos)-abs(current$target.pos) < min.gap
		current <- current[!lost,c(2,4)]
		rownames(current) <- NULL
		blah <- read.table(file.path(out, reported[x,3]), header=TRUE)
		if (!identical(current, blah)) { stop("mismatch after stripping out tags") }
		if (any(lost)) { 
			retained[[x]] <- data.frame(reported[x,1], -truth$target.pos[keep][lost], truth$anchor.pos[keep][lost], stringsAsFactors=FALSE)
		}
	}
	all.lost <- do.call(rbind, retained)
	if (is.null(all.lost)) {
		stopifnot(stripped==0L)
	} else {
		stopifnot(identical(stripped, nrow(all.lost)))
	}

	# Checking the discard location.
	if (stripped && !is.null(locale)) { 
		o <- do.call(order, all.lost)
		all.lost <- all.lost[o,]
		reported <- read.table(locale, colClasses=list("character", "integer", "integer"), stringsAsFactors=FALSE)
		o <- do.call(order, reported)
		reported <- reported[o,]
		colnames(all.lost) <- colnames(reported) <- c("chr", "start", "end")
		rownames(all.lost) <- rownames(reported) <- NULL
		if (!identical(all.lost, reported)) { stop("mismatch in the discarded tags") }
	} else if (stripped==0L) { 
		out <- tryCatch({
		   read.table(locale)	
		}, error=function(cond) { return(FALSE) })
		if (!is.logical(out)) { stop("non-empty discard list") }
	}
	return(stripped)
}

set.seed(8573910)
tmpdir <- tempfile(tmpdir=".")
discardee <- "whee.tsv" #tempfile(tmpdir=".", fileext=".tsv")

stripcomp(1000, out=NULL, min.gap=100, locale=NULL)
stripcomp(1000, out=tmpdir, min.gap=100, locale=NULL)
unlink(tmpdir, recursive=TRUE)
stripcomp(1000, out=NULL, min.gap=100, locale=discardee)

stripcomp(1000, out=NULL, min.gap=500, locale=NULL)
stripcomp(1000, out=tmpdir, min.gap=500, locale=NULL)
unlink(tmpdir, recursive=TRUE)
stripcomp(1000, out=NULL, min.gap=500, locale=discardee)

stripcomp(1000, out=NULL, min.gap=1000, locale=NULL)
stripcomp(1000, out=tmpdir, min.gap=1000, locale=NULL)
unlink(tmpdir, recursive=TRUE)
stripcomp(1000, out=NULL, min.gap=1000, locale=discardee)

stripcomp(100, out=NULL, min.gap=100, locale=NULL)
stripcomp(100, out=tmpdir, min.gap=100, locale=NULL)
unlink(tmpdir, recursive=TRUE)
stripcomp(100, out=NULL, min.gap=100, locale=discardee)

stripcomp(10, out=NULL, min.gap=1000, locale=NULL)
stripcomp(10, out=tmpdir, min.gap=1000, locale=NULL)
unlink(tmpdir, recursive=TRUE)
stripcomp(10, out=NULL, min.gap=1000, locale=discardee)

###########################################################################################################
# Testing extractMito.

mitocomp <- function(ntags) {
	truth <- simdir(dirpath, ntags)
	collected <- c(truth$anchor.chr[truth$target.chr=="chrM"],
					truth$target.chr[truth$anchor.chr=="chrM" & truth$target.chr!="chrM"])

	mito <- extractMito(dirpath, self=TRUE)
	for (x in names(mito)) { if (!identical(mito[[x]], sum(collected==x))) { stop("mismatch in mitocounts") } }

	mito <- extractMito(dirpath, self=TRUE, restrict=c("chrA", "chrB"))
	for (x in names(mito)) { if (!identical(mito[[x]], sum(collected==x))) { stop("mismatch in mitocounts") } }

	mito <- extractMito(dirpath)
	for (x in names(mito)) { if (!identical(mito[[x]], sum(collected==x & x!="chrM"))) { stop("mismatch in mitocounts") } }

	return(mito)
}

set.seed(475632)
mitocomp(1000)
mitocomp(1000)
mitocomp(100)
mitocomp(100)
mitocomp(10)
mitocomp(10)

###########################################################################################################
# Testing mergePET.

altpath <- "second.stuff"
finpath <- "final.stuff"
dir.create(altpath)
mergecomp <- function(ntags1, ntags2) {
	t1 <- simdir(dirpath, ntags1)
	t2 <- simdir(altpath, ntags2)
	mergePET(c(dirpath, altpath), out=finpath)

	# Checking by ensuring that everyone is present who is meant to be.
	ref <- read.table(dacpet:::.getLengths(dirpath), header=TRUE)
	sec <- read.table(dacpet:::.getLengths(altpath), header=TRUE)
	final <- read.table(dacpet:::.getLengths(finpath), header=TRUE)
	if (!identical(ref, sec) || !identical(ref, final)) { stop("mismatches in the length files") }
	
	# Running through the merged directory, and ensuring everyone is accounted for.
	indices <- read.table(dacpet:::.getIndex(finpath))
	hits1 <- integer(nrow(t1))
	hits2 <- integer(nrow(t2))
	for (x in 1:nrow(indices)) { 
		anchor <- indices[x,1]
		target <- indices[x,2]
		current <- read.table(file.path(finpath, indices[x,3]), header=TRUE)
		current <- current[do.call(order, current),]

		keep <- t1$anchor.chr==anchor & t1$target.chr==target
		current1 <- t1[keep,c(2,4)]
		hits1[keep] <- hits1[keep] + 1L
		keep <- t2$anchor.chr==anchor & t2$target.chr==target
		current2 <- t2[keep,c(2,4)]
		hits2[keep] <- hits2[keep] + 1L

		finality <- rbind(current1, current2)
		finality <- finality[do.call(order, finality),]
		rownames(finality) <- rownames(current) <- NULL
		if (!identical(current, finality)) { stop("merged coordinates are inconsistent with originals") }
	}

	if (!all(hits1==1L) || !all(hits2==1L)) { stop("not all original coordinates used") }
	return(head(t1))
}

set.seed(1002)
mergecomp(1000, 1000)
mergecomp(1000, 100)
mergecomp(1000, 10)
mergecomp(1000, 1)
mergecomp(100, 100)
mergecomp(10, 10)
mergecomp(1, 10)

###########################################################################################################
# Cleaning up the mess.

unlink(dirpath, recursive=TRUE)
unlink(altpath, recursive=TRUE)
unlink(finpath, recursive=TRUE)
unlink(discardee)

###########################################################################################################
