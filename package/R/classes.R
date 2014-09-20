# Defines the IList class, that will be used to hold various information.

setClass("IList", representation(counts="matrix", info="data.frame", anchor.id="integer", target.id="integer", region="GRanges"))

setValidity("IList", function(object) {
	if (nrow(object@counts)!=length(object@anchor.id)) {
		return('rows in count matrix not equal to length of anchor vector')
	} 
	if (nrow(object@counts)!=length(object@target.id)) { 
		return('rows in count matrix not equal to length of target vector')
	}
	if (ncol(object@counts)!=nrow(object@info)) { 
		return('columns of count matrix not equal to rows of information data.frame')
	}
	if (! "totals" %in% colnames(object@info) || !is.numeric(object@info$totals)) { 
		return("information data.frame should contain a numeric 'totals' field")
	}

	if (!all(object@anchor.id >= 1L)) { 
		return('not all anchors are positive integers')
	} 
	if (!all(object@target.id >= 1L)) {
		return('not all targets are positive integers')
	}
	if (!all(object@anchor.id <= length(object@region))) {
		return('not all anchors refer to valid regions')
	} 
	if (!all(object@target.id <= length(object@region))) { 
		return('not all targets refer to valid regions')
	}
	if (!all(object@anchor.id >= object@target.id)) { 
		return('target indices cannot be greater than anchor indices')
	}
	return(TRUE)
})

setMethod("initialize", signature("IList"), function(.Object, ...) {
	value <- callNextMethod()
	validObject(value)
	value
})

setMethod("show", signature("IList"), function(object) {
	total <- nrow(object@counts)
	leftover <- 0
	if (total > 10) {
		toshow <- 5
		leftover <- total - toshow
	} else {
		toshow <- total
	}
	nregs <- length(object@region)
	nlibs <- ncol(object@counts)
	cat("IList object for", nlibs, ifelse(nlibs==1L, "library", "libraries"), 
		"with", total, ifelse(total==1L, "pair", "pairs"), "across", 
		nregs, ifelse(nregs==1L, "region\n", "regions\n"))
	cat("\n")
	
	cat("Counts:\n")
	print(head(object@counts, toshow))
	if (leftover) { cat("... and", leftover, "more rows\n") }
	cat("\n")
	
	cat("Information:\n")
	print(object@info)
	cat("\n")

	a.show <- object@region[head(object@anchor.id, toshow)]
	cat("Anchors:\n")
	print(data.frame(Chr=seqnames(a.show), Start=start(a.show), End=end(a.show)))
	if (leftover) { cat("... and", leftover, "more rows\n") }
	cat("\n")

	cat("Targets:\n")
	t.show <- object@region[head(object@target.id, toshow)]
	print(data.frame(Chr=seqnames(t.show), Start=start(t.show), End=end(t.show)))
	if (leftover) { cat("... and", leftover, "more rows\n") }
})

# Assorted subsetting methods.
setMethod("[", "IList", function(x, i, j, ..., drop=TRUE) {
	if (missing(i)) {
		new.counts <- x@counts				
		new.anchors <- x@anchor.id
		new.targets <- x@target.id	
	} else {
		new.counts <- x@counts[i,,drop=FALSE]	
		new.anchors <- x@anchor.id[i]
		new.targets <- x@target.id[i]
	}

	if (missing(j)) { 
		new.info <- x@info
	} else {
		new.counts <- new.counts[,j,drop=FALSE]
		new.info <- x@info[j,]
	}
	initialize(x, counts=new.counts, info=new.info,
		anchor.id=new.anchors, target.id=new.targets,
		region=x@region)
})

# Some getters. No need for setters, really.
setGeneric("anchors", function(object, ...) { standardGeneric("anchors") })
setMethod("anchors", signature("IList"), function(object, id=FALSE) {
	if (id) { return(object@anchor.id) }
	object@region[object@anchor.id]
})

setGeneric("targets", function(object, ...) { standardGeneric("targets") })
setMethod("targets", signature("IList"), function(object, id=FALSE) {
	if (id) { return(object@target.id) }
	object@region[object@target.id]
})

setGeneric("counts", function(object) { standardGeneric("counts") })
setMethod("counts", signature("IList"), function(object) {
	object@counts
})

setGeneric("info", function(object) { standardGeneric("info") })
setMethod("info", signature("IList"), function(object) { 
	object@info
})

setGeneric("regions", function(object) { standardGeneric("regions") })
setMethod("regions", signature("IList"), function(object) {
	object@region
})

setMethod("dim", signature("IList"), function(x) {
	dim(x@counts)
})

setMethod("dimnames", signature("IList"), function(x) {
	dimnames(x@counts)
})

# Constructor object.
IList <- function(counts, info=NULL, anchors, targets, regions) {
	if (!is.integer(counts)) { storage.mode(counts) <- "integer" }
	anchors <- as.integer(anchors)
	targets <- as.integer(targets)
	if (is.null(info)) { info <- data.frame(totals=colSums(counts)) }
	new("IList", counts=counts, info=info, anchor.id=anchors, target.id=targets, region=regions)
}

setGeneric("asDGEList", function(object, ...) { standardGeneric("asDGEList") })
setMethod("asDGEList", signature("IList"), function(object, ...) {
	DGEList(counts(object), lib.size=info(object)$totals, ...)
})

