countstr <- 'counts'
lengthstr <- 'lengths'

.loadIndices <- function(y)
# A quick and dirty function for index loading with multiple libraries. This
# produces a list which describes the necessary HDF5 object corresponding to
# each chromosome combination for each library. 
{
	overall <- list()
	ni<-length(y)
	for (ix in 1:ni) {
		current <- h5ls(y[ix])
		keep <- grepl(paste0("^/", countstr), dirname(current$group)) & current$otype=="H5I_DATASET"
		all.anchors <- basename(current$group[keep])
		assorted.info <- current[keep,c("name", "dim")]

		current <- split(assorted.info, all.anchors)
		for (ac in names(current)) {
			if (is.null(overall[[ac]])) { overall[[ac]]<-list() }
			subcurrent <- current[[ac]]
			subcurrent <- split(subcurrent$dim, subcurrent$name)
			for (tc in names(subcurrent)) {
				if (is.null(overall[[ac]][[tc]])) { overall[[ac]][[tc]] <- integer(ni) }
				overall[[ac]][[tc]][ix] <- as.integer(subcurrent[[tc]]) # Dims will be convertible as it stores the number of rows in a data.frame.
			}
		}
	}
	return(overall)
}

.getPairs <- function(y, anchor, target) { h5read(y, file.path(countstr, anchor, target)) }

.getLengths <- function(y) { h5read(y, lengthstr) }

.initializeH5 <- function(y, len) {
	if (file.exists(y)) { unlink(y, recursive=TRUE) } 
	if (!h5createFile(y)) { stop(sprintf("failed to create '%s'", y)) }
	if (!h5createGroup(y, countstr)) { stop(sprintf("failed to add 'counts' group to '%s'", y)) }
	if (h5write(len, y, lengthstr)) { stop(sprintf("failed to add length data to '%s'", y)) }
	return(invisible(NULL))
}

.addGroup <- function(y, anchor) {
	if (!h5createGroup(y, file.path(countstr, anchor))) { stop("failed to add '%s' group to '%s'", anchor, y) }
	return(invisible(NULL))
}

.writePairs <- function(y, anchor, target, pairs) {
	if (h5write(pairs, y, file.path(countstr, anchor, target))) { stop("failed to add tag pair data to '%s'", y) }
	return(invisible(NULL))
}
