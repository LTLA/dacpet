mergePET <- function(files, file.out)
# Merges a bunch of count directories together, usually to deal with technical
# replicates.
#
# written by Aaron Lun
# 15 June, 2014
{
	if (length(files)==0L) { stop("must specify non-empty vector of input files") }
	if (length(files)==1L) { 
		file.copy(files, file.out)
		return(invisible(NULL))
	}
	overall <- .loadIndices(files)
	file.tmp <- tempfile(tmpdir=".", fileext=".h5")
	on.exit({ if (file.exists(file.tmp)) { unlink(file.tmp) } })
	
	# Checking the lengths
	my.lengths <- .getLengths(files[1])
	my.lengths <- my.lengths[order(my.lengths$chr),]
	for (x in files[-1]) {
		cur.lengths <- .getLengths(x)
		cur.lengths <- cur.lengths[order(cur.lengths$chr),]
		stopifnot(identical(my.lengths, cur.lengths))
	}

	# Merging for each combination, as necessary.
	.initializeH5(file.tmp, my.lengths)
	for (ac in names(overall)) {
		current<-overall[[ac]]
		.addGroup(file.tmp, ac)
		for (tc in names(current)) {
			is.okay <- current[[tc]]

			all.tab <- list()
			for (x in 1:length(is.okay)) { 
				current.tab <- NULL
				if (is.okay[x]) { current.tab <- .getPairs(files[x], ac, tc) }
				all.tab[[x]] <- current.tab
			}
			all.tab <- do.call(rbind, all.tab)
			.writePairs(file.tmp, ac, tc, all.tab)
		}
	}
	
	# Moving to the output.
	if (!file.rename(file.tmp, file.out)) { stop("failed to copy temporary file to the specified destination") }
	return(invisible(NULL))
}
