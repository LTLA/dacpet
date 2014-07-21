mergePET <- function(files, file.out)
# Merges a bunch of count directories together, usually to deal with technical
# replicates.
#
# written by Aaron Lun
# 15 June, 2014
{
	overall <- .loadIndices(files)
	file.tmp <- tempfile(tmpdir=".", fileext=".h5")
	on.exit({ if (file.exists(file.tmp)) { unlink(file.tmp) } })
	sofar<-1L
	h5createFile(file.tmp)
	h5createGroup(file.tmp, "counts")

	# Merging for each combination, as necessary.
	for (ac in names(overall)) {
		current<-overall[[ac]]
		h5createGroup(file.tmp, file.path('counts', ac))
		for (tc in names(current)) {
			is.okay <- current[[tc]]

			all.tab <- list()
			for (x in 1:length(is.okay)) { 
				current.tab <- NULL
				if (is.okay[x]) { current.tab <- .getPairs(files[x], ac, tc) }
				all.tab[[x]] <- current.tab
			}
			all.tab <- do.call(rbind, all.tab)
			h5write(all.tab, file.tmp, file.path("counts", ac, tc))						
		}
	}
	
	# Checking the lengths
	my.lengths <- NULL
	for (x in files) {
		cur.lengths <- h5read(x, "lengths")
		cur.lengths <- cur.lengths[order(cur.lengths$chr),]
		if (is.null(my.lengths)) { 
			my.lengths <- cur.lengths
		} else {
			stopifnot(identical(my.lengths, cur.lengths))						
		}
	}
	h5write(my.lengths, file.tmp, "lengths")

	# Moving to the output.
	file.rename(file.tmp, file.out)
	return(invisible(NULL))
}
