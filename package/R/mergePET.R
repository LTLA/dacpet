mergePET <- function(dirs, out)
# Merges a bunch of count directories together, usually to deal with technical
# replicates.
#
# written by Aaron Lun
# 15 June, 2014
{
	overall <- .loadIndices(unlist(dirs))
	tmpdir<-tempfile(tmpdir=".")
	dir.create(tmpdir)
	on.exit({ if (file.exists(tmpdir)) { unlink(tmpdir, recursive=TRUE) } })
	newindex<-.getIndex(tmpdir)
	sofar<-1L
	
	# Merging for each combination, as necessary.
	for (ac in names(overall)) {
		current<-overall[[ac]]
		for (tc in names(current)) {
			fnames<-current[[tc]]

			all.tab <- list()
			for (x in 1:length(fnames)) { 
				current.tab <- NULL
				if (nchar(fnames[x])) { 
					current.tab <-read.table(file.path(dirs[x], fnames[x]), header=TRUE, colClasses="integer", comment.char="")
				}
				all.tab[[x]] <- current.tab
			}
			all.tab <- do.call(rbind, all.tab)

			ofname<-paste0(sofar, ".gz")
			.saveExt(all.tab, file.path(tmpdir, ofname))
			write.table(file=newindex, data.frame(ac, tc, ofname),
				row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE, quote=FALSE)
			sofar<-sofar+1L
		}
	}
	
	# Checking the lengths
	my.lengths <- NULL
	for (x in dirs) {
		cur.lengths <- read.table(.getLengths(x), header=TRUE, stringsAsFactors=FALSE)
		cur.lengths <- cur.lengths[order(cur.lengths$chr),]
		if (is.null(my.lengths)) { 
			my.lengths <- cur.lengths
		} else {
			stopifnot(identical(my.lengths, cur.lengths))						
		}
	}
	file.copy(.getLengths(dirs[1]), .getLengths(tmpdir))

	# Moving to the output.
	if (file.exists(out)) { unlink(out, recursive=TRUE) }
	file.rename(tmpdir, out)
	return(invisible(NULL))
}
