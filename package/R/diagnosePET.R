diagnosePET <- function(file, restrict=NULL) 
# This function pulls up PETs and returns their orientation and insert size.
# The idea is to allow some sort of diagnostics to be performed on each
# directory based on strand balance and gap data.
#
# written by Aaron Lun
# 22 November 2013
{
    overall <- .loadIndices(file)
	all.ori <- list()
	all.gap <- list()
	all.inter <- integer(4)

	# Running through each pair of chromosomes.
	for (anchor in names(overall)) {
		if (!is.null(restrict) && !(anchor %in% restrict)) { next }
		current<-overall[[anchor]]
        for (target in names(current)) {
	
			if (!is.null(restrict) && !(target %in% restrict)) { next }
			if (!current[[target]]) { next }
			stuff <- .getPairs(file, anchor, target)
			flagged <- .getFlag(stuff$anchor.pos, stuff$target.pos)
	
			if (anchor==target) { 
				all.gap[[length(all.gap)+1L]] <- .getGap(stuff$anchor.pos, stuff$target.pos)
				all.ori[[length(all.ori)+1L]] <- flagged 
			} else {
				all.inter <- all.inter + tabulate(flagged + 1L, nbins=length(all.inter))
			}
		}
	}

	# Storing the output.
	all.gap<-unlist(all.gap)
	all.ori<-unlist(all.ori)
	if (!length(all.gap)) {
		all.ori<-integer(0)
		all.gap<-integer(0)
	}
	names(all.inter) <- 0:3
	return(list(flag=all.ori, gap=all.gap, inter=all.inter))
}

.getFlag <- function(ax, tx) { ifelse(ax > 0L, 0L, 1L)+ifelse(tx > 0L, 0L, 2L) }
.getGap <- function(ax, tx) { ifelse(ax > 0L, ax, -1L*ax) - ifelse(tx > 0L, tx, -1L*tx) }

