marginPET <- function(files, regions, ext=1L, restrict=NULL)
# This counts the number of PETs with each end overlapping one of 
# the specified regions. It then reports counts for each pair of
# regions with counts above the specified threshold.
#
# written by Aaron Lun
# 23 January, 2014
{
    nlibs<-length(files)
    if (!is.integer(ext)) { ext<-as.integer(ext) }
    if (nlibs==0) {
        stop("number of libraries must be positive")
    } else if (ext < 0) {
        stop("fragment length must be a non-negative integer")
    }
	gr <- split(ranges(regions), seqnames(regions))
	oridex <- split(1:length(regions), seqnames(regions))

	chromosomes <- .getChrs(files)
	chrs <- names(chromosomes)
	my.chrs <- names(oridex)

    # Running through each pair of chromosomes.
    overall <- .loadIndices(files)
    all.anchors <- all.targets <- list()
	all.counts <- list()
	ix <- 1L
	totals <- integer(nlibs)

	all.margins <- matrix(0L, nrow=length(regions), ncol=nlibs)
    for (anchor in names(overall)) {
		stopifnot(anchor %in% chrs)
		if (!is.null(restrict) && !(anchor %in% restrict)) { next }
        current<-overall[[anchor]]
        for (target in names(current)) {
			stopifnot(target %in% chrs)
			if (!is.null(restrict) && !(target %in% restrict)) { next }
            is.okay <- current[[target]]

			# Deciding whether or not to skip after computing the total counts.
			if (! (anchor %in% my.chrs) || ! (target %in% my.chrs)) {
				for (x in 1:length(is.okay)) { 
					if (is.okay[x]) { totals[x] <- totals[x] + nrow(.getPairs(files[x], anchor, target)) }
				}
				next
			}

			# Otherwise, collating the marginal counts.
			pulled <- list()
			for (x in 1:length(is.okay)) { 
                if (!is.okay[x]) {
					pulled[[x]] <- list(integer(0), integer(0), integer(0))
					next
				}
				stuff <- .getPairs(files[x], anchor, target)
				totals[x] <- totals[x] + nrow(stuff)
				ar <- .forgeInterval(stuff$anchor.pos, ext=ext, spacing=1L)
				tr <- .forgeInterval(stuff$target.pos, ext=ext, spacing=1L) 
                
				# Performing overlaps between the requested chromosomes. 
				alap <- findOverlaps(IRanges(ar$start, ar$end-1L), gr[[anchor]])
				tlap <- findOverlaps(IRanges(tr$start, tr$end-1L), gr[[target]])

				# Computing the support for each pair of regions.
				collected <- .Call(cxx_count_margins, 
					queryHits(alap), subjectHits(alap), 
					queryHits(tlap), subjectHits(tlap),
					anchor==target)
				if (is.character(collected)) { stop(collected) }

				my.anchors <- oridex[[anchor]][ collected[[1]][ collected[[2]] ] ]
				all.margins[my.anchors,x] <- all.margins[my.anchors,x] + collected[[3]][ collected[[2]] ]
				my.targets <- oridex[[target]][ collected[[1]][ !collected[[2]] ] ]
				all.margins[my.targets,x] <- all.margins[my.targets,x] + collected[[3]][ !collected[[2]] ]
			}
		}
	}

	# Cleaning up and cashing out.
	return(list(counts=all.margins, totals=totals, files=files))
}

# It's worth pointing out that the anchor/target definition is based on the 
# inbuild chromosome sorting order (i.e., later chromosomes will be the anchor)
# and, for pairs on the same chromosome, the sorting order of the start
# points of each region (i.e., later start points will be the anchor).
# This can be a bit weird if the supplied regions are not ordered, as 
# the anchor index can then be lower than the target index.
