recountPET <- function(files, regions, ext=1L, filter=20L, restrict=NULL)
# This counts the number of PETs with each end overlapping one of 
# the specified regions. It then reports counts for each pair of
# regions with counts above the specified threshold.
#
# written by Aaron Lun
# Created 23 January 2014. Last modified 19 September 2014.
{
    nlibs<-length(files)
    if (!is.integer(ext)) { ext<-as.integer(ext) }
    if (!is.integer(filter)) { filter<-as.integer(filter) }
    if (nlibs==0) {
        stop("number of libraries must be positive")
    } else if (ext < 0) {
        stop("fragment length must be a non-negative integer")
    }
	o <- GenomicRanges::order(regions)
	regions <- regions[o]

	# Setting up chromosome-by-chromosome data.
	gr <- offsets <- list()
	rv <- runValue(seqnames(regions))
	rl <- runLength(seqnames(regions))
	last.hit <- 0L
	for (x in 1:length(rv)) { 
		chr <- as.character(rv[x])
		gr[[chr]] <- ranges(regions)[last.hit + 1:rl[x]]
		offsets[[chr]] <- last.hit
		last.hit <- last.hit + rl[x]
	}
	chrs <- names(.getChrs(files))
	my.chrs <- names(gr)

    # Running through each pair of chromosomes.
    overall <- .loadIndices(files)
    all.anchors <- all.targets <- list()
	all.counts <- list()
	ix <- 1L
	totals <- integer(nlibs)

    for (anchor in names(overall)) {
		stopifnot(anchor %in% chrs)
		if (!is.null(restrict) && !(anchor %in% restrict)) { next }
        current<-overall[[anchor]]
        for (target in names(current)) {
			stopifnot(target %in% chrs)
			if (!is.null(restrict) && !(target %in% restrict)) { next }
            is.okay <- current[[target]]

			# Deciding whether to skip counting (but we still need the totals).
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
				collected <- .Call(cxx_count_all_pairs,
					queryHits(alap), subjectHits(alap), 
					queryHits(tlap), subjectHits(tlap),
					anchor==target)
				if (is.character(collected)) { stop(collected) }
				pulled[[x]] <- collected
			}
			
			# Combining counts between libraries, with some filtering to reduce memory usage.
			combined <- .Call(cxx_aggregate_pair_counts, pulled, filter)
			if (is.character(combined)) { stop(combined) }
			all.counts[[ix]] <- combined[[3]]

			combined[[1]] <- combined[[1]] + offsets[[anchor]]
			combined[[2]] <- combined[[2]] + offsets[[target]]
			if (offsets[[anchor]] >= offsets[[target]]) {  
				all.anchors[[ix]] <- combined[[1]]
 		   		all.targets[[ix]] <- combined[[2]]
			} else {
				all.anchors[[ix]] <- combined[[2]]
 		   		all.targets[[ix]] <- combined[[1]]
			}
			ix <- ix + 1L
		}
	}

	# Cleaning up and cashing out.
	all.counts <- do.call(rbind, all.counts)
	regions$original <- o
	return(IList(counts=all.counts, info=data.frame(totals=totals, files=files), 
		anchors=unlist(all.anchors), targets=unlist(all.targets), regions=regions))
}
