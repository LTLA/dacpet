recountPET <- function(dirs, regions, ext=1L, filter=20L)
# This counts the number of PETs with each end overlapping one of 
# the specified regions. It then reports counts for each pair of
# regions with counts above the specified threshold.
#
# written by Aaron Lun
# 23 January, 2014
{
    nlibs<-length(dirs)
    if (!is.integer(ext)) { ext<-as.integer(ext) }
    if (!is.integer(filter)) { filter<-as.integer(filter) }
    if (nlibs==0) {
        stop("number of libraries must be positive")
    } else if (ext < 0) {
        stop("fragment length must be a non-negative integer")
    }
	o <- GenomicRanges::order(regions)
	sregions <- regions[o]
	gr <- split(ranges(sregions), seqnames(sregions))
	oridex <- split(o, seqnames(sregions))
	chrs <- names(gr)

    # Running through each pair of chromosomes.
    overall <- .loadIndices(dirs)
    all.coords <- list()
	all.counts <- list()
	ix <- 1L
	totals <- integer(length(dirs))

    for (anchor in names(overall)) {
        if (! (anchor %in% chrs) ) { next }
        current<-overall[[anchor]]
        for (target in names(current)) {
			if (! (target %in% chrs)) { next }
            fnames<-current[[target]]
            
			pulled <- list()
			for (x in 1:length(fnames)) { 
                if (!nchar(fnames[x])) {
					pulled[[x]] <- list(integer(0), integer(0), integer(0))
					next
				}
				stuff<-read.table(file.path(dirs[x], fnames[x]), header=TRUE, colClasses="integer", comment.char="")
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
			all.coords[[ix]] <- data.frame(oridex[[anchor]][combined[[1]]], oridex[[target]][combined[[2]]])
			ix <- ix + 1L
		}
	}

	# Cleaning up and cashing out.
	all.counts <- do.call(rbind, all.counts)
	all.coords <- do.call(rbind, all.coords)
	colnames(all.coords) <- c("anchor", "target")
	return(list(counts=all.counts, pairs=all.coords, totals=totals, region=regions))
}

# It's worth pointing out that the anchor/target definition is based on the 
# inbuild chromosome sorting order (i.e., later chromosomes will be the anchor)
# and, for pairs on the same chromosome, the sorting order of the start
# points of each region (i.e., later start points will be the anchor).
# This can be a bit weird if the supplied regions are not ordered, as 
# the anchor index can then be lower than the target index.
