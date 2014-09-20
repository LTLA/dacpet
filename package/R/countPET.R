countPET <- function(files, ext=1L, shift=0L, width=5000L, spacing=width, filter=20L, restrict=NULL) 
# This function collates counts across multiple experiments to get the full set of results. 
# It takes a set of directories which it loads in and compiles into a set of weighted counts 
# for all sliding windows (or mats, in interaction space), in preparation for further analysis.
#
# written by Aaron Lun
# Created 15 January 2014. Lastmodified 19 September 2014.
{
	nlibs <- length(files)
	left <- shift
	right <- width - shift - 1L
	if (!is.integer(spacing)) { spacing <- as.integer(spacing) }
	if (!is.integer(ext)) { ext <- as.integer(ext) }
	if (!is.integer(left)) { left <- as.integer(left) }
	if (!is.integer(right)) { right <- as.integer(right) }
	if (!is.integer(filter)) { filter <- as.integer(filter) }
	if (nlibs==0L) { 
		stop("number of libraries must be positive")
	} else if (spacing < 1L) { 
 	   stop("spacing must be a positive integer")
   	} else if (ext < 0L) { 
		stop("fragment length must be a non-negative integer")
	} else if (shift < 0L) { 
		stop("left window shift must be a non-positive integer")
	} else if (width <= 0L) { 
		stop("width must be a positive integer")
	}
	chromosomes <- .getChrs(files)
	chrs <- names(chromosomes)

	# Checking out the incoming left/right, to decide whether we need a point at the start.
	if (left >= spacing) { stop("shift must be less than the spacing") }
   	at.start <- right >= 0L
	first.pt <- ifelse(at.start, 1L, spacing + 1L)

	# Setting up structures to store results.
    ptr <- .Call(cxx_create_counts, nlibs)
    if (is.character(ptr)) { stop(ptr) }
	full.sizes <- rep(0L, nlibs)

	# Constructing a GRanges object of the windows corresponding to each point.  We
	# need to know the number of points based on left/right shifts. Also noting the
	# offset for each chromosome, so we can merge all coordinates into a single object. 
	windows <- list()
	offsets <- list()
	npts <- list()
	last.hit <- 0L
	for (chr in names(chromosomes)) {	
		offsets[[chr]]<-last.hit
		cur.len <- chromosomes[[chr]]
        at.end <- spacing - (cur.len - 1L) %% spacing <= left # Figuring out if the shift adds something extra on the end.
		cur.npt <- at.start + as.integer((cur.len - 1L)/spacing) + at.end
		npts[[chr]] <- cur.npt

		if (cur.npt) {
        	center<-(0:(cur.npt - 1L))*spacing+first.pt
			current<-GRanges(chr, IRanges(pmax(1L, as.integer(center-left)), 
						pmin(chromosomes[[chr]], as.integer(center+right))))
			windows[[chr]]<-current
			last.hit<-last.hit+length(current)
		}
	}
	names(windows)<-NULL
	windows<-suppressWarnings(do.call(c, windows))
	seqlengths(windows)<-chromosomes

	# Running through each pair of chromosomes.
    overall <- .loadIndices(files)
	all.anchors <- all.targets <- coldex <- list()
	counter <- 1L

	for (anchor in names(overall)) {
		stopifnot(anchor %in% chrs)
		if (!is.null(restrict) && !(anchor %in% restrict)) { next }
		current <- overall[[anchor]]
        for (target in names(current)) {
			stopifnot(target %in% chrs)
			if (!is.null(restrict) && !(target %in% restrict)) { next }
			is.present <- current[[target]]

			pulled<-list()
			for (x in 1:nlibs) {
				if (!is.present[x]) { 
					pulled[[x]] <- data.frame(integer(0), integer(0), integer(0), integer(0))					
					next 
				}
				stuff <- .getPairs(files[x], anchor, target)
				full.sizes[x] <- full.sizes[x]+nrow(stuff)
				arange <- .forgeInterval(stuff$anchor.pos, ext=ext, spacing=spacing, left=left, right=right, maxed=npts[[anchor]] + 1L, is.first=at.start)
				trange <- .forgeInterval(stuff$target.pos, ext=ext, spacing=spacing, left=left, right=right, maxed=npts[[target]] + 1L, is.first=at.start)

				if (anchor==target) { 
					# Need to check who's the anchor, and who's the target, for intrachromosomal
					# read pairs after addition/subtraction of fragment sizes. If necessary, some
					# are flipped around; otherwise, the bulk of the area in the interaction space
					# corresponding to the extended reads might lie past the diagonal and thus be
					# inaccessible for some window pairs following the anchor>=target rule. Note
					# that these areas are squares, so it's sufficient to just switch each pair of
					# extended reads to cover the largest area under the diagonal; you don't need
					# to collect counts for both switched and unswitched read pairs, as it's 
					# impossible for a window pair following the anchor>=target rule to overlap one 
					# and not the other (as the window must lie on the diagonal, at its most extreme).
				
					a.as.anchor<-arange$start >= trange$start & arange$end >= trange$end
					if (!all(a.as.anchor)) { 
						ax<-ifelse(a.as.anchor, arange$start, trange$start)
						tx<-ifelse(a.as.anchor, trange$start, arange$start)
						arange$start<-ax
						trange$start<-tx
						ae<-ifelse(a.as.anchor, arange$end, trange$end)
						te<-ifelse(a.as.anchor, trange$end, arange$end)
						arange$end<-ae
						trange$end<-te
					}
					stopifnot(all(arange$start >= trange$start & arange$end >= trange$end))
				} 
							
				# Getting rid of useless combinations, reordering the remainders.
				combined<-data.frame(arange$start, trange$start, arange$end, trange$end)
				combined<-combined[combined[,1]!=combined[,3] & combined[,2]!=combined[,4],]
				combined<-combined[do.call(order, combined),]
				pulled[[x]]<-combined
			}

			# Collating the results.
            out<-.Call(cxx_count_chia, ptr, pulled, filter, anchor==target)
			if (is.character(out)) { stop(out) }
			out[[1]] <- out[[1]] + offsets[[anchor]]
			out[[2]] <- out[[2]] + offsets[[target]]
			if (offsets[[anchor]] >= offsets[[target]]) {
				all.anchors[[counter]] <- out[[1]]
 				all.targets[[counter]] <- out[[2]]
			} else {
				all.anchors[[counter]] <- out[[2]]
 				all.targets[[counter]] <- out[[1]]
			}
			coldex[[counter]] <- out[[3]]+1L
			counter <- counter + 1L
		}
	}

	# Storing the output.
	counts<-.Call(cxx_get_counts, ptr)
	if (is.character(counts)) { stop(counts) }
	return(IList(counts=counts[unlist(coldex),,drop=FALSE],  
		info=data.frame(totals=full.sizes, files=files), anchors=unlist(all.anchors), 
		targets=unlist(all.targets), regions=windows))
}

.forgeInterval <- function(pt, ext, spacing, left=0, right=0, maxed=NULL, is.first=TRUE) 
# This extends each tag BACKWARDs to the specified 'ext', as the bimodality is
# theoretically flipped around for ChIA-PET data (e.g., outward patterns
# observed for a binding site).  It also extends them to the left and right by
# 'right' and 'left', respectively, before division by the spacing to obtain
# subpoints. Capping is done to block all points above 'maxed', as there's no
# internal checks to protect against super-large chr's.
{
	# Adjusting coordinates so that extension recovers the correct values.
	is.pos <- pt > 0L
	pt[!is.pos] <- -1L*pt[!is.pos]
	pt[is.pos] <- pt[is.pos] - ext + 1L
	ender <- pt + ext + left
	pt <- pmax(pt - right, 1L)

	# Getting the actual coordinates after adjusting for spacing.
	pt <- ifelse(pt==1L, 1L, as.integer((pt-2L)/spacing)+1L+is.first)
	ender <- ifelse(ender==1L, 1L, as.integer((ender-2L)/spacing)+1L+is.first)
	if (!is.null(maxed)) { 
		pt <- pmin(pt, maxed) 
		ender <- pmin(ender, maxed)
	}
	return(list(start=pt, end=ender))
}

.getChrs <- function(files) {
	# Checking genome consistency.
	chromosomes <- NULL
	for (f in files) {
		temp <- .getLengths(f)
		temp2 <- as.vector(temp$length)
		names(temp2) <- temp$chr
		temp2 <- temp2[order(names(temp2))]
		if (is.null(chromosomes)) { chromosomes <- temp2 }
		else if (!identical(temp2, chromosomes)) { 
			stop("chromosome identities and lengths should be the same across files")
		}
	}
	return(chromosomes)
}
