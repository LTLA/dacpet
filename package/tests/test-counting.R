###################################################################################################
# This tests the interaction counting capabilities of, well, the interaction counter. 

set.seed(1001);
suppressPackageStartupMessages(require(dacpet))
chromos<-c(chrA=500, chrB=300)

simgen <- function(dir, num) {
	write.table(file=file.path(dir, "lengths.txt"), data.frame(chr=names(chromos), length=chromos), 
			row.names=FALSE, sep="\t", quote=FALSE)
	log<-file.path(dir, "index.txt")
	if (file.exists(log)) { unlink(log) }
	counter<-1L
	for (i in 1:length(chromos)) { 
		max.anchor<-chromos[[i]];
		for (j in 1:i) {
			max.target<-chromos[[j]];
			anchors<-as.integer(floor(runif(num, 1, max.anchor)));
			targets<-as.integer(floor(runif(num, 1, max.target)));
			if (i==j){
				anchor.1<-pmax(anchors, targets);
				target.1<-pmin(anchors, targets);
				anchors<-anchor.1;
				targets<-target.1;
			}
			astr<-rbinom(num, 1, 0.5)==1
			tstr<-rbinom(num, 1, 0.5)==1
			anchors[astr]<--anchors[astr]
			targets[tstr]<--targets[tstr]

			basic<-paste0(counter, ".gz")
			fname<-file.path(dir, basic)
			write.table(file=fname, data.frame(anchors, targets), row.names=FALSE,
					col.names=c("anchor.pos", "target.pos"), quote=FALSE, sep="\t")
			write.table(file=log, data.frame(names(chromos)[i], names(chromos)[j], basic),
					row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE, append=TRUE)
			counter<-counter+1L
		}
	}
	return(invisible())
}

getInterval <- function(pt, ext, left=0, right=0, maxed=NULL) {
	if (pt<0) {
		start<--1L*pt
		end<-start+ext - 1L
	} else {
		end<-pt
		start<-end -ext + 1L
	}
	start<-start-right
	end<-end+left
	start<-pmax(start, 1L)
	if (!is.null(maxed)) { end<-pmin(end, maxed) }
	return(c(start, end))
}

# We set up the comparison function to check our results. 

comp <- function(dir1, dir2, ext, spacing=10, left=0, right=0, filter=1L) {
	proposed<-countPET(dirs=c(dir1, dir2), ext=ext, shift=left, width=left+right+1, filter=filter, spacing=spacing)# restrict="chrA")
	stopifnot(all(rowSums(proposed$count)[proposed$pair$index] >= filter))

	# We check whether the regions make sense.
	checker <- list()
	for (i in 1:length(chromos)) {
		spaced<-0:ceiling(chromos[[i]]/spacing)*spacing+1L
		checker[[i]] <- GRanges(names(chromos)[i], IRanges(spaced-left, spaced+right), gunk=spaced)
	}
	suppressWarnings({
		checker <- do.call(c, checker)
		seqlengths(checker) <- chromos
		checker <- trim(checker)
		checker <- checker[width(checker) > 0L]
	})
	space.pts <- split(checker$gunk, as.character(seqnames(checker)))
	last.off <- 0L
	offsets <- list()
	for (i in 1:length(chromos)) { 
		cur.chr <- names(chromos)[i]
		offsets[[cur.chr]] <- last.off
		last.off <- length(space.pts[[cur.chr]]) + last.off
	}
	checker$gunk <- NULL
	if (!identical(checker, proposed$region)) { stop("mismatch in proposed regions") }
	
	# We need to determine who's who.
	x1<-read.table(file.path(dacpet:::.getIndex(dir1)))
	x2<-read.table(file.path(dacpet:::.getIndex(dir2)))
	for (k in 1:length(chromos)) {
		cur.k<-names(chromos)[k]
		for (l in 1:k) {
			cur.l<-names(chromos)[l]

			# Loading counts.
			x<-list(read.table(file.path(dir1, x1[x1[,1]==cur.k & x1[,2]==cur.l,3]), header=TRUE),
				read.table(file.path(dir2, x2[x2[,1]==cur.k & x2[,2]==cur.l,3]), header=TRUE))
			max.anchor<-chromos[k]
			aspace<-space.pts[[k]]
			aoff<-offsets[[k]]
			max.target<-chromos[l]
			tspace<-space.pts[[l]]
			toff<-offsets[[l]]

			original<-list()
			other<-0
			for (g in 1:length(x)) {
				mat1<-matrix(0L, nrow=max.anchor+left, ncol=max.target+left) # Need to +left, as 'space.pts' can exceed chromosome length.
				for (i in 1:nrow(x[[g]])) {
					arange<-getInterval(x[[g]][i,1], ext, left=left, right=right, maxed=max.anchor+left)
					if (arange[2] < arange[1]) { next }
					arange<-arange[1]:arange[2] 
					trange<-getInterval(x[[g]][i,2], ext, left=left, right=right, maxed=max.target+left)
					if (trange[2] < trange[1]) { next }
					trange<-trange[1]:trange[2]

					if (k!=l) {
						mat1[arange,trange]<-mat1[arange,trange]+1L
					} else {
						# Reflecting around the diagonal for intra-chromosomals.
						collected <- unique(c(outer((arange-1L)*nrow(mat1), trange, FUN="+"), 
											outer((trange-1L)*nrow(mat1), arange, FUN="+")))
						mat1[collected] <- mat1[collected] + 1L
					}
				}
				submat1<-mat1[aspace,tspace]
				if (k==l) { submat1[upper.tri(submat1)] <- 0L }
				original[[g]]<-submat1

				# Subtracting off the elements that we know.
				wascovered<-matrix(FALSE, nrow(submat1), ncol(submat1))
				keep<-as.logical(seqnames(proposed$region[proposed$pair$anchor])==cur.k & 
						seqnames(proposed$region[proposed$pair$target])==cur.l)
				if (any(keep)){ 
					kept<-proposed$pairs[keep,,drop=FALSE]
					kept.counts <- proposed$counts[keep,,drop=FALSE]
					for (i in 1:nrow(kept)) {
						arange<-kept$anchor[i] - aoff
						trange<-kept$target[i] - toff
 			   			submat1[arange,trange]<-submat1[arange,trange]-kept.counts[i,g] 
						if (any(wascovered[arange,trange])) { stop("overlapping boxes reported in output") }
						wascovered[arange,trange]<-TRUE
						
						if (any(submat1[arange,trange]<0L)) {
							stop("additional counts present in the proposed set which are not in the truth") 
						}
					}
				}
				other<-other+submat1
			}
			stopifnot(all(other < filter))  # Checks for any rows above the filter in the truth that are not in the proposed.
		}
	}
	return(head(proposed$pairs))
}

###################################################################################################
# Checking a vanilla count.

dir1<-"temp.out.1"
dir2<-"temp.out.2"
dir.create(dir1)
dir.create(dir2)

simgen(dir1, 20)
simgen(dir2, 10)
comp(dir1, dir2, ext=10)
comp(dir1, dir2, ext=10, left=9)
comp(dir1, dir2, ext=10, right=25)
comp(dir1, dir2, ext=10, spacing=23)
comp(dir1, dir2, ext=10, filter=2)
comp(dir1, dir2, ext=10, filter=5)
comp(dir1, dir2, ext=55)
comp(dir1, dir2, ext=55, left=15, spacing=30)
comp(dir1, dir2, ext=55, right=25)
comp(dir1, dir2, ext=55, spacing=46)
comp(dir1, dir2, ext=55, filter=5)
comp(dir1, dir2, ext=55, filter=20)

# Throwing more sequences in.
simgen(dir1, 50)
simgen(dir2, 30)
comp(dir1, dir2, ext=10)
comp(dir1, dir2, ext=10, left=12, spacing=20)
comp(dir1, dir2, ext=10, right=15)
comp(dir1, dir2, ext=10, spacing=13)
comp(dir1, dir2, ext=10, filter=5)
comp(dir1, dir2, ext=10, filter=10)
comp(dir1, dir2, ext=55)
comp(dir1, dir2, ext=55, left=35, spacing=50)
comp(dir1, dir2, ext=55, right=45)
comp(dir1, dir2, ext=55, spacing=36)
comp(dir1, dir2, ext=55, filter=15)
comp(dir1, dir2, ext=55, filter=30)

# Throwing a ridiculous number of sequences in.
simgen(dir1, 100)
simgen(dir2, 200)
comp(dir1, dir2, ext=10)
comp(dir1, dir2, ext=10, left=12, spacing=13)
comp(dir1, dir2, ext=10, right=15)
comp(dir1, dir2, ext=10, spacing=13)
comp(dir1, dir2, ext=10, filter=5)
comp(dir1, dir2, ext=10, filter=10)
comp(dir1, dir2, ext=55)
comp(dir1, dir2, ext=55, left=35, spacing=38)
comp(dir1, dir2, ext=55, right=45)
comp(dir1, dir2, ext=55, spacing=36)
comp(dir1, dir2, ext=55, filter=15)
comp(dir1, dir2, ext=55, filter=30)

# Easing off the pedal, and fiddling with some parameters
simgen(dir1, 10)
simgen(dir2, 20)
comp(dir1, dir2, ext=20)
comp(dir1, dir2, ext=23, left=12, spacing=29)
comp(dir1, dir2, ext=23, right=17)
comp(dir1, dir2, ext=23, spacing=18)
comp(dir1, dir2, ext=23, filter=2)
comp(dir1, dir2, ext=20, filter=5)
comp(dir1, dir2, ext=41)
comp(dir1, dir2, ext=41, left=9)
comp(dir1, dir2, ext=41, right=58)
comp(dir1, dir2, ext=41, spacing=26)
comp(dir1, dir2, ext=41, filter=3)
comp(dir1, dir2, ext=41, filter=10)

# Throwing another ridiculous number of sequences in.
simgen(dir1, 100)
simgen(dir2, 200)
comp(dir1, dir2, ext=10, left=8, right=-5)
comp(dir1, dir2, ext=10, left=12, spacing=13, right=-10)
comp(dir1, dir2, ext=10, right=-5, left=8)
comp(dir1, dir2, ext=10, spacing=13, left=5, right=-2)
comp(dir1, dir2, ext=10, filter=5, left=8, right=-1)
comp(dir1, dir2, ext=10, filter=10, left=8, right=-1)
comp(dir1, dir2, ext=55, left=6, right=-3)
comp(dir1, dir2, ext=55, left=35, spacing=40, right=-35)
comp(dir1, dir2, ext=55, right=-5, left=5, spacing=20)
comp(dir1, dir2, ext=55, spacing=36, left=10, right=-10)

##################################################################################################
# Cleaning up.

unlink(dir1, recursive=TRUE)
unlink(dir2, recursive=TRUE)

##################################################################################################
# End.
