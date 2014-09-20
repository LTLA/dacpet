####################################################################################################
# This tests the counting for known pairs.

suppressPackageStartupMessages(require(dacpet))
suppressPackageStartupMessages(require(rhdf5))

chromos <- c(chrA=1000, chrB=10000, chrC=100000)
simgen <- function(file, num) {
	if (file.exists(file)) { unlink(file) }
	h5createFile(file)
	h5createGroup(file, "counts")
	h5write(data.frame(chr=names(chromos), length=chromos, stringsAsFactors=FALSE), file, 'lengths')
    
	for (i in 1:length(chromos)) { 
        max.anchor<-chromos[[i]]
		cur.anchor <- names(chromos)[i]
		h5createGroup(file, file.path("counts", cur.anchor))

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
            astr<-rbinom(num, 1, 0.5)==1L
            tstr<-rbinom(num, 1, 0.5)==1L
            anchors[astr]<--anchors[astr]
            targets[tstr]<--targets[tstr]
			
			h5write(data.frame(anchor.pos=anchors, target.pos=targets), file,
				file.path("counts", cur.anchor, names(chromos)[j]))
        }
    }
    return(invisible())
}

reggen <- function(num, width) {
	output <- GRanges()
	for (x in names(chromos)) {
		starts <- runif(num, 1, chromos[[x]])
		widths <- runif(num, width[1], width[2])
		suppressWarnings(output <- c(output, GRanges(x, IRanges(starts, starts+widths))))
	}
	return(output)
}

countcomp <- function(alldirs, regs, ext, filter=1L, restrict=NULL) {
	observed <- recountPET(alldirs, regs, ext=ext, filter=1L, restrict=restrict)
	observed.interact <- paste0(anchors(observed, id=TRUE), ".", targets(observed, id=TRUE))
	dummycount <- counts(observed)
	dummytotes <- info(observed)$totals
	o <- order(regs)
	stopifnot(identical(o, regions(observed)$original))
	regs <- regs[o]

	# Picking out the truth.
	overall <- dacpet:::.loadIndices(alldirs)
	for (x in names(overall)) {
		cur1 <- overall[[x]]
		for (y in names(cur1)) { 
			cur2 <- cur1[[y]]
			if (!is.null(restrict) && !(x %in% restrict && y %in% restrict)) { next }

			collected <- list()
			for (z in 1:length(alldirs)) {
				if (cur2[z]) { 
            		stuff <- dacpet:::.getPairs(alldirs[z], x, y)
					attributes(stuff$anchor.pos) <- attributes(stuff$target.pos) <- NULL
				} else {
					stuff <- data.frame(anchor.pos=integer(0), target.pos=integer(0))
				}
				astart <- ifelse(stuff$anchor.pos > 0, stuff$anchor.pos - ext + 1L, -stuff$anchor.pos)
				aend <- astart + ext - 1L
				alap <- findOverlaps(GRanges(x, IRanges(astart, aend)), regs)

				tstart <- ifelse(stuff$target.pos > 0, stuff$target.pos - ext + 1L, -stuff$target.pos)
				tend <- tstart + ext - 1L
				tlap <- findOverlaps(GRanges(y, IRanges(tstart, tend)), regs)

				# Collating them to identify combinations.
				a.ok <- split(subjectHits(alap), queryHits(alap))
				t.ok <- split(subjectHits(tlap), queryHits(tlap))
				both.ok <- intersect(names(a.ok), names(t.ok))
				all.results <- list()
				for (pet in both.ok) {
					all.combos <- merge(a.ok[[pet]], t.ok[[pet]])
					interactions <- paste0(pmax(all.combos[,1], all.combos[,2]), ".",
							pmin(all.combos[,1], all.combos[,2]))
					all.results[[pet]] <- unique(interactions)
				}
				all.results <- unlist(all.results, use.names=FALSE)
				
				# Comparing them directly.
				final.counts <- table(all.results)
				comp <- match(names(final.counts), observed.interact)
				if (any(is.na(comp))) { stop("interaction not in canonical set") }
				if (!identical(as.integer(final.counts), counts(observed)[comp,z])) { 
					stop("mismatches in count values") }

				# Getting rid of them after loading.
				dummycount[comp,z] <- dummycount[comp,z] - as.integer(final.counts)
				dummytotes[z] <- dummytotes[z] - length(astart)
			}
		}
	}

	# Checking that there isn't any unID's odds and ends.
	if (any(dummycount!=0L)) { 
		stop("interaction in observed set not identified, or identified multiple times") 
	}
	if (any(dummytotes!=0L)) { 
		stop("totals don't match up")
	}

	# Checking what happens when I set the filter on.
	filtered <- recountPET(alldirs, regs, ext=ext, filter=filter, restrict=restrict)
	keep <- rowSums(counts(observed)) >= filter

	if (!identical(anchors(observed[keep,], id=TRUE), anchors(filtered, id=TRUE)) ||
		!identical(targets(observed[keep,], id=TRUE), targets(filtered, id=TRUE)) ||
		!identical(counts(observed[keep,]), counts(filtered))) {
		stop("filtering is not correct") }

	return(head(data.frame(anchor=anchors(observed, id=TRUE), target=targets(observed, id=TRUE))))
}

####################################################################################################
# Initializing the analysis.

set.seed(485632481)

dir.create("temp-pair")
dir1<-"temp-pair/out.1.h5"
dir2<-"temp-pair/out.2.h5"
dir3<-"temp-pair/out.3.h5"

simgen(dir1, 200)
simgen(dir2, 100)

myreg <- reggen(10, c(10, 50))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

myreg <- reggen(10, c(40, 100))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

myreg <- reggen(50, c(10, 25))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

# Again, with more reads.

simgen(dir1, 500)
simgen(dir2, 1000)

myreg <- reggen(10, c(10, 50))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

myreg <- reggen(10, c(40, 100))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

myreg <- reggen(50, c(10, 25))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

# Again, with more directories.

simgen(dir1, 200)
simgen(dir2, 100)
simgen(dir3, 400)

myreg <- reggen(10, c(10, 50))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

myreg <- reggen(10, c(40, 100))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

myreg <- reggen(50, c(10, 25))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

# Lastly, with a lot more regions.

simgen(dir1, 200)
simgen(dir2, 100)
simgen(dir3, 400)

myreg <- reggen(20, c(10, 50))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

myreg <- reggen(20, c(40, 100))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

myreg <- reggen(100, c(10, 25))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)
countcomp(c(dir1, dir2), myreg, 100, filter=5)

# Repeating with restriction.

myreg <- reggen(100, c(10, 25))
countcomp(c(dir1, dir2), myreg, 10, restrict="chrA")
countcomp(c(dir1, dir2), myreg, 20, restrict="chrB")
countcomp(c(dir1, dir2), myreg, 30, restrict=c("chrB", "chrC"))
countcomp(c(dir1, dir2), myreg, 100, filter=5, restrict=c("chrA", "chrC"))

####################################################################################################

unlink("temp-pair", recursive=TRUE)

####################################################################################################
# End.
