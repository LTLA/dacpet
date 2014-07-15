####################################################################################################
# This tests the counting for known pairs.

suppressPackageStartupMessages(require(dacpet))

chromos <- c(chrA=1000, chrB=10000, chrC=100000)
simgen <- function(dir, num) {
	write.table(file=file.path(dir, "lengths.txt"), data.frame(chr=names(chromos), length=chromos), row.names=FALSE, sep="\t", quote=FALSE)
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

reggen <- function(num, width) {
	output <- GRanges()
	for (x in names(chromos)) {
		starts <- runif(num, 1, chromos[[x]])
		widths <- runif(num, width[1], width[2])
		suppressWarnings(output <- c(output, GRanges(x, IRanges(starts, starts+widths))))
	}
	return(output)
}

countcomp <- function(alldirs, regs, ext, filter=1L) {
	observed <- recountPET(alldirs, regs, ext=ext, filter=1L)
	observed.interact <- paste0(observed$pairs$anchor, ".", observed$pairs$target)
	dummycount <- observed$counts
	o <- order(regs)
	sregs <- regs[o]

	# Picking out the truth.
	overall <- dacpet:::.loadIndices(alldirs)
	for (x in names(overall)) {
		cur1 <- overall[[x]]
		for (y in names(cur1)) { 
			cur2 <- cur1[[y]]

			collected <- list()
			for (z in 1:length(alldirs)) {
                stuff <- read.table(file.path(alldirs[z], cur2[z]), header = TRUE)
				astart <- ifelse(stuff$anchor.pos > 0, stuff$anchor.pos - ext + 1L, -stuff$anchor.pos)
				aend <- astart + ext - 1L
				alap <- findOverlaps(GRanges(x, IRanges(astart, aend)), sregs)

				tstart <- ifelse(stuff$target.pos > 0, stuff$target.pos - ext + 1L, -stuff$target.pos)
				tend <- tstart + ext - 1L
				tlap <- findOverlaps(GRanges(y, IRanges(tstart, tend)), sregs)

				# Collating them to identify combinations.
				a.ok <- split(subjectHits(alap), queryHits(alap))
				t.ok <- split(subjectHits(tlap), queryHits(tlap))
				both.ok <- intersect(names(a.ok), names(t.ok))
				all.results <- list()
				for (pet in both.ok) {
					all.combos <- merge(a.ok[[pet]], t.ok[[pet]])
					interactions <- paste0(o[pmax(all.combos[,1], all.combos[,2])], ".",
							o[pmin(all.combos[,1], all.combos[,2])])
					all.results[[pet]] <- unique(interactions)
				}
				all.results <- unlist(all.results, use.names=FALSE)
				
				# Comparing them directly.
				final.counts <- table(all.results)
				comp <- match(names(final.counts), observed.interact)
				if (any(is.na(comp))) { stop("interaction not in canonical set") }
				if (!identical(as.integer(final.counts), observed$counts[comp,z])) { 
					stop("mismatches in count values") }

				# Getting rid of them after loading.
				dummycount[comp,z] <- dummycount[comp,z] - as.integer(final.counts)
			}
		}
	}

	# Checking that there isn't any unID's odds and ends.
	if (any(dummycount!=0L)) { 
		stop("interaction in observed set not identified, or identified multiple times") 
	}

	# Checking what happens when I set the filter on.
	filtered <- recountPET(alldirs, regs, ext=ext, filter=filter)
	keep <- rowSums(observed$counts) >= filter

	fpairs <- observed$pairs[keep,,drop=FALSE]
	rownames(fpairs) <- NULL
	fcounts <- observed$counts[keep,,drop=FALSE]

	if (!identical(fpairs, filtered$pairs) || !identical(fcounts, filtered$counts)) { 
		stop("filtering is not correct") }

	return(head(observed$pairs))
}

####################################################################################################
# Initializing the analysis.

set.seed(485632481)

dir1<-"blah.out.1"
dir2<-"blah.out.2"
dir3<-"blah.out.3"
dir.create(dir1)
dir.create(dir2)
dir.create(dir3)

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

####################################################################################################

unlink(dir1, recursive=TRUE)
unlink(dir2, recursive=TRUE)
unlink(dir3, recursive=TRUE)

####################################################################################################
# End.
