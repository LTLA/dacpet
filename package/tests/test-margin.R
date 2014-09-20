####################################################################################################
# This tests the marginal counting for supplied regions.

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

countcomp <- function(alldirs, regs, ext, restrict=NULL) {
	observed <- marginPET(alldirs, regs, ext=ext, restrict=restrict)
	dummycount <- counts(observed)
	dummytotes <- info(observed)$totals

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
				if (x==y) {
					all.s <- c(subjectHits(alap), subjectHits(tlap))
					all.q <- c(queryHits(alap), queryHits(tlap))
					o <- order(all.s, all.q)
					is.diff <- c(TRUE, diff(all.s[o])!=0L | diff(all.q[o])!=0L)
					dummycount[,z] <- dummycount[,z] - tabulate(all.s[o][is.diff], nbins=length(regs))
				} else {
					dummycount[,z] <- dummycount[,z] - tabulate(c(subjectHits(alap), subjectHits(tlap)), nbins=length(regs)) 
				}
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
	return(head(counts(observed)))
}

####################################################################################################
# Initializing the analysis.

set.seed(485632481)

dir.create("temp-marg")
dir1<-"temp-marg/out.1.h5"
dir2<-"temp-marg/out.2.h5"
dir3<-"temp-marg/out.3.h5"

simgen(dir1, 200)
simgen(dir2, 100)

myreg <- reggen(10, c(10, 50))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

myreg <- reggen(10, c(40, 100))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

myreg <- reggen(50, c(10, 25))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

# Again, with more reads.

simgen(dir1, 500)
simgen(dir2, 1000)

myreg <- reggen(10, c(10, 50))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

myreg <- reggen(10, c(40, 100))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

myreg <- reggen(50, c(10, 25))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

# Again, with more directories.

simgen(dir1, 200)
simgen(dir2, 100)
simgen(dir3, 400)

myreg <- reggen(10, c(10, 50))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

myreg <- reggen(10, c(40, 100))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

myreg <- reggen(50, c(10, 25))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

# Lastly, with a lot more regions.

simgen(dir1, 200)
simgen(dir2, 100)
simgen(dir3, 400)

myreg <- reggen(20, c(10, 50))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

myreg <- reggen(20, c(40, 100))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

myreg <- reggen(100, c(10, 25))
countcomp(c(dir1, dir2), myreg, 10)
countcomp(c(dir1, dir2), myreg, 100)

# Repeating with restriction.

myreg <- reggen(100, c(10, 25))
countcomp(c(dir1, dir2), myreg, 10, restrict="chrA")
countcomp(c(dir1, dir2), myreg, 20, restrict="chrB")
countcomp(c(dir1, dir2), myreg, 30, restrict=c("chrB", "chrC"))

####################################################################################################

unlink("temp-marg", recursive=TRUE)

####################################################################################################
# End.
