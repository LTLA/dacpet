preparePET<-function(bam, dir, dedup=TRUE, yield=1e7, minq=0)
# This function prepares ChIA-PET data by stripping out valid pairs from 
# the BAM file and storing a set of tables describing the interacting loci 
# (if it passes quality controls). We use an anchor/target set-up whereby 
# the later chromosome is designated as the anchor. 
#
# written by Aaron Lun
# 7 November, 2013
{
	# Running through all pairs.
	bf<-open(BamFile(bam, yieldSize=yield, obeyQname=TRUE, index=character(0)))
	chromosomes<-scanBamHeader(bf)[[1]]
	chrs<-names(chromosomes)
	if (file.exists(dir)) { unlink(dir, recursive=TRUE) }
	dir.create(dir)
	allfiles<-list()
	file.count<-1L

	# Storing some diagnostics.
	singles <- 0L
	multies <- 0L
	total <- 0L
	marked <- 0L
	filtered <- 0L
	valid <- 0L
	
	while (1) {
		out<-scanBam(bf, param=ScanBamParam(what=c("qname", "flag", "rname", "pos", "qwidth", "mapq")))[[1]]
#			flag=scanBamFlag(isUnmappedQuery=FALSE, isPaired=TRUE, hasUnmappedMate=FALSE)))[[1]]
		if (!length(out[[1]])) { break; }

		# Converting chromosome ids.
		rematched<-match(levels(out$rname), chrs)
		if (any(is.na(rematched))) { stop("unrecognised chromosomes in the BAM file") }
		cur.chrs<-rematched[as.integer(out$rname)]

		# Collating them into read pairs.
		name.rle <- rle(out$qname)
		is.single <- name.rle$lengths==1L
		is.other <- name.rle$lengths>2L
		is.unmapped <- bitwAnd(out$flag, 0x4)!=0L | out$mapq < minq
		is.marked <- bitwAnd(out$flag, 0x400)!=0L
		mate.read <- cumsum(name.rle$lengths)[!is.single & !is.other] 
		first.read <- mate.read - 1L

		# Assembling statistics.
		singles <- singles + sum(is.single)
		multies <- multies + sum(is.other)
		total <- total + sum(!is.single)
		either.marked <-  is.marked[first.read] | is.marked[mate.read]
		marked <- marked + sum(either.marked)
		either.unmapped <- is.unmapped[first.read] | is.unmapped[mate.read]
		filtered <- filtered + sum(either.unmapped)
		keep <- (!dedup | !either.marked) & !either.unmapped
		first.read <- first.read[keep]
		mate.read <- mate.read[keep]
		
		extra.valid <- sum(keep)
		if (!extra.valid) { next }
		valid <- valid + extra.valid

		# Getting first, second reads.
		is.reverse<-bitwAnd(out$flag, 0x10)!=0L
		chr1<-cur.chrs[first.read]
		pos1<-out$pos[first.read]
		str1<-is.reverse[first.read]
		if (!all(str1)) { pos1[!str1]<-pos1[!str1]+out$qwidth[first.read][!str1]-1L } # Storing 3' ends.
		chr2<-cur.chrs[mate.read]
		pos2<-out$pos[mate.read]
		str2<-is.reverse[mate.read]
		if (!all(str2)) { pos2[!str2]<-pos2[!str2]+out$qwidth[mate.read][!str2]-1L }

		# Now, actually figuring out who's the target, and who is the anchor.
		is1anchor<-chr1 > chr2
		equalities<-chr1==chr2
		if (any(equalities)) { is1anchor[equalities]<-pos1[equalities] > pos2[equalities] }
		anchor.chr<-ifelse(is1anchor, chr1, chr2)
		target.chr<-ifelse(is1anchor, chr2, chr1)

		pos1[str1]<--1L*pos1[str1]  # We mark reads on the reverse strand with '-1'.
		pos2[str2]<--1L*pos2[str2]
		anchor.pos<-ifelse(is1anchor, pos1, pos2)
		target.pos<-ifelse(is1anchor, pos2, pos1)

		# Ordering them and putting each pair into individual files.
		o<-order(anchor.chr, target.chr)
		anchor.chr<-anchor.chr[o]
		target.chr<-target.chr[o]
		anchor.pos<-anchor.pos[o]
		target.pos<-target.pos[o]
		is.diff<-which(c(TRUE, diff(anchor.chr)>0 | diff(target.chr)>0))
		total.runs<-length(is.diff)
		is.diff<-c(is.diff, length(o)+1L)

		# Dumping it into a set of temporary files; avoid multiple 'rbind' calls.
		for (i in 1:total.runs) {
			start<-is.diff[i]
			end<-is.diff[i+1L]-1L
			anchor<-chrs[anchor.chr[start]]
			target<-chrs[target.chr[start]]
			if (is.null(allfiles[[anchor]])) { allfiles[[anchor]]<-list() }
			
			current.file<-allfiles[[anchor]][[target]]	
			print.names<-FALSE
			if (is.null(current.file)) { 
				print.names<-c("anchor.pos", "target.pos")
				current.file<-file.path(dir, paste0(file.count, ".gz"))
				allfiles[[anchor]][[target]]<-current.file
				file.count<-file.count+1L
			}	
			fout<-gzfile(current.file, open="ab")
			write.table(file=fout, cbind(anchor.pos[start:end], target.pos[start:end]), 
				sep="\t", quote=FALSE, col.names=print.names, row.names=FALSE)
			close(fout)
		}	
	}
	close(bf)

	# We print out an index file for reference to the chromosomes corresponding to each file.
	indices<-lapply(names(allfiles), FUN=function(anchor) {
		tfiles<-allfiles[[anchor]]
		data.frame(anchor, names(tfiles), basename(unlist(tfiles)))
	})	
	names(indices)<-NULL
	indices<-do.call(rbind, indices)
	write.table(file=.getIndex(dir), indices, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
	
	# We also need to write in how long the sequences are.
	write.table(file=.getLengths(dir), data.frame(chr=chrs, length=chromosomes), 
		row.names=FALSE, sep="\t", quote=FALSE)	
	return(list(other=c(singles=singles, multi=multies), 
		pairs=c(total=total, marked=marked, filtered=filtered, valid=valid)));
}

.getIndex<- function(dir) { file.path(dir, "index.txt") }
.getLengths<- function(dir) { file.path(dir, "lengths.txt") }
.saveExt<-function(x, fname) {
    fopen<-gzfile(fname, open="wb")
	write.table(file=fopen, x, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	close(fopen)
}

