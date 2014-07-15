stripOutwardPET <- function(dir, out=NULL, min.gap=1e4, discard.to=NULL)
# This strips outward-facing read pairs with a gap below a minimum size. The
# idea is to get rid of read pairs which are formed by self-ligation.  This is
# results in lots of "interactions" which have a significant homo-/hetero-count
# fold change, which are really just self-ligated binding sites which may or
# may not interact. Loss of power is mitigated as we only remove outward facing
# pairs.
#
# written by Aaron Lun
# 22 January, 2014
{
	if (!is.null(out)) {
		if (file.exists(out)) { unlink(out, recursive=TRUE) }
 	    dir.create(out)
		allfiles <- dir(dir)
		file.copy(file.path(dir, allfiles), file.path(out, allfiles))
		dir <- out	
	}
	discarded <- 0L 
	if (!is.null(discard.to)) {
		handle <- file(discard.to, open="w")
		close(handle)
	}
	
	# Running through each chromosome.
    overall<-.loadIndices(dir)
	reindex <- NULL
    for (curchr in names(overall)) {
		current <- overall[[curchr]]
		if (! (curchr %in% names(current)) ) { next }
	
		# Loading all variables, and choosing which ones to throw out.
		curfile <- file.path(dir, current[[curchr]])
		mytab <- read.table(curfile, header=TRUE)
		flagged <- .getFlag(mytab$anchor.pos, mytab$target.pos)
		gapped <- .getGap(mytab$anchor.pos, mytab$target.pos)
		discard <- gapped < min.gap & flagged==2L

		if (all(discard)) { 
			unlink(curfile)
			reindex <- append(reindex, curchr)
		} else if (any(discard)) {
			mytab.alt <- mytab[!discard,,drop=FALSE]
			.saveExt(x=mytab.alt, fname=curfile)

			# Saving it into a discard file, in a chr:start:end format.
			discarded <- discarded + sum(discard)
			if (!is.null(discard.to)) {
				write.table(file = discard.to, data.frame(curchr, -mytab$target.pos[discard], mytab$anchor.pos[discard]),
					col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE, append=TRUE)
			}
		}
	}
	
	# Deleting the absent entries from the index file.
	if (length(reindex)) { 
		curdex <- .getIndex(dir)
		collected <- read.table(curdex,stringsAsFactors=FALSE)
		discard <- collected[,1]==collected[,2] & collected[,1] %in% reindex
		collected <- collected[!discard,,drop=FALSE]
        write.table(file=curdex, collected, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
	}	
	return(invisible(discarded))
}
