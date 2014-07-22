stripOutwardPET <- function(file.in, file.out=file.in, min.gap=1e4, discard.to=NULL)
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
	file.tmp <- tempfile(tmpdir=".", fileext=".h5")
	on.exit({if (file.exists(file.tmp)) { unlink(file.tmp) }})

	.initializeH5(file.tmp, .getLengths(file.in))
	discarded <- 0L 
	if (!is.null(discard.to)) {
		dhandle <- file(discard.to, open="w")
		on.exit(close(dhandle), add=TRUE)
	}
	
	# Running through each chromosome.
    overall<-.loadIndices(file.in)
	reindex <- NULL
    for (anchor in names(overall)) {
		current <- overall[[anchor]]
		launcher <- FALSE
		for (target in names(current)) { 
			reads <- .getPairs(file.in, anchor, target)

			if (anchor==target) {
				# Loading all variables, and choosing which ones to throw out.
				flagged <- .getFlag(reads$anchor.pos, reads$target.pos)
				gapped <- .getGap(reads$anchor.pos, reads$target.pos)
				discard <- gapped < min.gap & flagged==2L

				# Saving it into a discard file, in a chr:start:end format.
				discarded <- discarded + sum(discard)
				if (!is.null(discard.to) && any(discard)) { 
					write.table(file = dhandle, data.frame(anchor, -reads$target.pos[discard], reads$anchor.pos[discard]),
						col.names=FALSE, quote=FALSE, sep="\t", row.names=FALSE)
				}
				reads <- reads[!discard,]
			}
			
			if (!nrow(reads)) { next }
			if (!launcher) { 
				.addGroup(file.tmp, anchor)
				launcher <- TRUE
			}
			.writePairs(file.tmp, anchor, target, reads)
		}
	}

	# Moving it to the destination.
	if (!file.rename(file.tmp, file.out)) { stop("failed to move temporary file to the specified destination") }
	return(invisible(discarded))
}
