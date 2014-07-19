##########################################################
# This script tests the various functions in the linker splitting machinery.

suppressPackageStartupMessages(require(dacpet))
suppressPackageStartupMessages(require(Biostrings))

compfun <- function(alpha, bravo, m=1, mm=-1, go=-3, ge=-1) {
	out <- .Call(dacpet:::cxx_test_align, alpha, bravo, as.integer(m),
		as.integer(mm), as.integer(go), as.integer(ge))
	if (is.character(out)) { stop(out) }

	mat <- nucleotideSubstitutionMatrix(match = m, mismatch = mm, baseOnly = TRUE)
	ref <- pairwiseAlignment(toupper(alpha), toupper(bravo), type="local", substitutionMatrix=mat, gapOpening=go-ge, gapExtension=ge, scoreOnly=TRUE) 
	# A Gotcha, here; gapOpening needs to subtract out gapExtension.

	if (ref!=out) { 
		stop(sprintf("mismatch in alignment scores: %i vs %i", out, ref)) 
	}
	return(ref)
}

original <- modified <- "GTTGGAATGTATATCG"
compfun(original, modified) # Perfect match
compfun(original, modified, mm=-3) 
compfun(original, modified, m=2) 
compfun(original, modified, go=-2) 
compfun(original, modified, ge=-2) 

modified <- "GTTGGAATGTAaATCG" # Mismatch at back
compfun(original, modified) 
compfun(original, modified, mm=-3) 
compfun(original, modified, m=2) 
compfun(original, modified, go=-2) 
compfun(original, modified, ge=-2) 

modified <- "GTcGGAATGTATATCG" # Mismatch at front
compfun(original, modified) 
compfun(original, modified, mm=-3) 
compfun(original, modified, m=2) 
compfun(original, modified, go=-2) 
compfun(original, modified, ge=-2) 

modified <- "GTTGGATAAGATATCG" # Wrong linker
compfun(original, modified) 
compfun(original, modified, mm=-3) 
compfun(original, modified, m=2) 
compfun(original, modified, go=-2) 
compfun(original, modified, ge=-2) 

modified <- "GTTGAATGTATATCG" # Deletion in the front
compfun(original, modified) 
compfun(original, modified, mm=-3) 
compfun(original, modified, m=2) 
compfun(original, modified, go=-2) 
compfun(original, modified, ge=-2) 

modified <- "GTTGGAAATATCG" # Deletion in the middle
compfun(original, modified) 
compfun(original, modified, mm=-3) 
compfun(original, modified, m=2) 
compfun(original, modified, go=-2) 
compfun(original, modified, ge=-2) 

modified <- "GTTGGAATGTATATG" # Deletion in the back
compfun(original, modified) 
compfun(original, modified, mm=-3) 
compfun(original, modified, m=2) 
compfun(original, modified, go=-2) 
compfun(original, modified, ge=-2) 

modified <- "GTTGcGAATGTATATCG" # Insertion in front
compfun(original, modified) 
compfun(original, modified, mm=-3) 
compfun(original, modified, m=2) 
compfun(original, modified, go=-2) 
compfun(original, modified, ge=-2) 

modified <- "GTTGGAATGagTATATCG" # Insertion in middle
compfun(original, modified) 
compfun(original, modified, mm=-3) 
compfun(original, modified, m=2) 
compfun(original, modified, go=-2) 
compfun(original, modified, ge=-2) 

modified <- "GTTGGAATGTATATcCG" # Insertion in middle
compfun(original, modified) 
compfun(original, modified, mm=-3) 
compfun(original, modified, m=2) 
compfun(original, modified, go=-2) 
compfun(original, modified, ge=-2) 

try(
compfun(original, "GTTCGAATCTATACCG") # No matching subwords.
)
try(
compfun(original, "CAGGGACTAGCATGTCAT") # No matching subwords.
)

##########################################################
# Now, testing the actual linker splitting itself. We make some tags of varying quality
# and we throw them in, and see whether they show up in the right file.

outdir <- "link-test"
if (file.exists(outdir)) { unlink(outdir, recursive=TRUE) }
dir.create(outdir)

genqual <- function(n) {
	qual <- sample(c("H", "I", "#", "^", "-", "A", "B", "C", 0:9), n, replace=TRUE)
	paste(qual, collapse="")
}

fastq1 <- file.path(outdir, "out_1.fastq")
fastq2 <- file.path(outdir, "out_2.fastq")
generator <- function(seq1, seq2, readname) {
	write(c(paste0("@", readname), seq1, "+", genqual(nchar(seq1))), 
		file=fastq1, append=TRUE)
	write(c(paste0("@", readname), seq2, "+", genqual(nchar(seq2))),
		file=fastq2, append=TRUE)
}

linkA <- "GTTGGAATGTATATCG"
linkB <- "GTTGGATAAGATATCG"
runNcomp <- function(start, prefix=file.path(outdir, "gunkitorious"), rename=NULL) {
	out <- splitLinkers(fastq1, fastq2, linkerA=linkA, linkerB=linkB, prefix=prefix, start=start, read.prefix=rename)

	# Pulling reference lines.
	lines1 <- readLines(fastq1)
	is.core <- 0:(length(lines1)/4 - 1L)*4L
	is.seq <- is.core + 2L
	is.qual <- is.core + 4L
	ender <- ifelse(is.na(start), nchar(lines1[is.seq][1]) - nchar(linkA), start-1L)
	lines1[is.seq] <- substr(lines1[is.seq], 1L, ender)
	lines1[is.qual] <- substr(lines1[is.qual], 1L, ender)

	lines2 <- readLines(fastq2)
	lines2[is.seq] <- substr(lines2[is.seq], 1L, ender)
	lines2[is.qual] <- substr(lines2[is.qual], 1L, ender)

	in.aa <- t(outer(grep("^@AA", lines1), 0:3, FUN="+"))
	in.ab <- t(outer(grep("^@(AB|BA)", lines1), 0:3, FUN="+"))
	in.bb <- t(outer(grep("^@BB", lines1), 0:3, FUN="+"))
	in.other <- t(outer(grep("^@other", lines1), 0:3, FUN="+"))

	# Renaming, if necessary.
	if (!is.null(rename)) { 
		is.name <- is.core + 1L
		lines1[is.name] <- lines2[is.name] <- paste0("@", rename, ".", 1:length(is.core))
	}

	# Comparing.
	if (!identical(lines1[in.aa], readLines(out[grep("AA_1", out)]))) { stop("mismatch for AA, read 1") }
	if (!identical(lines1[in.ab], readLines(out[grep("AB_1", out)]))) { stop("mismatch for AB, read 1") }
	if (!identical(lines1[in.bb], readLines(out[grep("BB_1", out)]))) { stop("mismatch for AA, read 1") }
	if (!identical(lines1[in.other], readLines(out[grep("other_1", out)]))) { stop("mismatch for AA, read 1") }

	if (!identical(lines2[in.aa], readLines(out[grep("AA_2", out)]))) { stop("mismatch for AA, read 2") }
	if (!identical(lines2[in.ab], readLines(out[grep("AB_2", out)]))) { stop("mismatch for AB, read 2") }
	if (!identical(lines2[in.bb], readLines(out[grep("BB_2", out)]))) { stop("mismatch for AA, read 2") }
	if (!identical(lines2[in.other], readLines(out[grep("other_2", out)]))) { stop("mismatch for AA, read 2") }

	return(out)
}

### Running through an example.
generator(paste0("XUIXIArUXIZOCIEDH", linkA), paste0("ASaVDIOSDFGGDGOGS", linkA), "AA.one") 
generator(paste0("XasdadXgvdgryufD1", linkA), paste0("ASVDI12s3ifosjvkc", linkA), "AA.two") 
generator(paste0("asdasdavhoiujbod3", linkA), paste0("ALSKMJsas;ldkasln", linkB), "AB.one") 
generator(paste0("asdasdavouiujbodd", linkA), paste0("ALSaKxc;vkjxcvvBG", linkB), "AB.two") 
generator(paste0("aiiajdoajsjdbfbDe", linkB), paste0("kSqodwiusovcnvq9S", linkA), "BA.one") 
generator(paste0("aiiajdoajsdgbfbdf", linkB), paste0("dSqowiusovcnvq9rS", linkA), "BA.two") 
generator(paste0("XUdfsdfsjdfpssdpa", linkB), paste0("VUDHIghFVSHDISras", linkB), "BB.one") 
generator(paste0("XUdfsdfsjdfpssadj", linkB), paste0("VUDHIFVSHDIsdSFas", linkB), "BB.two") 

generator(paste0("asdfafdvnjghCIEbH", "GTTGGAACCCATATCG"), paste0("Asdffsidflenmkfbf", linkA), "AA.x.one") # (score == 10)
generator(paste0("asdadsdvnjghCdEDH", linkB), paste0("sAsdfsidflenmkfbf", "GTTGGATCCCATATCG"), "BB.x.one") 
generator(paste0("asassssvnjghCIEDu", "GTTGGACCCCATATCG"), paste0("Asdfbsidflenmkfbf", linkB), "other.one") # Fails with an extra mismatch
generator(paste0("XUIsdlkvjslvmyfdH", linkA), paste0("ASad;oasjdsnapodk", "GTTGGAACCCCTATCG"), "other.two") 
generator(paste0("PSojfiobjmmghmggh", "GACTGTTATGCTATAC"), paste0("dfflbimdflbmlvkGS", linkA), "other.three") # Something completely wrong
generator(paste0("XsdpspvsdvospnovH", linkB), paste0("S;poqkwepfdfokoGS", "GTCACGATGGAGCGGG"), "other.four") 

runNcomp(start=NA) 
runNcomp(start=18, prefix=file.path(outdir, "blah")) # position 18 is where the linker starts (tag is 17 bp long here).
runNcomp(start=NA, prefix=file.path(outdir, "xcxcxa0_da"), rename="phosphorous") 

###########################################################

unlink(outdir, recursive=TRUE)

###########################################################
