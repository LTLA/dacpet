splitLinkers <- function(fastq1, fastq2, linkerA, linkerB, 
    min.score=10, start=NA, prefix="out", read.prefix=NULL)
# This function takes a couple of FASTQ files and pulls the linker
# out from them. It then stores the tags in files starting with
# 'prefix'. Linker identification is performed with SW alignment,
# starting from 'start' and with a minimum score.
{
	min.score <- as.integer(min.score)
	start <- as.integer(start)
	read.prefix <- as.character(read.prefix)

	# Splitting, as requested.
	returned <- .Call(cxx_split_linkers, fastq1, fastq2, linkerA, linkerB, 
		prefix, min.score, start, read.prefix)
	if (is.character(returned)) { stop(returned) }

	# Returning file names.
	all.files <- NULL
	for (x in 1:2) { 
		all.files <- append(all.files, 
			paste0(prefix, "_", c("AA", "AB", "BB", "other"), "_", x, ".fastq"))
	}
	return(all.files)
}
