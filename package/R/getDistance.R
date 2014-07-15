getDistance <- function(data) 
# Gets the gap distance between all interacting pairs in a data object.
# NA values indicate interchromosomals, whereas negative values indicate
# an interaction between overlapping regions.
{
	output <- integer(nrow(data$counts))
	is.intra <- as.logical(seqnames(data$region[data$pairs$anchor])==seqnames(data$region[data$pairs$target]))
	output[!is.intra] <- NA

	as <- start(data$region)[data$pairs$anchor[is.intra]]
	ts <- start(data$region)[data$pairs$target[is.intra]]
	ae <- end(data$region)[data$pairs$anchor[is.intra]]
	te <- end(data$region)[data$pairs$target[is.intra]]
	output[is.intra] <- pmax(as, ts) - pmin(ae, te) - 1L
	return(output)
}
