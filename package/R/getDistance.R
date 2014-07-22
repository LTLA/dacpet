getDistance <- function(data) 
# Gets the gap distance between all interacting pairs in a data object.
# NA values indicate interchromosomals, whereas negative values indicate
# an interaction between overlapping regions.
{
	output <- integer(nrow(data$pairs))
	ax <- getAnchor(data)
	tx <- getTarget(data)
	is.intra <- as.logical(seqnames(ax)==seqnames(tx))
	output[!is.intra] <- NA

	ax <- ax[is.intra]
	tx <- tx[is.intra]
	output[is.intra] <- pmax(start(ax), start(tx)) - pmin(end(ax), end(tx)) - 1L
	return(output)
}

getAnchor <- function(data) 
# Gets the requested anchor regions.
{
	data$region[data$pairs$anchor]
}

getTarget <- function(data) 
# Gets the requested target regions.
{
	data$region[data$pairs$target]
}
