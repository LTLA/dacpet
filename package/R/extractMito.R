extractMito <- function(file, self=FALSE, restrict=NULL, name="chrM")
# Pretty much does what it says. Extracts all counts for interactions
# involving the mitochondrial genome. By default, it doesn't extract
# counts involving interactions between the mitochondrial genome and 
# itself. It can also be regulated further by specifying a non-trivial
# restrict argument.
# 
# written by Aaron Lun
# 21 January, 2014
{
	# Setting up the output list.
	collected <- list()
	all.len <- h5read(file, 'lengths')
	for (i in 1:nrow(all.len)) { 
		cur.chr <- all.len$chr[i]
		if (is.null(restrict) || cur.chr %in% restrict) { collected[[cur.chr]] <- 0L }
	}
	
	indices <- .loadIndices(file)
	for (anchor in names(indices)) { 
		for (target in names(indices[[anchor]])) {
			if (!self && anchor==name && target==name) { next }
			if (anchor!=name && target!=name) { next }
			if (anchor==name) {
				if (target %in% names(collected)) { collected[[target]] <- indices[[anchor]][[target]] }
			} else if (target==name) { 
				if (anchor %in% names(collected)) { collected[[anchor]] <- indices[[anchor]][[target]] }
			}
		}
	}

	return(unlist(collected))
}
