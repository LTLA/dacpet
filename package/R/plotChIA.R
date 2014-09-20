plotChIA <- function(file, anchor, target=anchor, 
   width=1000, cap=50, xlab=NULL, ylab=NULL, col="red", diag=TRUE, ...)
# This generates a ChIA-PET plot, akin to the plaid plots for Hi-C data. 
#
# written by Aaron Lun
# 2 December 2013	
{
	target.lim<-c(start(target), end(target))
	target<-as.character(seqnames(target))
	anchor.lim<-c(start(anchor), end(anchor))
	anchor<-as.character(seqnames(anchor))
   	if (length(anchor)!=1L || length(anchor.lim)!=2L) { stop("exactly one anchor range is required for plotting") }
    if (length(target)!=1L || length(target.lim)!=2L) { stop("exactly one target range is required for plotting") }

	# Choosing which file to open (accounting for mixed specification of anchor/target, relative to the definition in the file).
	indexed <- .loadIndices(file)
    if (!is.null(indexed[[anchor]][[target]])) {
		allpts <- .getPairs(file, anchor, target)
	} else if (!is.null(indexed[[target]][[anchor]])) {
		allpts <- .getPairs(file, target, anchor)
        allpts[,1:2] <- allpts[,2:1] 
	} else {
		allpts <- data.frame(anchor.pos=integer(0), target.pos=integer(0))
	}

	# Setting up the plot.
	if (is.null(xlab)) { xlab <- anchor }
	if (is.null(ylab)) { ylab <- target }
	plot(0,0,type="n",xlim=anchor.lim, ylim=target.lim, xlab=xlab, ylab=ylab, ...) 

	# Loading points and lengths, picking which reads to keep (auto-extending here to get things on the plot boundaries).
	# It'll also plot things past the diagonal in intra-chromosomal plots, if any are requested.
	apos <- abs(allpts$anchor.pos)
	tpos <- abs(allpts$target.pos)
	a.ext <- 0.1*(anchor.lim[2] - anchor.lim[1])
	t.ext <- 0.1*(target.lim[2] - target.lim[1])
	end.a <- anchor.lim[2]+a.ext 
	start.a <- anchor.lim[1]-a.ext 
	end.t <- target.lim[2]+t.ext
	start.t <- target.lim[1]-t.ext

	decol <- col2rgb(col)
	for (it in 1:2) {
		# Figuring out which bin they belong to.
		kept <- apos <= end.a & apos >= start.a & tpos >= start.t & tpos <= end.t
	    aix <- as.integer(apos[kept]/width)
		tix <- as.integer(tpos[kept]/width)
		o <- order(aix, tix)
		indexed <- c(which(diff(aix[o])!=0L | diff(tix[o])!=0L), length(o))
	    counts <- c(1L, diff(indexed))
	    
		coord.a <- aix[o][indexed]
		coord.t <- tix[o][indexed]
		if (it==2L) { 
			# Get rid of the guys on the diagonal, don't want to plot them twice.
			retain <- coord.a!=coord.t
			coord.a <- coord.a[retain]
			coord.t <- coord.t[retain]
			counts <- counts[retain]
		}
		
	    # You'll need to define start and end separately for anchor and target.
		col <- rgb(decol[1,], decol[2,], decol[3,], alpha = 255*pmin(1, counts/cap), maxColorValue=255)
		rect(coord.a*width, coord.t*width, (coord.a+1)*width, (coord.t+1)*width, border=NA, col=col)

		# Deciding whether to repeat the plotting with unplotted elements past the diagonal.
		if (anchor==target && diag) { 
			temp <- apos
			apos <- tpos
			tpos <- temp
		} else { break }
	}
	
	return(invisible(NULL))
}



