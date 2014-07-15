###################################################################################################
# This script is designed to test the pair-identifying capabilities of dacpet, and whether it
# saves the results in a proper format.

suppressPackageStartupMessages(require("dacpet"))
source("simsam.R")

comp<-function (fname, chromos, npairs, singles=0, rlen=10, spacer=10, 
		yield=max(1L, round(npairs/runif(1, 2, 10)))) {
	rlen<-as.integer(rlen+0.5)
	spacer<-as.integer(spacer+0.5)

	# Randomly generating reads.
    names<-paste('x', rep(1:npairs, 2), sep=".");
    chrs<-sample(length(chromos), length(names), replace=TRUE);
    pos<-integer(length(names));

    # Assigning positions to all of them.
    for (i in 1:length(chromos)) {
        current<-chrs==i;
        pos[current]<-as.integer(round(runif(sum(current), 1, chromos[[i]])))
    }

    # Throwing them into the SAM file generator. 
    str<-rbinom(length(names), 1, 0.5)==1; 
    sstr<-rbinom(singles, 1, 0.5)==1; 
	reversi<-c(1:npairs+npairs, 1:npairs)
    out<-simsam(fname, names(chromos)[chrs], pos, str, chromos, names=names, len=rlen, 
			is.first=c(rep(TRUE, npairs), rep(FALSE, npairs)), is.paired=TRUE,
			mate.chr=names(chromos)[chrs][reversi], mate.pos=pos[reversi], mate.str=str[reversi])

	if (singles) { 
	    # Adding some singles.
    	snames<-schrs<-spos<-NULL
        snames<-paste('y', 1:singles, sep=".");
        schrs<-sample(length(chromos), singles, replace=TRUE);
        spos<-integer(singles);
		for (i in 1:length(chromos)) {
	       	scurrent<-schrs==i;
        	spos[scurrent]<-as.integer(round(runif(sum(scurrent), 1, chromos[i])))
		}

    	tempname<-file.path(dir, "temp")
		out2<-simsam(tempname, names(chromos)[schrs], spos, sstr, chromos, names=snames, len=rlen)
		more.temp<-file.path(dir, "temp2")
		out<-mergeBam(c(out, out2), more.temp, indexDestination=TRUE, overwrite=TRUE)
		file.rename(more.temp, out)
	}

	# Resorting by name.
	temp<-sortBam(out, "temp", byQname=TRUE)
	file.rename(temp, out)

	################ THEORETICAL MATCH ###################

	# Actually assembling the theoretical values.
	primary<-1:npairs
	pchrs<-chrs[primary]
	pstr<-str[primary]
	ppos<-pos[primary]
	ppos[pstr]<-ppos[pstr]+rlen-1L
	
	secondary<-npairs+primary
	schrs<-chrs[secondary]
	sstr<-str[secondary]
	spos<-pos[secondary]
	spos[sstr]<-spos[sstr]+rlen-1L

	modppos<-ifelse(pstr, ppos, -ppos)
	modspos<-ifelse(sstr, spos, -spos)

	################ ACTUAL MATCH ###################
	# Assembling the output list for comparison.

	tmpdir<-paste0(fname, "_temp")
	preparePET(out, dir=tmpdir, yield=yield)
	indices<-dacpet:::.loadIndices(tmpdir)
	used<-indices

	for (i in 1:length(chromos)) {
		for (j in 1:i) {
			stuff<-(chrs[primary]==i & chrs[secondary]==j) | (chrs[primary]==j & chrs[secondary]==i) 
			if (i > j) {
				is.anchor<-pchrs[stuff]==i 
			} else {
				is.anchor<-ppos[stuff] > spos[stuff]
			}
			cur.anchor<-ifelse(is.anchor, modppos[stuff], modspos[stuff])
			cur.target<-ifelse(is.anchor, modspos[stuff], modppos[stuff])
			o<-order(cur.anchor, cur.target)

			# Checking anchor/target/length/orientation/gap statistics (sorting to ensure comparability).
			achr<-names(chromos)[i]
			tchr<-names(chromos)[j]
            used[[achr]][[tchr]]<-NULL
			if (!(achr%in%names(indices)) || !(tchr %in% names(indices[[achr]]))) { 
				if (length(o)) { stop("true interactions are missing"); }
				next; 
			}
			current<-read.table(file.path(tmpdir, indices[[achr]][[tchr]]), header=TRUE, colClasses="integer")
			current<-current[order(current$anchor.pos, current$target.pos),]
			if (!identical(current$anchor.pos, cur.anchor[o])) { stop("mismatch in anchor elements") }
			if (!identical(current$target.pos, cur.target[o])) { stop("mismatch in target elements") }
		}
	}

	# Checking there's nothing left.
	if (!is.null(unlist(used))) { stop("files left unused in the directory") }
	return(head(read.table(file.path(tmpdir, indices[[1]][[1]]), header=TRUE)))
}

####################################################################################################
# Initiating testing with something fairly benign.

set.seed(348753485)
dir<-"prep-test"
dir.create(dir)
fname<-file.path(dir, "out");
chromos<-c(chrA=20000L, chrB=1000L, chrC=50000L);

comp(fname, npairs=20, chromos=chromos);
comp(fname, npairs=50, chromos=chromos);
comp(fname, npairs=100, chromos=chromos);

# Ramping up the aggression in terms of overlap density

comp(fname, npairs=20, chromos=c(chrA=300));
comp(fname, npairs=50, chromos=c(chrA=300));
comp(fname, npairs=100, chromos=c(chrA=300));

# Increasing the number of reads all round.

comp(fname, npairs=200, chromos=chromos);
comp(fname, npairs=1000, chromos=chromos);

# Adding some singletons.

comp(fname, npairs=200, chromos=chromos, singles=10)
comp(fname, npairs=200, chromos=chromos, singles=50)

# Making the fragments smaller.

comp(fname, npairs=200, rlen=50, chromos=chromos)
comp(fname, npairs=500, rlen=50, chromos=chromos)

# Trying it out with some more elements in a more restricted space.

comp(fname, npairs=500, chromos=c(chrA=20L))
comp(fname, npairs=200, chromos=c(chrA=20L, chrB=10L))
comp(fname, npairs=1000, chromos=c(chrA=10L))

# Adding lots of chromosomes.

chromos<-c(chrA=500L, chrB=6000L, chrC=70L, chrD=400L, chrE=3000L, chrF=10L, chrG=200L)
comp(fname, npairs=200, chromos=chromos);
comp(fname, npairs=1000, chromos=chromos);
comp(fname, npairs=200, chromos=chromos);
comp(fname, npairs=1000, chromos=chromos);

###################################################################################################

unlink(dir, recursive=TRUE)
