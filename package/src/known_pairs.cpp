#include "dacpet.h"

/* This function computes the counts for each pairwise combination
 * of known regions, given the overlaps for each read with said regions.
 * The indices are derived from findOverlaps, so we assume they're sorted
 * with respect to the queries (i.e. aq and tq).
 */

typedef std::pair<int, int> coord;

SEXP count_all_pairs (SEXP aq, SEXP as, SEXP tq, SEXP ts, SEXP samechr) try {
	if (!isInteger(aq)) { throw std::runtime_error("anchor query index must be integer"); }
	if (!isInteger(as)) { throw std::runtime_error("anchor subject index must be integer"); }
	if (!isInteger(tq)) { throw std::runtime_error("target query index must be integer"); }
	if (!isInteger(ts)) { throw std::runtime_error("target subject index must be integer"); }

	const int na=LENGTH(aq), nt=LENGTH(tq);
	if (na!=LENGTH(as) || nt!=LENGTH(ts)) { throw std::runtime_error("subject and query vectors must have equal length"); }
	const int* aqptr=INTEGER(aq),
		  *asptr=INTEGER(as),
		  *tqptr=INTEGER(tq),
		  *tsptr=INTEGER(ts);

	/* If we're on the same chromosome, some extra protection is
	 * necessary to avoid counting some things twice, if both 
	 * PETs overlap a subset of the same features.
	 */
	if (!isLogical(samechr) || LENGTH(samechr)!=1) { throw std::runtime_error("same chromosome specifier must be a logical scalar"); }
	const bool issame=asLogical(samechr);
	std::set<coord> sorter;

	// Setting up some variables.
	std::map<coord, int> collected;
	int i=0, j=0, jlast, jcopy;
	
	// Running through and collecting all region pairs for each PET.
	while (i < na && j<nt) {
		const int& currenta=aqptr[i];
		const int& currentt=tqptr[j];
		if (currenta > currentt) {
			++j;
			continue;
		} else if (currenta < currentt) {
			++i;
			continue;
		}
		
		// Searching for the last 'j'.
		jlast=j+1;
		while (jlast < nt && currentt==tqptr[jlast]) { ++jlast; }

		// Running through and computing all combinations.
		while (i < na && currenta==aqptr[i]) {
			for (jcopy=j; jcopy < jlast; ++jcopy) {
				const int& anchor=asptr[i];
				const int& target=tsptr[jcopy];
				if (!issame) { 
					++(collected[coord(anchor, target)]);
				} else {
					if (anchor >= target) {
						sorter.insert(coord(anchor, target));
					} else {
						sorter.insert(coord(target, anchor));
					}
				}
			}
			++i;
		}
		j=jlast; 

		/* This approach ensures we never count one PET twice for any given interaction.
 		 * I don't think we can avoid iterating over them, using a rule like "break on 
 		 * anchor < target". This is because you can get unique features in the middle
 		 * of a common subset. For eample, the target can overlap [1,2,3] and the anchor
 		 * might overlap only [1,3] if "2" starts between "1" and "3' but is shorter. 
 		 * If you use the rule above, you just end up with [1,1], [3,1], [3,2] and [3,3];
 		 * this misses out on [2,1].
 		 */
		if (issame) {
			for (std::set<coord>::const_iterator its=sorter.begin(); its!=sorter.end(); ++its)  { ++(collected[*its]); }
			sorter.clear();
		}
	}

	SEXP output=PROTECT(allocVector(VECSXP, 3));
	try {
		const int nout=collected.size();
		std::deque<int*> optrs;
		for (int x=0; x<3; ++x) {
			SET_VECTOR_ELT(output, x, allocVector(INTSXP, nout));
			optrs.push_back(INTEGER(VECTOR_ELT(output, x)));
		}
		for (std::map<coord, int>::const_iterator itc=collected.begin(); 
				itc!=collected.end(); ++itc) {
			optrs[0][0]=(itc->first).first;
			++optrs[0];
			optrs[1][0]=(itc->first).second;
			++optrs[1];
			optrs[2][0]=(itc->second);
			++optrs[2];
		}
	} catch (std::exception & e){
		UNPROTECT(1);
		throw;
	}

	UNPROTECT(1);
	return output;
} catch (std::exception &e) {
	return mkString(e.what());
}

// A comparator, for use below.

struct comp_pairs {
	comp_pairs(const int& nlib) : aptrs(nlib), tptrs(nlib) {}
	bool operator() (const coord& l, const coord& r) const {
		if (aptrs[l.first][l.second] > aptrs[r.first][r.second]) { 
			return true; 
		} else if (aptrs[l.first][l.second]==aptrs[r.first][r.second]) { 
			return (tptrs[l.first][l.second] > tptrs[r.first][r.second]);
		} else {
			return false;
		}
	}
    std::deque<const int*> aptrs, tptrs;
};

/* This function aggregates results from multiple libraries. 
 * It does some filtering to avoid aggregation of trivial
 * results.
 */

SEXP aggregate_pair_counts (SEXP collected, SEXP filter) try {
	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value should be an integer scalar"); }
	const int filt=asInteger(filter);
	if (!isNewList(collected)) { throw std::runtime_error("input must be a list"); }
	const int nlib=LENGTH(collected);
	comp_pairs ptr_holder(nlib);
	std::deque<const int*>& aptrs=ptr_holder.aptrs;
	std::deque<const int*>& tptrs=ptr_holder.tptrs;
	std::deque<const int*> cptrs(nlib);
	
	// Checking the status of each library.
	std::deque<int> nrows(nlib);
	for (int i=0; i<nlib; ++i) {
		SEXP current=VECTOR_ELT(collected, i);
		if (!isNewList(current) || LENGTH(current)!=3) { throw std::runtime_error("each library must have a list of 3 vectors"); }

		for (int j=0; j<3; ++j) {
			SEXP curvec=VECTOR_ELT(current, j);
			if (!isInteger(curvec)) { throw std::runtime_error("vectors should be integer"); }
			if (j==0) { 
				nrows[i]=LENGTH(curvec);
				aptrs[i]=INTEGER(curvec);
			} else {
				if (nrows[i]!=LENGTH(curvec)) { throw std::runtime_error("vector lengths should be consistent within a library"); }
				(j==1 ? tptrs[i] : cptrs[i])=INTEGER(curvec);
			}
		}
	}
	
	// Filling the priority queue.
	std::priority_queue<coord, std::deque<coord>, comp_pairs> sofar(ptr_holder);
	for (int i=0; i<nlib; ++i) {
		if (nrows[i]) { sofar.push(coord(i, 0)); }
	}

	/* Aggregating them together. We assume that the anchor/target pairs have
	 * been sorted by anchor, then target. This should be true if they're coming 
	 * fresh out of the previous function (due to enforced sorting in <map>).
	 */
	std::deque<int> anchor, target, counts;
	std::deque<int> currentcounts(nlib);
	int lib, dex;
	while (sofar.size()) {
		lib=sofar.top().first;
		dex=sofar.top().second;
		const int& currenta=aptrs[lib][dex];
		const int& currentt=tptrs[lib][dex];
		do {
			currentcounts[lib]+=cptrs[lib][dex];
			sofar.pop();
			if (++dex != nrows[lib]) {
				sofar.push(coord(lib, dex));
			}
			if (sofar.empty()) { 
				break; 
			} else {
				lib=sofar.top().first;
				dex=sofar.top().second;
			}
		} while (aptrs[lib][dex]==currenta && tptrs[lib][dex]==currentt);
		
		// Checking if we're above the filter.
	 	int totes=0;
		for (int i=0; i<nlib; ++i){ totes+=currentcounts[i]; }
		if (totes >= filt) { 
			counts.insert(counts.end(), currentcounts.begin(), currentcounts.end());
			anchor.push_back(currenta);
			target.push_back(currentt);
		}
		std::fill(currentcounts.begin(), currentcounts.end(), 0);
	}

	// Storing results.
	SEXP output=PROTECT(allocVector(VECSXP, 3));
	try {
		const int nrow_out=anchor.size();
		SET_VECTOR_ELT(output, 0, allocVector(INTSXP, nrow_out));
		std::copy(anchor.begin(), anchor.end(), INTEGER(VECTOR_ELT(output, 0)));
		SET_VECTOR_ELT(output, 1, allocVector(INTSXP, nrow_out));
		std::copy(target.begin(), target.end(), INTEGER(VECTOR_ELT(output, 1)));

		SET_VECTOR_ELT(output, 2, allocMatrix(INTSXP, nrow_out, nlib));
		std::deque<int*> optrs(nlib);
		optrs[0]=INTEGER(VECTOR_ELT(output, 2));
		for (int i=1; i<nlib; ++i) { optrs[i]=optrs[i-1]+nrow_out; }
		
		std::deque<int>::const_iterator itc=counts.begin();
		for (int i=0; i<nrow_out; ++i) {
			for (int j=0; j<nlib; ++j) {
				optrs[j][i]=*itc;
				++itc;
			}			
		}
	} catch (std::exception& e) {
		UNPROTECT(1);
		throw;
	}
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}
