#include "dacpet.h"

/* This function computes the counts for each region. Some subtlety 
 * is needed whenever the anchor and target reads overlap with the 
 * same genomic region, in which case that pair should only be counted once.
 */

typedef std::pair<int, bool> region;

SEXP count_margins (SEXP aq, SEXP as, SEXP tq, SEXP ts, SEXP samechr) try {
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
	std::set<region> sorter;

	// Setting up some variables.
	std::map<region, int> collected;
	int i=0, j=0, current;
	
	// Running through and collecting all region pairs for each PET.
	while (1) {
		if (i<na) {
 		   	if (j<nt) {
				current=(aqptr[i] <= tqptr[j] ? aqptr[i] : tqptr[j]);
			} else {
				current=aqptr[i];
			}
		} else {
			if (j < nt) {
				current=tqptr[j];				
			} else {
				break;
			}
		}

		// Adding anchors (and advancing, if current==currenta).
		while (i < na && current==aqptr[i]) {
			if (!issame) {
				++(collected[region(asptr[i], true)]);
			} else {
				sorter.insert(region(asptr[i], true));
			}
			++i;
		}

		// Adding targets (and advancing, if current==currentt).
		while (j < nt && current==tqptr[j]) {
			if (!issame) { 
				++(collected[region(tsptr[j], false)]);				
			} else {
				sorter.insert(region(tsptr[j], true)); // Deliberate, so only unique ones are kept.
			}
			++j;
		}

		// This approach ensures we never count one PET twice for any given interaction.
		if (issame) {
			for (std::set<region>::const_iterator its=sorter.begin(); its!=sorter.end(); ++its)  { ++(collected[*its]); }
			sorter.clear();
		}
	}

	SEXP output=PROTECT(allocVector(VECSXP, 3));
	try {
		const int nout=collected.size();
		std::deque<int*> optrs;
		SET_VECTOR_ELT(output, 0, allocVector(INTSXP, nout));
		optrs.push_back(INTEGER(VECTOR_ELT(output, 0)));
		SET_VECTOR_ELT(output, 1, allocVector(LGLSXP, nout));
		optrs.push_back(LOGICAL(VECTOR_ELT(output, 1)));
		SET_VECTOR_ELT(output, 2, allocVector(INTSXP, nout));
		optrs.push_back(INTEGER(VECTOR_ELT(output, 2)));

		for (std::map<region, int>::const_iterator itc=collected.begin(); itc!=collected.end(); ++itc) {
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


