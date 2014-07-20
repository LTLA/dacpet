#include "chia_counter.h"

SEXP count_chia(SEXP extptr, SEXP pulled, SEXP filter, SEXP issame) try {
   	// Launching pointers.
   	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value should be an integer scalar"); }
   	if (!isLogical(issame) || LENGTH(issame)!=1) { throw std::runtime_error("same chromosome specification should be a logical scalar"); }
	const int filt=asInteger(filter);
	const bool is_same_chr=asLogical(issame);

	// Loading the count combinations.
	count_manager* xptr=(count_manager*)R_ExternalPtrAddr(extptr);
	if (xptr==NULL) { throw std::runtime_error("external pointer cannot point to an invalid address"); }

	// Running through the input lists to get the counts.	
	chia_counter cc(pulled, is_same_chr, filt, xptr);
	cc.process();

	SEXP output=PROTECT(allocVector(VECSXP, 3));
	try {
		const int totaln=cc.access_index().size();
		for (int i=0; i<3; ++i) {
			SET_VECTOR_ELT(output, i, allocVector(INTSXP, totaln));
			int* optr=INTEGER(VECTOR_ELT(output, i));
			const std::deque<int> * curcopy=NULL;
			switch(i) {
				case 0:
					curcopy=&(cc.access_anchor());
					break;
				case 1:
					curcopy=&(cc.access_target());
					break;
				case 2:
					curcopy=&(cc.access_index());
					break;
			}
			std::copy(curcopy->begin(), curcopy->end(), optr);
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
