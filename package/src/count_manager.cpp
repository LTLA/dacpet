#include "count_manager.h"

/* Managing the count combinations. */

count_manager::count_manager (const int& n) : nlibs(n), current(0), 
	count_store(n, std::deque<int>(1, 0)), unique_combos(comp(this)) {}

// Accessors to manipulate the counts.

int& count_manager::back(const int& i) { return count_store[i].back(); }

const int& count_manager::access (const int& i, const int& j) const { return count_store[i][j]; }

int count_manager::get_rows() const { return unique_combos.size(); }

/* Checking if the current count is in our list already. If it is, we don't flag it as unique
 * (and we actually are obliged to reset the current thingo to zero, as we're not copying
 * during re-initialization). If it isn't in our list, then we add it. 
 */

const int& count_manager::advance() {
	std::set<int, comp>::iterator ituc=unique_combos.lower_bound(current);
	if (ituc==unique_combos.end() || unique_combos.key_comp()(current, *ituc)) {
		ituc=unique_combos.insert(ituc, current);
		for (int m=0; m<nlibs; ++m) { count_store[m].push_back(0); }
		++current;
	} else {
		for (int m=0; m<nlibs; ++m) { count_store[m].back()=0; }
	}	
	return *ituc;
}

/* Some trickery to compare 'counts' i.e. count indices. */

count_manager::comp::comp(const count_manager* t) : parent_ptr(t) {}

bool count_manager::comp::operator() (const int& lhs, const int& rhs) const {
	const std::deque<std::deque<int> >& current=parent_ptr->count_store;
	const int& n=parent_ptr->nlibs;
    for (int i=0; i<n; ++i) {
        const int& mine=current[i][lhs];
        const int& yours=current[i][rhs];
        if (mine < yours) { return true; }
        else if (mine!=yours) { return false; }
    }
    return false;
}

/* Function to actually return the counts. Note that this must be immediately assigned to something
 * that it protected whenever it is used. We are also assuming 32-bit integers in order to hold the count
 * combinations, though that should be safe enough as R implementations specify INTSXP is 32-bit.
 */

SEXP count_manager::get_counts () const {
	const int& nrow=get_rows();
	SEXP output=PROTECT(allocMatrix(INTSXP, nrow, nlibs));
	try {
		std::deque<int*> optrs(1, INTEGER(output));
		for (int i=1; i<nlibs; ++i) { optrs.push_back(optrs[i-1]+nrow); }
		for (std::set<int, comp>::const_iterator ituc=unique_combos.begin(); ituc!=unique_combos.end(); ++ituc) {
			const int& index=*ituc;
			for (int i=0; i<nlibs; ++i) { optrs[i][index]=count_store[i][index]; }
		}		
	} catch (std::exception& e) {
		UNPROTECT(1);
		throw;
	}
	UNPROTECT(1);
	return output;
}

// Clearing the memory during garbage collection.
extern "C" { 

static void destroy_counts (SEXP ext) {
    count_manager* target=(count_manager*) R_ExternalPtrAddr(ext);
    if (target!=NULL) {
        delete target;
        R_ClearExternalPtr(ext);
    } 
    return;
}

}

// Initializing the pointer to the count_chip object.
SEXP create_counts(SEXP num_libs) try {
    if (!isInteger(num_libs) || LENGTH(num_libs)!=1) { throw std::runtime_error("number of libraries must be an integer scalar"); }
    count_manager* gunk=new count_manager(asInteger(num_libs));
    SEXP ext=PROTECT(R_MakeExternalPtr(gunk, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(ext, destroy_counts, ((Rboolean) TRUE));
    UNPROTECT(1);
    return ext;
} catch (std::exception& e) { 
    return mkString(e.what());
}

SEXP get_counts(SEXP ext) {
    count_manager* target=(count_manager*) R_ExternalPtrAddr(ext);
	if (target==NULL) {
		return mkString("inappropriate external pointer supplied");
	}
	return target->get_counts();	
}

