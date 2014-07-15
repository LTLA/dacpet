#ifndef COUNT_MANAGER_H
#define COUNT_MANAGER_H
#include "dacpet.h"

/* 
 * A class to manage and classify unique counts. This is necessary to store count combinations 
 * w/o the excessive overhead associated with individual deques.
 */

class count_manager {
public:
    count_manager(const int&);
	int& back(const int&);
	const int& access(const int&, const int&) const;
	const int& advance();
	SEXP get_counts() const;
	int get_rows() const;
	const int nlibs;
private:
	// Assorted structures in use.
	int current;
	std::deque<std::deque<int> > count_store;

	// Internal count typedef.
	struct comp {
		comp(const count_manager*);
		const count_manager* parent_ptr;
		bool operator() (const int&, const int&) const;
	};
   	std::set<int, comp> unique_combos;
};

extern "C" {

SEXP create_counts(SEXP);

SEXP get_counts(SEXP);

}

#endif
