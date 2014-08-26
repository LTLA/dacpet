#ifndef DACPET_H
#define DACPET_H

#include <deque>
#include <map>
#include <set>
#include <queue>
#include <stdexcept>
#include <algorithm>

#include "R.h"
#include "Rinternals.h"

extern "C" {

SEXP cluster_2d (SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP count_chia(SEXP, SEXP, SEXP, SEXP);

//SEXP decompress_pairs (SEXP, SEXP, SEXP, SEXP, SEXP,
//		        SEXP, SEXP, SEXP, SEXP);

SEXP count_all_pairs (SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP aggregate_pair_counts (SEXP, SEXP);

SEXP count_margins (SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP split_linkers(SEXP, SEXP, SEXP, SEXP, 
				SEXP, SEXP, SEXP, SEXP);

SEXP test_align(SEXP, SEXP, 
			SEXP, SEXP, SEXP, SEXP);

}

#endif
