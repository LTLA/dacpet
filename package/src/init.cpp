#include "dacpet.h"
#include "count_manager.h"
#include "R_ext/Rdynload.h"
#include "R_ext/Visibility.h"

extern "C" { 

static const R_CallMethodDef all_call_entries[] = {
	{"cluster_2d", (DL_FUNC) &cluster_2d, 5},
	{"count_chia", (DL_FUNC) &count_chia, 4},
//	{"decompress_pairs", (DL_FUNC) &decompress_pairs, 9},
	{"count_all_pairs", (DL_FUNC) &count_all_pairs, 5},
	{"aggregate_pair_counts", (DL_FUNC) &aggregate_pair_counts, 2},
	{"split_linkers", (DL_FUNC) &split_linkers, 8},

	{"create_counts", (DL_FUNC) &create_counts, 1},
	{"get_counts", (DL_FUNC) &get_counts, 1},
  	{NULL, NULL, 0}
};

void attribute_visible R_init_dacpet(DllInfo *info)
{
	R_registerRoutines(info, NULL, all_call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
}

}
