#ifndef CHIA_COUNTER_H
#define CHIA_COUNTER_H

#include "dacpet.h"
#include "count_manager.h"

enum action { START_ADD, 
	STOP_ADD, 
	START_RM, 
	STOP_RM, 
	EXISTING }; // EXISTING must be last for some loop conditions to hold, w.r.t direct checking at the top of the queue.

struct coord {
	coord(const int& a, const int& t, const int& l, const action& s) :
		anchor(a), target(t), library(l), status(s) { }
	bool operator> (const coord& other) const {
		if (anchor==other.anchor) {
			if (target==other.target) {
				return status>other.status;				
			} 
			return target>other.target;
		} 
		return anchor>other.anchor;
	}
	int anchor, target, library;
	action status;
};

class chia_counter {
public:
	chia_counter(SEXP, bool, int, count_manager*);
	void process();
	const std::deque<int>& access_anchor() const;
	const std::deque<int>& access_target() const;
	const std::deque<int>& access_index() const;
private:
	// Input parameters.
	int nlibs;
	std::deque<const int*> asptrs, tsptrs, aeptrs, teptrs;
	std::deque<int> num, indices;

	// 'next' holds the incoming edits to the 'landscape'.
	std::priority_queue<coord, std::deque<coord>, std::greater<coord> > next;
	bool anchor_is_finished;

	// 'landscape' holds the target:(anchor depth, combination) pairs.
	std::map<int, std::pair<int, int> > landscape;
	std::map<int, std::pair<int, int> >::iterator itl;
	int underlying;
	bool on_existing, allzeros;

	// Choosing how to modify existing combinations.
	std::deque<int> modifier;
	count_manager * cmptr;
	
	// Choosing whether to report.
	const int filter;
	const bool is_same_chr;
	int total;
	bool recorded;
	std::deque<int> all_anchors, all_targets, index;

	// Assorted temporaries.
	int cur_anchor, cur_target, cur_lib;
	action todo;
	int as, ts;
   	bool all_zeros;

	// Class methods, to split up the monster call in 'progress'.
	bool start_in_landscape();
	bool jump_in_landscape();
	void get_modifier();
	void record_current();
	void update_landscape();
};

#endif
