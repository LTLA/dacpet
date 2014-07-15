#ifndef QUICKALIGN_H
#define QUICKALIGN_H

#include <deque> 
#include <cctype>
#include <stdexcept>

/* An alignment class, that matches substrings, tries to expand the alignments
 * as much as possible, and then performs SW alignment on the rest.
 */

class quickalign {
public:
	quickalign(const char*, int=1, int=-1, int=-3, int=-1);
	~quickalign();
	int score_incoming(const char*);
private:
	char* reference;
	int reflen;
	
	std::deque<int> substrings;
	std::deque<std::deque<int> > matched_pos;

	const int matched, mismatch, gapopen, gapext;
	int semilocal(int, int, const char*, int, int, bool, bool);
	std::deque<int> primary_score, old_score, 
		ref_opened, inc_opened, old_inc_opened;
};

const int WORDSIZE=6;
#endif
