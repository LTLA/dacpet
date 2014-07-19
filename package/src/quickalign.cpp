#include "quickalign.h"
//#define DEBUG 1

#ifdef DEBUG
#include <iostream>
#endif

/* An encoder class, to convert DNA sequences into k-mers. */

struct encoder {
	encoder(const char*);
	int get_code() const;
	int get_position() const;
	int advance ();
private:
	const char* target;
	int curstart, curcode;
	static int encode(const char);
	int mult;
};

encoder::encoder(const char* x) : target(x), curstart(0), curcode(0), mult(1) {
	for (int i=0; i<WORDSIZE; ++i) { 
		if (x[i]=='\0') { throw std::runtime_error("supplied string is less than the specified word size"); }
		curcode+=encode(x[i])*mult;
		mult*=4;
	}
	mult/=4;
	return;
}

int encoder::get_code() const { return curcode; }

int encoder::get_position() const { return curstart; }

int encoder::advance () {
	const char& current=target[curstart+WORDSIZE];
	if (current=='\0') { return 0; }
	curcode -= encode(target[curstart]);
	curcode /= 4;
	curcode += encode(current) * mult;
	++curstart;
	return 1;
}

int encoder::encode (const char x) {
	switch (x) {
		case 'A': case 'a': case 'N': case 'n': // N's are just treated as A's.
			return 0;
		case 'C': case 'c':
			return 1;
		case 'G': case 'g':
			return 2;
		case 'T': case 't':
			return 3;
		default:
			throw std::runtime_error("invalid character requested"); 
	}
}

/* Aligner class that performs k-mer assisted SW alignment. */
		
quickalign::quickalign(const char* ref, int ma, int mis, int gapo, int gape) : reflen(0), 
		reference(NULL), matched(ma), mismatch(mis), gapopen(gapo), gapext(gape) {
	encoder refen(ref);
	do {
		const int& current=refen.get_code();
		if (substrings.size() <= current) { substrings.resize(current+1, -1); }
		if (substrings[current] < 0) { 
			substrings[current]=matched_pos.size(); 
			matched_pos.resize(substrings[current]+1);
		}
		matched_pos[substrings[current]].push_back(refen.get_position());
	} while (refen.advance());	

	while (ref[reflen]!='\0') { ++reflen; }
	primary_score.resize(reflen);
	old_score.resize(reflen);
	ref_opened.resize(reflen);
	inc_opened.resize(reflen);
	old_inc_opened.resize(reflen);

	// Replacing the old copy with a new copy (throwing in the null for safety).
	reference = new char[reflen+1];
	for (int i=0; i<reflen; ++i) { reference[i]=std::toupper(ref[i]); }
	reference[reflen]='\0';
	return;
}

quickalign::~quickalign() {
	delete [] reference;
	return;
}

int quickalign::score_incoming(const char* target) {
	encoder incoming(target);
	int inclen=0;
	while (target[inclen]!='\0') { ++inclen; }

	// Identifies the longest matching subsequence, based on the words.
	std::deque<std::pair<int, int> > last_match(reflen, std::make_pair(-1, -1));
	int refstart, incstart, matchlen=0;
	do {
		const int& curcode=incoming.get_code();
		const int& curpos=incoming.get_position();
		if (curcode >= substrings.size() || substrings[curcode]<0) { continue; }

 		const std::deque<int>& refposes=matched_pos[substrings[curcode]];
		for (int r=refposes.size()-1; r>=0; --r) {
			const int& refpos=refposes[r];
			if (!refpos) {
				last_match[0]=std::make_pair(curpos, 1);
				if (!matchlen) { 
					refstart=refpos;
					incstart=curpos;
					matchlen=1;
				}
			} else {
				int& existing_len=last_match[refpos].second;
				if (last_match[refpos-1].first+1==curpos) {
					existing_len=last_match[refpos-1].second+1;
					if (existing_len > matchlen) {
						matchlen=existing_len;
						refstart=refpos-matchlen+1;
						incstart=curpos-matchlen+1;
					}
				} else {
					existing_len=1;
					if (!matchlen) { 
						refstart=refpos;
						incstart=curpos;
						matchlen=1;
					}
				}
				last_match[refpos].first=curpos;
			}
		}
	} while (incoming.advance());
	
	// Breaking if no subwords matched.
	if (!matchlen) { 
		return 0;
	}

	// Performs SW alignment on either end. 
	int actual_matchlen=WORDSIZE+matchlen-1;
	int base_score=actual_matchlen*matched;
	base_score += semilocal(refstart + actual_matchlen, reflen, target, incstart + actual_matchlen, inclen, true, false);
	base_score += semilocal(0, refstart, target, 0, incstart, false, true);

#ifdef DEBUG
	std::cout << "Choosing alignment score" << std::endl;
	std::cout << reference << " vs " << target << std::endl;
	std::cout << "\t" << std::string(reference + refstart).substr(0, actual_matchlen) << std::endl;
	std::cout << "\t" << std::string(target + incstart).substr(0, actual_matchlen) << std::endl;
	std::cout << "Base score is " << base_score << std::endl;
	std::cout << "####################" << std::endl;
#endif

	return base_score;
}

int quickalign::semilocal(int rstart, int rend,
						  const char* inc, int istart, int iend,
				  		  bool force_left, bool force_right) {
	// Taking a short cut if either strings are empty.
	const int nref=rend-rstart, ninc=iend-istart;
	if (nref==0 || ninc==0) { 
		if (force_left && force_right) { 
			if (nref || ninc) { 
				return (std::max(nref, ninc)-1)*gapext + gapopen;
			}
		} 
		return 0; 
	}

#ifdef DEBUG
	std::cout << "Forcing left " << force_left << " and right " << force_right << std::endl;
#endif

	/* Initializing starting values for those that need them. Imagine these
 	 * guys as the '-1' column of the DP matrix. Technically, old_inc_opened
 	 * shouldn't exist as it's impossible to have a opened gap before the start
 	 * of the sequence! Nonetheless, we give it a value so that the loops below
	 * will calculate the right value.
	 */
	if (force_left) { 
		old_score[0] = gapopen;
		old_inc_opened[0] = 2*gapopen - gapext;
		for (int j=1; j<nref; ++j) {
			old_score[j] = old_score[j-1]+gapext;
			old_inc_opened[j] = old_inc_opened[j-1]+gapext;
		}
	} else {
		old_score[0] = 0;
		old_inc_opened[0] = gapopen - gapext; // Ignoring gaps on the reference sequence.
		for (int j=1; j<nref; ++j) { 
			old_score[j] = 0;
			old_inc_opened[j] = old_inc_opened[0];
		}
		ref_opened[0] = gapopen; // Ignoring gaps on the incoming sequence (this is constant across the incoming sequence).
	}

#ifdef DEBUG
	std::cout << "Old starts are:" << std::endl;
	for (int j=0; j <nref; ++j) { 
			std::cout << old_score[j] << "/" << old_inc_opened[j] <<"\t";
	}
	std::cout << "\n" << std::endl;
#endif
	
	/* If force_right is used, maxscore won't be, so don't worry about it.
	 * If force_right is not used, then it's impossible to get a negative
	 * score as we could just choose to not align the right edge. So, a value
 	 * of zero is most appropriate. 
	 */
	int maxscore=0;
	
	/* Running through the DP matrix. 'ref_opened' refers to the best score 
	 * where every cell is the result of an opened gap on the reference sequence.
	 * 'inc_opened" refers to the best score where every cell results from an
	 * opened gap on the incoming sequence. These are necessary to distinguish
	 * between penalties for gap opening and extension.
	 */
	for (int i=0; i<ninc; ++i) {
		const char incbase=std::toupper(inc[istart+i]);
	
		// Setting up values at 'j=0', as 'j-1' makes no sense.	
		inc_opened[0] = std::max(old_inc_opened[0] + gapext, old_score[0] + gapopen);
		if (force_left) {
			ref_opened[0] = gapopen*2 + gapext*i; // Manhattan traversal onto the mat.
			primary_score[0] = std::max((incbase==reference[rstart] ? matched : mismatch) + (i!=0 ? gapopen + gapext * (i-1) : 0), // Corner is zero value, otherwise, indel.
					std::max(ref_opened[0], inc_opened[0]));
		} else {
			primary_score[0] = std::max((incbase==reference[rstart] ? matched : mismatch),
					std::max(ref_opened[0], inc_opened[0]));
			if (primary_score[0] < 0) { primary_score[0] = 0; }
		}
		if (!force_right && primary_score[0] > maxscore) { maxscore=primary_score[0]; }
#ifdef DEBUG
			std::cout << primary_score[0] << "(" << incbase << ", "<< reference[rstart] << ")";
#endif
				
		// Iterating across the remainders.
		for (int j=1; j<nref; ++j) {
			ref_opened[j]=std::max(ref_opened[j-1] + gapext, primary_score[j-1] + gapopen);
			inc_opened[j]=std::max(old_inc_opened[j] + gapext, old_score[j] + gapopen);
			primary_score[j]=std::max((incbase==reference[rstart+j] ? matched : mismatch) + old_score[j-1], 
					std::max(ref_opened[j], inc_opened[j]));

			if (!force_left && primary_score[j] < 0) { primary_score[j]=0; }
			if (!force_right && primary_score[j] > maxscore) { maxscore=primary_score[j]; }

#ifdef DEBUG
			std::cout << "\t" << primary_score[j] << "(" << incbase << ", "<< reference[rstart+j] << ")";
#endif
		}
#ifdef DEBUG		
		std::cout << std::endl;
#endif
		primary_score.swap(old_score);
		inc_opened.swap(old_inc_opened);
	}

	if (force_right) {
		return old_score[nref-1];
	} else {
		return maxscore;
	}
}

