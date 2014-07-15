#include "count_manager.h"
#include "dacpet.h"
//#define DEBUG 1

#ifdef DEBUG
#include <iostream>
#endif

enum action { START_ADD, 
	STOP_ADD, 
	START_RM, 
	STOP_RM, 
	EXISTING }; // EXISTING must be last for some loop conditions to hold.

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

SEXP count_chia(SEXP extptr, SEXP pulled, SEXP filter, SEXP issame) try {
   	// Launching pointers.
   	if (!isInteger(filter) || LENGTH(filter)!=1) { throw std::runtime_error("filter value should be an integer scalar"); }
   	if (!isLogical(issame) || LENGTH(issame)!=1) { throw std::runtime_error("same chromosome specification should be a logical scalar"); }
	const int filt=asInteger(filter);
	const bool is_same_chr=asLogical(issame);
	if (!isNewList(pulled)) { throw std::runtime_error("data on interacting PETs must be contained within a list"); }
	const int nlibs=LENGTH(pulled);
	std::deque<const int*> asptrs(nlibs), tsptrs(nlibs), aeptrs(nlibs), teptrs(nlibs);
	std::deque<int> num(nlibs), indices(nlibs);
	std::priority_queue<coord, std::deque<coord>, std::greater<coord> > next;

    for (int i=0; i<nlibs; ++i) {
		SEXP current=VECTOR_ELT(pulled, i);
		if (!isNewList(current) || LENGTH(current)!=4) { 
			throw std::runtime_error("interactions must be supplied as a data.frame with anchor and target starts and ends"); }
		
		for (int j=0; j<4; ++j) {
			SEXP current_col=VECTOR_ELT(current, j);
			if (!isInteger(current_col)) { throw std::runtime_error("interaction data must be in integer format"); }
			int* ptr=INTEGER(current_col);
			switch (j) {
				case 0: 
					asptrs[i]=ptr; 
					num[i]=LENGTH(current_col);
					break;
				case 1: tsptrs[i]=ptr; break;
				case 2: aeptrs[i]=ptr; break;
				case 3: teptrs[i]=ptr; break;
				default: break;
			}
		}
		
		// Populating the priority queue, with support for pan-diagonal elements.
        if (num[i]) { 
			next.push(coord(asptrs[i][0], tsptrs[i][0], i, START_ADD));
			next.push(coord(aeptrs[i][0], tsptrs[i][0], i, START_RM));
 		}	
    }

	// Loading the count combinations.
	count_manager* xptr=(count_manager*)R_ExternalPtrAddr(extptr);
	if (xptr==NULL) { throw std::runtime_error("external pointer cannot point to an invalid address"); }
	count_manager& target=*xptr;
	if (target.nlibs!=nlibs) { throw std::runtime_error("number of libraries is not consistent between runs"); }

	/******************************************************
 	 *
 	 * Running through all of them and returning count combinations for blocks of windows as they become available and 
 	 * resolved. Everything is already sorted, so the idea is just to run through all lists, introduce points into the 
 	 * landscape as <target, current combo index> (which is possible, as all combinations are guaranteed to exist).
 	 *   Everytime we add a new point and the area on which we are adding to isn't an empty stretch, we are obliged to 
 	 * report the combination and coordinates of the original area before the modification. So, we just pop it out and 
 	 * that goes into the storage containers. We pretty much do the same when we perform deletions.
 	 *
 	 ******************************************************/

	std::map<int, std::pair<int, int> > landscape;
    std::map<int, std::pair<int, int> >::iterator itl;
	std::deque<int> modifier(nlibs);
	int underlying, total;
	int as, ts;
	bool failed, recorded, on_existing;
	std::deque<int> all_anchors, all_targets, index;

	while (!next.empty()) {
		const int cur_anchor=next.top().anchor;
		itl=landscape.begin();		
		if (itl!=landscape.end()) { 
			next.push(coord(cur_anchor, itl->first, (itl->second).second, EXISTING)); 
#ifdef DEBUG
			std::cout << "### Setting up new existing point for " << cur_anchor << ", " << itl->first << std::endl;
#endif
		} 
		underlying=-1;
		on_existing=false;
		recorded=true;
		failed=false;
		
		/* Inner loop, cycling over all points with the same anchor values. If the first
 		 * thing in the queue is an EXISTING element and there are no modifications to be applied
 		 * to it, we kill it and we jump to something that does have modifications lining up for it.
 		 */
		do {
			if (next.top().status==EXISTING) {
   				bool allzeros=true;
				for (int i=0; i<nlibs; ++i) { if (modifier[i]) { allzeros=false; break; } }
				if (allzeros) { 
					next.pop();
					if (next.empty() || next.top().anchor!=cur_anchor) { break; }
					itl=landscape.lower_bound(next.top().target); // cannot equal the first one, as we just killed a landscape point in the queue.
					--itl;
					underlying=(itl->second).second;
					recorded=(underlying==-1);
					++itl;
#ifdef DEBUG
					std::cout << "\tJumping to new element at " << itl->first << " (queue is at "<<next.top().target <<")" << std::endl;
#endif
					if (itl!=landscape.end()) { next.push(coord(cur_anchor, itl->first, (itl->second).second, EXISTING)); }
				}
			}
			const int curt=next.top().target;

			/* Inner-inner loop; processing those with exactly the same anchor/target start points, so as 
			 * to include all relevant points into this combination (as it will be non-modifiable past here).
			 * We add a stop order (and a new start order) for each start order removed.
			 */ 
			do {
				const int lib=next.top().library;
				const action todo=next.top().status;
#ifdef DEBUG
				std::cout << "Processing: " << next.top().anchor << "\t" << next.top().target << "\t" << next.top().library << "\t" << next.top().status <<std::endl;
#endif
				next.pop();
				switch (todo) {
					case START_ADD: 
						{
							/* This is the only safe place to add actions, as inputs are sorted by 
							 * anchor and target start; this is the only place where new additions
							 * are guaranteed to be after elements already in the queue.
							 */
							int& curdex=indices[lib];
							next.push(coord(asptrs[lib][curdex], teptrs[lib][curdex], lib, STOP_ADD));
							next.push(coord(aeptrs[lib][curdex], teptrs[lib][curdex], lib, STOP_RM));
							++curdex;
							if (curdex<num[lib]) {
								next.push(coord(asptrs[lib][curdex], tsptrs[lib][curdex], lib, START_ADD)); 
								next.push(coord(aeptrs[lib][curdex], tsptrs[lib][curdex], lib, START_RM)); 
							}
						}
						++modifier[lib];
						break;
					case STOP_RM:
						++modifier[lib];
						break;
					case START_RM: 
						--modifier[lib];
						break;
					case STOP_ADD:
						--modifier[lib];
						break;
					case EXISTING:
						underlying=lib;
						recorded=(underlying==-1);
						on_existing=true;
						break;
					default:
						throw std::runtime_error("invalid command specified");  
				}
				if (next.empty() || cur_anchor!=next.top().anchor) {
					failed=true;
					break;
				}
			} while (curt==next.top().target);

#ifdef DEBUG
			std::cout << "\tRecorded: " << recorded << ", underlying: " << underlying << ", on existing: " << on_existing << std::endl;
#endif			
			if (!recorded) {
				/* Storing the component that needs to be recorded (note that itl only points to 
				 * the overlapped element when curt is sitting on top of said element; otherwise,
		 		 * itl points to the current EXISTING element in the queue which is the next
		 		 * element from the perspective of the 'underlying' element).
		 		 */
				if (!on_existing) { 
					if (itl==landscape.begin()) { throw std::runtime_error("synchronisation fault, attempting to find element before start");} 					
					--itl; 
				}
				if (itl==landscape.end()) { throw std::runtime_error("synchronisation fault, attempting to use end element"); }
				int astart=(itl->second).first;
				const int& tstart=itl->first;
				const int& aend=cur_anchor;
				if (astart>=aend) { throw std::runtime_error("synchronisation fault, attempting to record an empty box"); }								
				(itl->second).first=cur_anchor; // As the underlying box is now recorded, so we don't want to count it twice.
				++itl;
				if (itl==landscape.end()) { throw std::runtime_error("synchronisation fault, attempting to find element after end"); }
				const int& tend=itl->first;
				if (on_existing) { --itl; }

				// Deciding whether to actually add it.
 			    total=0;
				for (int i=0; i<nlibs; ++i) { total+=target.access(i, underlying); }
				if (filt <= total) {
#ifdef DEBUG
					std::cout << "\tTrying to record " << astart <<", " << aend << ", " << tstart << ", " << tend << " (";
					for (int i=0; i<nlibs; ++i) { std::cout << target.access(i, underlying) << ","; }
					std::cout << ")" << std::endl;

#endif					
					if (!is_same_chr || astart >= tend) {
						for (as=astart; as < aend; ++as) {
							for (ts=tstart; ts < tend; ++ts) {
								index.push_back(underlying);
								all_anchors.push_back(as);
								all_targets.push_back(ts);
							}
						}
					} else if (aend > tstart) {
						/* If a box runs across a diagonal in a same-chromosome comparison, we split it up so that anything 
 						 * past the diagonal is removed. Technically, additions (and thus, boxes) should be reflected across 
 						 * the diagonal. However, as the anchor end is guaranteed to be at least as big as the target end, 
 						 * any reflection will already have been covered by additions which are already present. So, no need.
	 	 	 	 	 	 */
						if (tstart > astart) { astart=tstart; }
						for (as=astart; as < aend; ++as) {
							for (ts=tstart; ts < tend; ++ts) {
								if (ts > as) { break; }
								index.push_back(underlying);
								all_anchors.push_back(as);
								all_targets.push_back(ts);
							}
						}
					}
				}
				recorded=true;
			}

			// Applying the modification to the current point.
			total=0;
			if (underlying==-1) {
				for (int i=0; i<nlibs; ++i) { 
					target.back(i)=modifier[i]; 
					total+=target.back(i);
					if (target.back(i)<0) { throw std::runtime_error("synchronization error, count cannot be negative"); 
					}
				}
			} else {
				for (int i=0; i<nlibs; ++i) { 
					target.back(i)=modifier[i]+target.access(i, underlying); 
					total+=target.back(i);
					if (target.back(i)<0) { throw std::runtime_error("synchronization error, count cannot be negative"); }
				}
			}
#ifdef DEBUG
			std::cout << "\tCurrent values are: ";
			for (int i=0; i<nlibs; ++i) { std::cout << "\t" << target.back(i); }
			std::cout << std::endl;
#endif

			/* Checking whether we need to save this element, or delete existing elements. We don't bother to
 			 * add new '-1' markers if the previous element is also a '-1' (or we're at the start of the landscape).
 			 * We also choose not to add markers if the new marker is the same as the old marker. For existing 
 			 * points, we erase it if it becomes the same as the previous marker or if it becomes a '-1' at the start.
 			 */
			const int& combodex=(total ? target.advance() : -1);
			if (!on_existing) {
				if (itl!=landscape.begin()) {
					--itl;
					if ((itl->second).second!=combodex || (total && (itl->second).first!=cur_anchor)) {
						++itl;
						landscape.insert(itl, std::make_pair(curt, std::make_pair(cur_anchor, combodex)));
					} else { ++itl; }
				} else if (total) {
					landscape.insert(itl, std::make_pair(curt, std::make_pair(cur_anchor, combodex)));
				}
			} else {
				(itl->second).first=cur_anchor;
				(itl->second).second=combodex;
				if (itl==landscape.begin()) {
					if (!total) { landscape.erase(itl++); } 
					else { ++itl; }
				} else {
					--itl;
					if (combodex==(itl->second).second && (!total || cur_anchor==(itl->second).first)) { 
						++itl;
						landscape.erase(itl++);
					} else { ++itl; ++itl; }
				}

				/* Need to re-add the point in the priority queue (note the free increment, above).
 				 * No point doing so if we're just going to break out again, though.
 				 */
				on_existing=false;
#ifdef DEBUG
				std::cout << "Currently sitting on ";
 			    if (itl==landscape.end()) { std::cout << "the end"; 
				} else { std::cout << itl->first; } 
				std::cout << " with " << failed << std::endl;
#endif
				if (itl!=landscape.end() && !failed) { 	next.push(coord(cur_anchor, itl->first, (itl->second).second, EXISTING)); }
			}

#ifdef DEBUG
		// Printing stuff at every anchor cycle, to have a look at what's going on.
		std::cout << "### LAUNCH" << std::endl;
		for (std::map<int, std::pair<int, int> >::iterator itl2=landscape.begin(); itl2!=landscape.end(); ++itl2) {
			std::cout << itl2->first << " ("<< (itl2->second).first << ", " 
				<< (itl2->second).second << ")" << std::endl; 
		}
#endif
		} while (!failed);
	}

	/*****************************************************
 	 *
 	 * Filling up output. We return, at least from this function, a list
 	 * of integer vector specifying the anchor start/end, target start/end
 	 * and the index of the count combinations. Note that the actual
 	 * count combinations are only returned later.
 	 *
 	 ******************************************************/

	SEXP output=PROTECT(allocVector(VECSXP, 3));
	try {
		const int totaln=index.size();
		for (int i=0; i<3; ++i) {
			SET_VECTOR_ELT(output, i, allocVector(INTSXP, totaln));
			int* optr=INTEGER(VECTOR_ELT(output, i));
			const std::deque<int> * curcopy=NULL;
			switch(i) {
				case 0:
					curcopy=&all_anchors;
					break;
				case 1:
					curcopy=&all_targets;
					break;
				case 2:
					curcopy=&index;
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
