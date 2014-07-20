//#define DEBUG 1

#ifdef DEBUG
#include <iostream>
#endif

#include "chia_counter.h"

chia_counter::chia_counter(SEXP pulled, bool issame, int filt, count_manager* ptr) : nlibs(0), filter(filt), 
		is_same_chr(issame), cmptr(ptr) {
	if (!isNewList(pulled)) { throw std::runtime_error("data on interacting PETs must be contained within a list"); }
	nlibs=LENGTH(pulled);
	asptrs.resize(nlibs);
	tsptrs.resize(nlibs);
   	aeptrs.resize(nlibs);
   	teptrs.resize(nlibs);
	num.resize(nlibs);
	indices.resize(nlibs); // Refers to the element currently IN (or has just been popped out of) the priority queue.
	modifier.resize(nlibs);
	
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
	
	if (cmptr->nlibs!=nlibs) { throw std::runtime_error("number of libraries is not consistent between runs"); }
	return;
}

void chia_counter::process () 
/**************************************
 * A function that brings everything together. The idea is to begin at the
 * start of 'landscape' (with jumps if necessary, if there is nothing
 * happening). We then determine if there's any modifications to be implemented
 * from the presence of incoming read pairs (or from the closure of existing
 * read pairs). We decide whether we want to record the underlying box. We then
 * figure out whether we need to update the landscape to reflect the mods.
 ************************************/
{
	while (start_in_landscape()) { 
		do {
			if (next.top().status==EXISTING) { // Deciding whether we should jump to a new spot, if there's no changes to implement.
			   	allzeros=true;
				for (int i=0; i<nlibs; ++i) { if (modifier[i]) { allzeros=false; break; } }
				if (allzeros) { 
					if (!jump_in_landscape()) { break; } // i.e., nothing more in 'next'.
				}
			}
			get_modifier();

#ifdef DEBUG
	std::cout << "\tRecorded: " << recorded << ", underlying: " << underlying << ", on existing: " << on_existing << std::endl;
#endif			
			if (!recorded) { 
				record_current(); 
				recorded=true; // set to true, so that we don't re-record the same box.
			}
			update_landscape();

#ifdef DEBUG
			std::cout << "### LAUNCH" << std::endl;
			for (std::map<int, std::pair<int, int> >::iterator itl2=landscape.begin(); itl2!=landscape.end(); ++itl2) {
				std::cout << itl2->first << " ("<< (itl2->second).first << ", " 
					<< (itl2->second).second << ")" << std::endl; 
			}
#endif
		} while (!anchor_is_finished); 
	}
	return;
}

bool chia_counter::start_in_landscape() {
	if (next.empty()) { return false; }
	cur_anchor=next.top().anchor;
	itl=landscape.begin();		
	if (itl!=landscape.end()) { 
		next.push(coord(cur_anchor, itl->first, (itl->second).second, EXISTING)); // Adding the first landscape element.
#ifdef DEBUG
		std::cout << "### Setting up new existing point for " << cur_anchor << ", " << itl->first << std::endl;
#endif
	} 
	underlying=-1;
	on_existing=false;
	recorded=true;
	anchor_is_finished=false;
	return true;
}

bool chia_counter::jump_in_landscape() {
	next.pop();
	if (next.empty() || next.top().anchor!=cur_anchor) { return false; }
	itl=landscape.lower_bound(next.top().target); // cannot equal the first one, as we just killed a landscape point in the queue.
	--itl;
	underlying=(itl->second).second;
	recorded=(underlying==-1);
	++itl;
#ifdef DEBUG
	std::cout << "\tJumping to new element at " << itl->first << " (queue is at "<<next.top().target <<")" << std::endl;
#endif
	if (itl!=landscape.end()) { next.push(coord(cur_anchor, itl->first, (itl->second).second, EXISTING)); }
	return true;
}

void chia_counter::get_modifier() 
/************************************
 * Here, we collect all incoming commands with exactly the same anchor/target
 * start points, so as to include all relevant points into the modifier. We add
 * a stop order (and a new start order) for each start order removed. 
 ************************************/ 
{
	cur_target=next.top().target;
	do {
		cur_lib=next.top().library;
		todo=next.top().status;
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
					int& curdex=indices[cur_lib];
					next.push(coord(asptrs[cur_lib][curdex], teptrs[cur_lib][curdex], cur_lib, STOP_ADD));
					next.push(coord(aeptrs[cur_lib][curdex], teptrs[cur_lib][curdex], cur_lib, STOP_RM));
					++curdex;
					if (curdex<num[cur_lib]) {
						next.push(coord(asptrs[cur_lib][curdex], tsptrs[cur_lib][curdex], cur_lib, START_ADD)); 
						next.push(coord(aeptrs[cur_lib][curdex], tsptrs[cur_lib][curdex], cur_lib, START_RM)); 
					}
				}
				++modifier[cur_lib];
				break;
			case STOP_RM:
				++modifier[cur_lib];
				break;
			case START_RM: 
				--modifier[cur_lib];
				break;
			case STOP_ADD:
				--modifier[cur_lib];
				break;
			case EXISTING:
				underlying=cur_lib;
				recorded=(underlying==-1);
				on_existing=true;
				break;
			default:
				throw std::runtime_error("invalid command specified");  
		}
		if (next.empty() || cur_anchor!=next.top().anchor) {
			anchor_is_finished=true; // Nothing more to be added at this anchor.
			break;
		}
	} while (cur_target==next.top().target);
	return;
}

void chia_counter::record_current() 
/******************************************************
 * Everytime we add a new point and the area on which we are adding to isn't an
 * empty stretch, we are obliged to report the combination and coordinates of
 * the original area before the modification. So, we just pop it out and that
 * goes into the storage containers. We pretty much do the same when we perform
 * deletions.
 ******************************************************/
{
	/* 'itl' only points to the overlapped/underlying element when cur_target is sitting on top of said element; otherwise,
	 * it points to the current EXISTING in 'next' which is the one right after the 'underlying' element of interest.
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
	for (int i=0; i<nlibs; ++i) { total+=cmptr->access(i, underlying); }
	if (filter <= total) {
#ifdef DEBUG
		std::cout << "\tTrying to record " << astart <<", " << aend << ", " << tstart << ", " << tend << " (";
		for (int i=0; i<nlibs; ++i) { std::cout << cmptr->access(i, underlying) << ","; }
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
	return;
}

void chia_counter::update_landscape () 
/***************************************
 * Applying the modification to the 'underlying' value, to potentially generate
 * a new combination.  If so, we need to decide whether or not we need to add
 * another point in the landscape to mark this as so. Alternatively, we could
 * just modify an existing point on the landscape if we're directly on top of
 * one. Some cleaning may be necessary to avoid redundant landscape points. 
 ***************************************/
{
	total=0;
	for (int i=0; i<nlibs; ++i) {
		int& current=(cmptr->back(i)=modifier[i]); 	
		if (underlying!=-1) { current+=cmptr->access(i, underlying); }
		total+=current;
		if (current < 0) { throw std::runtime_error("synchronization error, count cannot be negative"); } 
	}
#ifdef DEBUG
	std::cout << "\tCurrent values are: ";
	for (int i=0; i<nlibs; ++i) { std::cout << "\t" << cmptr->back(i); }
	std::cout << std::endl;
#endif

	/* Checking whether we need to save this element, or delete existing elements. We don't bother to
 	 * add new '-1' markers if the previous element is also a '-1' (or we're at the start of the landscape).
 	 * We also choose not to add markers if the new marker is the same as the old marker. For existing 
 	 * points, we erase it if it becomes the same as the previous marker or if it becomes a '-1' at the start.
 	 */
	const int& combodex=(total ? cmptr->advance() : -1);
	if (!on_existing) {
		if (itl!=landscape.begin()) {
			--itl;
			if ((itl->second).second!=combodex || (total && (itl->second).first!=cur_anchor)) {
				++itl;
				landscape.insert(itl, std::make_pair(cur_target, std::make_pair(cur_anchor, combodex)));
			} else { ++itl; }
		} else if (total) {
			landscape.insert(itl, std::make_pair(cur_target, std::make_pair(cur_anchor, combodex)));
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
 		 * No point doing so if there's no actual addition/subtraction commands present, though.
 		 */
		on_existing=false;
#ifdef DEBUG
		std::cout << "Currently sitting on ";
 		if (itl==landscape.end()) { std::cout << "the end"; 
		} else { std::cout << itl->first; } 
		std::cout << " with " << anchor_is_finished << std::endl;
#endif
		if (itl!=landscape.end() && !anchor_is_finished) { next.push(coord(cur_anchor, itl->first, (itl->second).second, EXISTING)); }
	}
	return;
}

const std::deque<int>& chia_counter::access_anchor() const { return all_anchors; }

const std::deque<int>& chia_counter::access_target() const { return all_targets; }

const std::deque<int>& chia_counter::access_index() const { return index; }

