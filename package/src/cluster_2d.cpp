#include "dacpet.h"
//#define DEBUG 0

#ifdef DEBUG
#include <iostream>
#endif

const int nothing=-1;

SEXP cluster_2d (SEXP start_a, SEXP start_t, SEXP end_a, SEXP end_t, SEXP tol) try {
	if (!isInteger(start_a) || !isInteger(start_t) || !isInteger(end_a) || !isInteger(end_t)) {
		throw std::runtime_error("anchor or target start and ends must be integer"); }
	const int npts=LENGTH(start_a);
	if (npts!=LENGTH(start_t) || npts!=LENGTH(end_a) || npts!=LENGTH(end_t)) { 
		throw std::runtime_error("lengths of coordinate vectors are not equal"); }
	if (!isInteger(tol) || LENGTH(tol)!=1) { 
		throw std::runtime_error("tolerance should be an integer scalar"); }

	const int width=asInteger(tol);
	const int* asptr=INTEGER(start_a),
		* aeptr=INTEGER(end_a),
		* tsptr=INTEGER(start_t),
		* teptr=INTEGER(end_t);
	
	// Setting up the output construct, which specifies the cluster ID of each point.
	SEXP output=PROTECT(allocVector(INTSXP, npts));
try {
	int* optr=INTEGER(output);
	
	/* We assume all points are sorted by start_a, start_t already. The idea is to
 	 * raster over the interaction space, assimilating all those which lie within 'tol'
 	 * of previous points. Overlaps with multiple points result in the formation of
 	 * synonyms which will be resolved at the end of the function.
 	 */
	std::map<int, std::pair<int, int> > landscape;
	std::map<int, std::pair<int, int> >::iterator itl, ito;
	std::set<std::pair<int, int> > synonyms;
	int ptdex=0, numids=0;

	while (ptdex < npts) {
		const int& cura=asptr[ptdex];
		const int& curt=tsptr[ptdex];
		const int& endt=teptr[ptdex];
		const int& enda=aeptr[ptdex];

		const int baset=curt-width;
		const int basea=cura-width;
		const int finisht=endt+width;
		itl=landscape.lower_bound(baset);
		int& myid=(optr[ptdex]=nothing);

		// Identifying the landscape element before or at the current target
		if (itl==landscape.end()) {
			; 
		} else {
			if ((itl->first)!=baset && itl!=landscape.begin()) {
				--itl;
				if ((itl->second).first > basea) { myid=(itl->second).second; }
			}
			
			// Running through everything in range, seeing if it overlaps anything else.
			while (itl!=landscape.end() && itl->first < finisht) {
 				if ((itl->second).first > basea) {
					const int& alternative=(itl->second).second;
					if (alternative==nothing) { 
						;
					} else if (myid==nothing) {
						myid=alternative;
					} else if (myid!=alternative) {
						if (alternative > myid) {
 	 					   	synonyms.insert(std::make_pair(myid, alternative));
						} else {
 	 					   	synonyms.insert(std::make_pair(alternative, myid));
						}
					}
 				}
				++itl;
			} 
		}
		if (myid==nothing) {
			// Incrementing to the next ID, if the current box doesn't overlap with anything.
			myid=numids;
			++numids;
		}

		/* Now, adding the current box to the landscape. We delete all points below it;
		 * these shouldn't affect synonym identification later as anything matched 
		 * by something with a higher anchor at the same target should also be matched
		 * by the current box i.e. all synonyms accounted for.
		 */
		if (itl==landscape.begin()) {
 		   itl=landscape.insert(itl, std::make_pair(endt, std::make_pair(enda, nothing)));
		} else {
			--itl; 
			while (itl!=landscape.begin() && itl->first > endt) { --itl; }
			if (itl->first==endt) { 
				;
			} else {
				const int& lastid=(itl->second).second;
				const int& lastheight=(itl->second).first;
				++itl;
				itl=landscape.insert(itl, std::make_pair(endt, std::make_pair(lastheight, lastid)));
			}
		}
	
		if (itl==landscape.begin()) {
			itl=landscape.insert(itl, std::make_pair(curt, std::make_pair(enda, myid)));
		} else {
			--itl;
			while (itl!=landscape.begin() && itl->first > curt) { landscape.erase(itl--); }
			if (itl->first==curt) {
				(itl->second).first=enda;
				(itl->second).second=myid;
			} else {
				++itl;
				itl=landscape.insert(itl, std::make_pair(curt, std::make_pair(enda, myid)));
			}
		}

//#ifdef DEBUG
//		std::cout << "#### New landscape!" << std::endl;
//		for (std::map<int, std::pair<int, int> >::const_iterator itx=landscape.begin(); itx!=landscape.end(); ++itx) {
//			std::cout << itx->first << ": " << (itx->second).first << " (" << (itx->second).second << ")" << std::endl;
//		}
//#endif		
		++ptdex;
	}
#ifdef DEBUG
	std::cout << "Synonyms are: " << std::endl;
	for (std::set<std::pair<int, int> >::const_iterator itx=synonyms.begin(); itx!=synonyms.end(); ++itx) {
		std::cout << itx->first << "\t" << itx->second << std::endl;;
	}
#endif


	// Resolving synonyms using a recursive-ish algorithm. It's a bit slow but at least it'll work.
	std::deque<int> newids(numids, nothing);
	numids=1;
	int curid=0;

	while (!synonyms.empty()) {
		std::priority_queue<int, std::deque<int>, std::greater<int> > next;
		std::set<std::pair<int, int> >::iterator itx=synonyms.begin();
		curid=itx->first;
		newids[curid]=numids;

		while (1) {
			do {
				next.push(itx->second);
				newids[itx->second]=numids;
				synonyms.erase(itx++);
			} while (itx!=synonyms.end() && curid==itx->first);

			while (!next.empty() && itx!=synonyms.end() && itx->first!=curid) {
				curid=next.top();
				itx=synonyms.lower_bound(std::make_pair(curid, nothing)); 
				next.pop();
			}
			if (next.empty()) { break; }
			
			if (itx==synonyms.end()) { break; }
		}
		++numids;
	}

	// Mopping up anything which doesn't have any synonyms.
	for (size_t i=0; i<newids.size(); ++i) {
		if (newids[i]==nothing) { 
			newids[i]=numids;
			++numids;
		}
	}

	// Now, going through and correcting the synonyms.
	for (int i=0; i<npts; ++i) { optr[i]=newids[optr[i]]; }
} catch (std::exception &e) {
	UNPROTECT(1);
	throw;
}
	UNPROTECT(1);
	return output;
} catch (std::exception& e) {
	return mkString(e.what());
}
