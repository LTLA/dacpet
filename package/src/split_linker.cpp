#include "quickalign.h"
#include <string>
#include <fstream>
#include <sstream>

#define R_NO_REMAP
#include "dacpet.h"

/* Interrupt eating functions, because this function might take a while. */

static void chkIntFn(void *dummy) { R_CheckUserInterrupt(); }

bool checkInterrupt() { return (R_ToplevelExec(chkIntFn, NULL) == FALSE); } // this will call the above in a top-level context so it won't longjmp-out of your context

/* Checking linker validity. */

int check_linker(const char* in) {
	int i=0;
	while (in[i]!='\0') {
		switch(in[i]) { 
			case 'A': case 'a': case 'C': case 'c': case 'G': case 'g': case 'T': case 't':
				break;
			default:
				throw std::runtime_error("invalid character found in the linker sequence"); 
		}
		++i;
	}
	return i;
}

int find_mode(const int& A, const int& B, const int& minscore) {
	if (A > B && A >= minscore) { return 0; }
	else if (A <B && B >= minscore) { return 1; }
	return -1;
}

SEXP split_linkers(SEXP file1, SEXP file2, SEXP linkerA, SEXP linkerB, 
		SEXP prefix, SEXP mscore, SEXP starting, SEXP readfix) try {
	// Checking input files.
	if (!Rf_isString(file1) || LENGTH(file1)!=1 || !Rf_isString(file2) || LENGTH(file2)!=1) { 
		throw std::runtime_error("each FastQ file path must be a single character string"); 
	} 
	const char* infile1=CHAR(STRING_ELT(file1, 0));
	const char* infile2=CHAR(STRING_ELT(file2, 0));

	// Checking linkers.
	if (!Rf_isString(linkerA) || LENGTH(linkerA)!=1 || !Rf_isString(linkerB) || LENGTH(linkerB)!=1) { 
		throw std::runtime_error("each linker sequence must be a single character string"); 
	} 
	const char* linkA=CHAR(STRING_ELT(linkerA, 0));
	const char* linkB=CHAR(STRING_ELT(linkerB, 0));
	const int nA=check_linker(linkA);
	const int nB=check_linker(linkB);
	if (nA!=nB) {  throw std::runtime_error("length of linker sequences should really be equal"); }

	// Checking output prefixes.
	if (!Rf_isString(prefix) || LENGTH(prefix)!=1) {
		throw std::runtime_error("file prefix should be a single character string"); 
	}
	const char* oprefix=CHAR(STRING_ELT(prefix, 0));
	const char* newprefix=NULL;
	if (Rf_isString(readfix) && LENGTH(readfix)==1) { 
		newprefix=CHAR(STRING_ELT(readfix, 0));
	} else if (LENGTH(readfix) > 1) {
		throw std::runtime_error("read prefix should be a single character string"); 
	}

	// Checking minimum score, starting position.
	if (!Rf_isInteger(starting) || LENGTH(starting)!=1) { throw std::runtime_error("starting position should be an integer scalar"); }
	const int startpos=Rf_asInteger(starting);
	if (!Rf_isInteger(mscore) || LENGTH(mscore)!=1) { throw std::runtime_error("minimum alignment score should be an integer scalar"); }
	const int minscore=Rf_asInteger(mscore);

	// Opening all input files.
	std::ifstream input_1(infile1), input_2(infile2);
	{
		std::stringstream err;
		if (! input_1.good()) { 
 	    	err << "failed to read from `" << infile1 << "'";
			throw std::runtime_error(err.str().c_str());
		} else if (! input_2.good()) { 
			err << "failed to read from `" << infile2 << "'";
			throw std::runtime_error(err.str().c_str());
		}
	}

	// Setting up output file names.
	std::string prefname(oprefix);
	std::ofstream output_AA_1((prefname+"_AA_1.fastq").c_str());
	std::ofstream output_AB_1((prefname+"_AB_1.fastq").c_str());
	std::ofstream output_BB_1((prefname+"_BB_1.fastq").c_str());
	std::ofstream output_other_1((prefname+"_other_1.fastq").c_str());
	std::ofstream output_AA_2((prefname+"_AA_2.fastq").c_str());
	std::ofstream output_AB_2((prefname+"_AB_2.fastq").c_str());
	std::ofstream output_BB_2((prefname+"_BB_2.fastq").c_str());
	std::ofstream output_other_2((prefname+"_other_2.fastq").c_str());
	if (output_AA_1.is_open() && output_AB_1.is_open() && output_BB_1.is_open() && output_other_1.is_open() &&
		output_AA_2.is_open() && output_AB_2.is_open() && output_BB_2.is_open() && output_other_2.is_open() ) {
	} else {
		throw std::runtime_error("failed to open output files, check output prefix"); 
	}
	std::ofstream* optr1=NULL, *optr2=NULL;

	/* Setting up the alignment. */
	quickalign align_A(linkA);
	quickalign align_B(linkB);

	std::string line1, line2, header1, header2;
	bool left, right;
	int counter=0, storage_mode, all_count=0;
	int sA1=startpos, sB1=startpos, sA2=startpos, sB2=startpos;
	while (1) {
		left=getline(input_1, line1);
		right=getline(input_2, line2);
		if (left!=right) {
 		    throw std::runtime_error("unbalanced numbers of reads in the paired FastQ files");
 		} else if (!left) { 
			if (counter!=0) { throw std::runtime_error("FastQ files are truncated before EOF"); }
			break; 
		}

		if (counter==0) {
			header1=line1;
			header2=line2;
			++all_count;
		} else {
			if (counter==1) {
				// Taking the last subsequence and assuming it's the linker.
				if (startpos<0) { 
					sA1 = line1.length() - nA;
					sB1 = line1.length() - nB;
					sA2 = line2.length() - nA;
					sB2 = line2.length() - nB;
				} else {
					sA1=sA2=sB1=sB2=startpos;
				}
				int left_store = find_mode(	align_A.score_incoming(line1.substr(sA1, nA).c_str()),
											align_B.score_incoming(line1.substr(sB1, nB).c_str()),
											minscore);
				int right_store =find_mode(	align_A.score_incoming(line2.substr(sA1, nA).c_str()),
											align_B.score_incoming(line2.substr(sB1, nB).c_str()),
											minscore);
				if (left_store < 0 || right_store < 0) {
					storage_mode=4;
				} else {
					storage_mode=left_store+right_store*2;
				}
			}

			switch(storage_mode) {
				case 0:
					optr1=&output_AA_1;
					optr2=&output_AA_2;
					break;
				case 1:
					optr1=&output_AB_1;
					optr2=&output_AB_2;
					break;
				case 2:
					optr1=&output_AB_1;
					optr2=&output_AB_2;
					break;
				case 3:
					optr1=&output_BB_1;
					optr2=&output_BB_2;
					break;
				case 4:
					optr1=&output_other_1;
					optr2=&output_other_2;
					break;
				default:
					throw std::runtime_error("invalid storage mode specified"); 
			}

			if (counter==1) {
				if (newprefix==NULL) {
					(*optr1) << header1 << std::endl;
					(*optr2) << header2 << std::endl;
				} else {
					(*optr1) << "@" << newprefix << "." << all_count << std::endl;
					(*optr2) << "@" << newprefix << "." << all_count << std::endl;
				}
			} 
			if (counter==2) {
				(*optr1) << line1 << std::endl;
				(*optr2) << line2 << std::endl;
			} else {
				// Editing the seuqnece and quality stromgs
				(*optr1) << line1.substr(0, sA1) << std::endl;
				(*optr2) << line2.substr(0, sA2) << std::endl;
			}
		}
		if (counter==3) { counter=0; }
		else { ++counter; }

		// Checking whether the user has requested an interrupt.
		if (all_count % 1000 == 0) { 
			if (checkInterrupt()) { throw std::runtime_error("interrupt requested"); }
		}
	}

	/* This is a function that has no purpose but to run through a bunch of possibilities
	 * in order to check the aligner function. Probably good to set #define DEBUG 1 in
	 * quickalign.cpp to examine the output when doing this.
	quickalign rlign_A("GTTGGAATGTATATCG");
	
	Rprintf("Perfection:\n");
	rlign_A.score_incoming("GTTGGAATGTATATCG"); 

	Rprintf("Mismatch at front:\n");
	rlign_A.score_incoming("GTGGGAATGTATATCG"); 
	Rprintf("Mismatch at back:\n");
	rlign_A.score_incoming("GTTGGAATGTAAATCG"); 
	Rprintf("Wrong linker:\n");
	rlign_A.score_incoming("GTTGGATAAGATATCG"); 
	
	Rprintf("Deletion at front:\n");
	rlign_A.score_incoming("GTTGAATGTATATCG");
	Rprintf("Deletion in middle:\n");
	rlign_A.score_incoming("GTTGGAAATATCG"); 
	Rprintf("Deletion in back:\n");
	rlign_A.score_incoming("GTTGGAATGTATATG");
	
	Rprintf("Insertion in front:\n");
	rlign_A.score_incoming("GTTGCGAATGTATATCG"); 
	Rprintf("Insertion in middle:\n");
	rlign_A.score_incoming("GTTGGAATGAGTATATCG");
	Rprintf("Insertion at back:\n");
	rlign_A.score_incoming("GTTGGAATGTATATCCG"); 
	*/	
	
	return Rf_ScalarInteger(0);
} catch (std::exception& e) {
	return Rf_mkString(e.what());
}


void checkery () {
	return;
}
 
