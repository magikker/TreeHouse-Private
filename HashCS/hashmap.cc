///*****************************************************/
//
// Copyright (C) 2006-2009 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS IMPLEMENTATION
//		HashMap: Class for hashing bit
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details (www.gnu.org).
//
//*****************************************************************/

#include "hashmap.hh"
#include "hashfunc.hh"

// For populating hash table
// It does not store bitstrings.
void
HashMap::hashing_bs_MJ(
    unsigned long long hv1,
    unsigned long long hv2,
    unsigned int res_numtrees) 
{
    if (_hashtab_hashCS[hv1].empty()) { 
        BUCKET_STRUCT_T_HASHCS bk;
        bk._hv2 = hv2;
//        bk._bs = NULL;
        bk._bs2 = NULL;
        bk._count = 1;
        _hashtab_hashCS[hv1].push_back(bk);
    } else {
        bool found = false;
        for (unsigned int i = 0; i < _hashtab_hashCS[hv1].size(); ++i) {
            if (_hashtab_hashCS[hv1][i]._hv2 == hv2) {
                found = true;
                if (_hashtab_hashCS[hv1][i]._count <= res_numtrees) { // filter for increasement
                    _hashtab_hashCS[hv1][i]._count++; // Increase the occurrence count
                }
                break;
            }
        }
        if (!found) {
            BUCKET_STRUCT_T_HASHCS bk;
            bk._hv2 = hv2;
//            bk._bs = NULL;
            bk._bs2 = NULL;
            bk._count = 1;
            _hashtab_hashCS[hv1].push_back(bk);
        }
    }
}

// For comparing hv1 and hv2 and storing bitstrings.
//void
//HashMap::hashing_bs_MJ(
//    BitSet& bs,
//    unsigned int numTaxa,
//    unsigned long long hv1,
//    unsigned long long hv2,
//    vector<BitSet*> & vec_bs,
//    unsigned int res_numtrees) 
//{
//    if (!_hashtab_hashCS[hv1].empty()) {
//        bool found = false;
//        for (unsigned int i = 0; i < _hashtab_hashCS[hv1].size(); ++i) {
//            if (_hashtab_hashCS[hv1][i]._hv2 == hv2) {
//                if (_hashtab_hashCS[hv1][i]._count <= res_numtrees) { 
//                    _hashtab_hashCS[hv1][i]._count++; // Increase the occurrence count
//                    
//                    // now we find a majority bipartition
//                    // save the bitstring in vec_bs
//                    if (_hashtab_hashCS[hv1][i]._count > res_numtrees) { 
//                        BitSet* tempbs = new BitSet(numTaxa);
//                        *tempbs = bs;
//                        vec_bs.push_back(tempbs);
//                    }
//                }
//                found = true;
//                break;
//            }
//        }
//    }
//}

// For comparing hv1 and hv2 and storing bitstrings.
void
HashMap::hashing_bs_MJ(
    bool * bs,
    unsigned int numTaxa,
    unsigned long long hv1,
    unsigned long long hv2,
    vector<bool *> &vec_bs,
    unsigned int res_numtrees) 
{
    if (!_hashtab_hashCS[hv1].empty()) {
        bool found = false;
        for (unsigned int i = 0; i < _hashtab_hashCS[hv1].size(); ++i) {
            if (_hashtab_hashCS[hv1][i]._hv2 == hv2) {
                if (_hashtab_hashCS[hv1][i]._count <= res_numtrees) { 
                    _hashtab_hashCS[hv1][i]._count++; // Increase the occurrence count
                    
                    // now we find a majority bipartition
                    // save the bitstring in vec_bs
                    if (_hashtab_hashCS[hv1][i]._count > res_numtrees) { 
//                        BitSet* tempbs = new BitSet(numTaxa);
		      bool *tempbs = new bool[numTaxa];
		      for (unsigned int i = 0; i < numTaxa; i++)
                        tempbs[i] = bs[i];
                        vec_bs.push_back(tempbs);
                    }
                }
                found = true;
                break;
            }
        }
    }
}



// if the first tree, just collect the bipartition
//void
//HashMap::hashing_bs_SC(
//    BitSet& bs,
//    unsigned int numTaxa,
//    unsigned long long hv1,
//    unsigned long long hv2) 
//{    		
//		
//		BUCKET_STRUCT_T_HASHCS bk;
//		bk._hv2 = hv2;
//		BitSet* tempbs = new BitSet(numTaxa);
//		*tempbs = bs;
//		bk._bs = tempbs;
//		bk._count = 1;
//		_hashtab_hashCS[hv1].push_back(bk);
//}


void
HashMap::hashing_bs_SC(
    bool * bs,
    unsigned int numTaxa,
    unsigned long long hv1,
    unsigned long long hv2) 
{    				
		BUCKET_STRUCT_T_HASHCS bk;
		bk._hv2 = hv2;
//		BitSet* tempbs = new BitSet(numTaxa);
		bool *tempbs = new bool[numTaxa];
		for (unsigned int i = 0; i < numTaxa; i++)
		  tempbs[i] = bs[i];
		bk._bs2 = tempbs;
		bk._count = 1;
		_hashtab_hashCS[hv1].push_back(bk);
}


// to process the other trees than the 1st tree.
//void
//HashMap::hashing_bs_SC(
//    unsigned int treeIdx,
//    unsigned int numTaxa,
//    unsigned int numTrees,
//    unsigned long long hv1,
//    unsigned long long hv2,
//    vector<BitSet*> & vec_bs) 
//{
// 
//    // The other trees: just compare the biaprtitions in the hashtable
//    if ((treeIdx > 0) & (treeIdx < numTrees - 1)) { 
//        if (!_hashtab_hashCS[hv1].empty()) { // if the hash table location has already been occupied
//	
//			// get the number of items in the linked list			
//            for (unsigned int i = 0; i < _hashtab_hashCS[hv1].size(); ++i) { 
//                if (_hashtab_hashCS[hv1][i]._hv2 == hv2) { // if the same bp is found
//                    _hashtab_hashCS[hv1][i]._count++;
//                    break;
//                }
//            }
//        }
//    } 
//    // if treeIdx == numTrees ?????????????????
//    else { 
//        if (!_hashtab_hashCS[hv1].empty()) {
//            for (unsigned int i = 0; i < _hashtab_hashCS[hv1].size(); ++i) {
//                if (_hashtab_hashCS[hv1][i]._hv2 == hv2) {
//                	
//                		// now we find strict bipartition
//                		// store the bitstring in vec_bs
//                    if (_hashtab_hashCS[hv1][i]._count == numTrees - 1)
//                        vec_bs.push_back(_hashtab_hashCS[hv1][i]._bs);
//                    break;
//                }
//            }
//        }
//    }
//}

void
HashMap::hashing_bs_SC(
    unsigned int treeIdx,
    unsigned int numTaxa,
    unsigned int numTrees,
    unsigned long long hv1,
    unsigned long long hv2,
    vector<bool *> &vec_bs) 
{
 
    // The other trees: just compare the biaprtitions in the hashtable
    if ((treeIdx > 0) & (treeIdx < numTrees - 1)) { 
        if (!_hashtab_hashCS[hv1].empty()) { // if the hash table location has already been occupied
	
			// get the number of items in the linked list			
            for (unsigned int i = 0; i < _hashtab_hashCS[hv1].size(); ++i) { 
                if (_hashtab_hashCS[hv1][i]._hv2 == hv2) { // if the same bp is found
                    _hashtab_hashCS[hv1][i]._count++;
                    break;
                }
            }
        }
    } 
    // if treeIdx == numTrees ?????????????????
    else { 
        if (!_hashtab_hashCS[hv1].empty()) {
            for (unsigned int i = 0; i < _hashtab_hashCS[hv1].size(); ++i) {
                if (_hashtab_hashCS[hv1][i]._hv2 == hv2) {
                	
                		// now we find strict bipartition
                		// store the bitstring in vec_bs
                    if (_hashtab_hashCS[hv1][i]._count == numTrees - 1)
                        vec_bs.push_back(_hashtab_hashCS[hv1][i]._bs2);
                    break;
                }
            }
        }
    }
}


void
HashMap::HashMap_clear() 
{
    for (unsigned int i = 0; i < _hashtab_hashCS.size(); ++i) {
        for (unsigned int j = 0; j < _hashtab_hashCS[i].size(); ++j) {
//            if (_hashtab_hashCS[i][j]._bs) {
//                delete _hashtab_hashCS[i][j]._bs;
//                _hashtab_hashCS[i][j]._bs = NULL;
//            }
            if (_hashtab_hashCS[i][j]._bs2) {    
                delete [] _hashtab_hashCS[i][j]._bs2;
                _hashtab_hashCS[i][j]._bs2 = NULL;
            }
            _hashtab_hashCS[i].clear();
        }
    }

    _hashtab_hashCS.clear();
}

void
HashMap::UHashfunc_init(
    unsigned int t, // number of trees
    unsigned int n, // number of taxa
    unsigned int c) // c value in m2 > c*t*n
{
    _HF.UHashfunc_init(t, n, c);
}


void
HashMap::UHashfunc_init(
    unsigned int t, // number of trees
    unsigned int n, // number of taxa
    unsigned int c, // c value in m2 > c*t*n
    int32 newseed)  // user specified seed
{
    _HF.UHashfunc_init(t, n, c, newseed);
}



