//*****************************************************************/
//
// Copyright (C) 2006-2009 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS DEFINITION
//		HashMap: Class for hashing
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

#ifndef HASHMAP_HH
#define HASHMAP_HH

#include <vector>

//#include "bitset.hh"
#include "hashfunc.hh"

#include <bitset>

using namespace std;


// For HashConsense

typedef struct {
    unsigned long long _hv2;
//    BitSet *_bs;
    bool *_bs2;
    unsigned long _count;
} BUCKET_STRUCT_T_HASHCS;

typedef vector<BUCKET_STRUCT_T_HASHCS> V_BUCKET_T_HASHCS;
typedef vector<V_BUCKET_T_HASHCS> HASHTAB_T_HASHCS;

class HashMap {
public:

    HashMap() {
    }

    ~HashMap() {
    }

    HASHTAB_T_HASHCS _hashtab_hashCS;
    CHashFunc _HF;

    void hashing_bs_MJ(unsigned long long hv1, unsigned long long hv2, unsigned int res_numtrees);
//    void hashing_bs_MJ(BitSet& bs, unsigned int numTaxa, unsigned long long hv1, unsigned long long hv2, vector<BitSet*> & vec_bs, unsigned int res_numtrees);
    void hashing_bs_MJ(bool* bs, unsigned int numTaxa, unsigned long long hv1, unsigned long long hv2, vector<bool *> & vec_bs, unsigned int res_numtrees);
    
//    void hashing_bs_SC(BitSet& bs, unsigned int treeIdx, unsigned int numTaxa, unsigned int numTrees, unsigned long long hv1, unsigned long long hv2, vector<BitSet*> & vec_bs);
//    void hashing_bs_SC(unsigned int treeIdx, unsigned int numTaxa, unsigned int numTrees, unsigned long long hv1, unsigned long long hv2, vector<BitSet*> & vec_bs);
//    void hashing_bs_SC(BitSet& bs, unsigned int numTaxa, unsigned long long hv1, unsigned long long hv2);
    void hashing_bs_SC(unsigned int treeIdx, unsigned int numTaxa, unsigned int numTrees, unsigned long long hv1, unsigned long long hv2, vector<bool *> & vec_bs);
    void hashing_bs_SC(bool * bs, unsigned int numTaxa, unsigned long long hv1, unsigned long long hv2);

    void UHashfunc_init(unsigned int t, unsigned int n, unsigned int c);
    void UHashfunc_init(unsigned int t, unsigned int n, unsigned int c, int32 newseed);
    void HashMap_clear();
};


#endif


