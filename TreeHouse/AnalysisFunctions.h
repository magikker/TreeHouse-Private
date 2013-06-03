#ifndef _ANALYSIS_FUNCTIONS_HH_
#define _ANALYSIS_FUNCTIONS_HH_

#include <vector>

//Global Vars like the hashtable.
#include "THGlobals.h"
#include "UtilityFunctions.h"
#include "global.h"

//The label map class. 
#include "label-map.hh"

using namespace std;

void generate_random_bt();

std::vector<string> distinguishing_taxa(set<unsigned int> inputtrees1, set<unsigned int> inputtrees2);

std::vector<int> distinguishing_bipart(set<unsigned int> inputtrees1, set<unsigned int> inputtrees2);

void bitpartitions_by_frequency(set<unsigned int> inputtrees, float threshold, vector< bool * > &consensus_bs, vector< float > &consensus_branchs, vector< unsigned int> &consensus_bs_sizes);

unsigned int compute_threshold(unsigned int numberofTrees, float threshold);

string consen(set<unsigned int> inputtrees, float percent);
  
BipartitionTable least_conflict_bt(set<unsigned int> inputtrees);

string least_conflict(set<unsigned int> inputtrees);

BipartitionTable get_consen_bt(set<unsigned int> inputtrees, float percent);

BipartitionTable greedy_consen_bt(set<unsigned int> inputtrees, float percent);

string greedy_consen(set<unsigned int> inputtrees, float percent);

void test_trait_correlation(int t1ind, int t1val, int t2ind, int t2val, unsigned int tree, int iterations, string folder);

#endif
