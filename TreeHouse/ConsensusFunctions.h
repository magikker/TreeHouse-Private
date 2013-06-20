#ifndef _CONSENSUS_FUNCTIONS_HH_
#define _CONSENSUS_FUNCTIONS_HH_

#include <vector>
#include <stack>

//Global Vars like the hashtable.
#include "THGlobals.h"
#include "UtilityFunctions.h"
//#include "global.h"

//The label map class. 
#include "label-map.hh"

#include "newick.h"
#include "SearchFunctions.h"

using namespace std;

void bitpartitions_by_frequency(set<unsigned int> inputtrees, float threshold, vector< bool * > &consensus_bs, vector< float > &consensus_branchs, vector< unsigned int> &consensus_bs_sizes);

unsigned int compute_threshold(unsigned int numberofTrees, float threshold);

string consen(set<unsigned int> inputtrees, float percent);

float consensus_reso_rate(set<unsigned int> inputtrees, float percent);

float reso_rate(string inputtree);
  
BipartitionTable least_conflict_bt(set<unsigned int> inputtrees);

string least_conflict(set<unsigned int> inputtrees);

BipartitionTable get_consen_bt(set<unsigned int> inputtrees, float percent);

BipartitionTable greedy_consen_bt(set<unsigned int> inputtrees, float percent);

string greedy_consen(set<unsigned int> inputtrees, float percent);

#endif
