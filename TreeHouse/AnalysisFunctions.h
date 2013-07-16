#ifndef _ANALYSIS_FUNCTIONS_HH_
#define _ANALYSIS_FUNCTIONS_HH_

#include <vector>
#include <stack>
#include <iomanip> //Formatting output

//Global Vars like the hashtable.
#include "THGlobals.h"
#include "UtilityFunctions.h"

#include "global.h"
#include "distance.h"

//The label map class. 
#include "label-map.hh"

#include "newick.h"
#include "SearchFunctions.h"

using namespace std;

void generate_random_bt();

std::vector<unsigned int> biparts_in_tree(unsigned int);

std::vector<string> distinguishing_taxa(set<unsigned int> inputtrees1, set<unsigned int> inputtrees2);

std::vector<int> distinguishing_bipart(set<unsigned int> inputtrees1, set<unsigned int> inputtrees2);

vector < vector < unsigned int>> compute_distances(set <unsigned int> treeset, string dist_type);

vector < vector <unsigned int>> compute_distances(vector <unsigned int> treevect, string dist_type);

void test_trait_correlation(int t1ind, int t1val, int t2ind, int t2val, unsigned int tree, int iterations, string folder);
   

#endif
