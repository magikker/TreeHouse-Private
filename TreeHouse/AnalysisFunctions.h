#ifndef _ANALYSIS_FUNCTIONS_HH_
#define _ANALYSIS_FUNCTIONS_HH_

#include <vector>
#include <stack>
#include <iomanip> //Formatting output
#include <math.h>


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

void print_summary_stats();

void print_window_stats();

std::vector<unsigned int> biparts_in_tree(unsigned int);

std::vector<string> distinguishing_taxa(set<unsigned int> inputtrees1, set<unsigned int> inputtrees2);

std::vector<int> distinguishing_bipart(set<unsigned int> inputtrees1, set<unsigned int> inputtrees2);

void test_trait_correlation(int t1ind, int t1val, int t2ind, int t2val, unsigned int tree, int iterations, string folder);

string psupport(vector < set < unsigned int > > treesets);

void dTree(vector <set < unsigned int >> treesets);

set <unsigned int> get_outgroup_ids(unsigned int tree);

#endif
