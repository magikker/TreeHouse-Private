#ifndef _SEARCH_FUNCTIONS_H
#define _SEARCH_FUNCTIONS_H

#include <set>
#include <vector>
#include <string>
#include <map>
#include <random>

#include "UtilityFunctions.h"
#include "THGlobals.h"
//#include "global.h"
#include "newick.h"

using namespace std;

set<unsigned int> clade_size_search(vector<int> required, int size);

set<unsigned int> clade_size_search(vector<string> RequiredTaxa, int size);

set<unsigned int> get_trees_with_taxa(vector<int> required);

set<unsigned int> get_trees_without_taxa(vector<int> excluded);

set<unsigned int> get_trees_without_taxa(vector<string> ExcludedTaxa);

set<unsigned int> get_trees_with_taxa(vector<string> RequiredTaxa);

set<unsigned int> get_trees_by_taxa(vector<string> RequiredTaxa, vector<string> ExcludedTaxa);

set<unsigned int> search_hashtable_ktets(vector < vector < int > > subtrees);

set<unsigned int> get_trees_by_subtree(string subtree);

set<unsigned int> search_hashtable_strict(vector<int> leftside, vector<int> rightside, int side);

set<unsigned int> search_hashtable_strict_and_timed(vector<int> leftside, vector<int> rightside, int side);

set<unsigned int> search_hashtable_strict_old(vector<int> leftside, vector<int> rightside, int side);

bool * dfs_compute_bitstrings(NEWICKNODE* startNode, NEWICKNODE* parent, vector< vector < int > > &solution );

int random_search(int left, int right, int side, int iterations);

int random_search2(int left, int right, int side, int iterations);

#endif
