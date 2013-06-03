#include <iostream>
#include <set>
#include "label-map.hh"
#include "newick.h"
#include "pql.h"

using namespace std;

#ifndef _HASHTABLESEARCH_HH_
#define _HASHTABLESEARCH_HH_

void handle_newick_error(int err);

set<unsigned int> subtree_to_ktets(string subtree, vector< bool * > list_bs);
bool * dfs_compute_bitstrings(NEWICKNODE* startNode, NEWICKNODE* parent, vector< string > &solution );

void print_vector_of_strings(vector< string > bitstrings);
void print_bitstring(bool * bitstring, unsigned int length);
void print_vector_of_bs(vector< bool * > bitstrings, unsigned int length_of_bitstrings);
void print_list_bs(vector< bool * > list_bs);
void print_hashtable();
void printBipartition(vector<string> leftside, vector<string> rightside);
void printBipartition(vector<unsigned int> leftside, vector<unsigned int> rightside);

set<unsigned int> taxaRequirements(vector<string> RequiredTaxa, vector<string> ExcludedTaxa);
//void search_hashtable(vector< bool * > list_bs, vector< bool * > &search_bs);
//set<unsigned int> search_hashtable2(vector< bool * > list_bs, vector<int> oneside, vector<int> otherside);
set<unsigned int> search_hashtable_strict(vector< bool * > list_bs, vector<unsigned int> leftside, vector<unsigned int> rightside, int side);

void LookUpLabels(vector<string> leftside, vector<string> rightside, vector<unsigned int> &lside, vector<unsigned int> &rside);
set<unsigned int> getTrees(vector<string> leftside, vector<string> rightside, vector< bool * > &list_bs, int strictFlag);

vector<unsigned int> returnAllTrees(unsigned int numberOfTaxa);

int pANTLR3_COMMON_TOKEN_to_int(pANTLR3_COMMON_TOKEN tok);

string pANTLR3_COMMON_TOKEN_string_lit_to_string(pANTLR3_COMMON_TOKEN tok);

string pANTLR3_COMMON_TOKEN_to_string(pANTLR3_COMMON_TOKEN tok);

double pANTLR3_COMMON_TOKEN_to_double(pANTLR3_COMMON_TOKEN tok);

vector<int> pANTLR3_COMMON_TOKEN_to_intvect(pANTLR3_COMMON_TOKEN tok);

#endif

