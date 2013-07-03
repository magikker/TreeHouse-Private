//file TreeParsing.h
#ifndef _TREE_PARSING_HH_
#define _TREE_PARSING_HH_

//Global Vars like the hashtable.
#include "THGlobals.h"
#include "SCNode.h"
#include "SCTree.h"
#include "Bipartition.h"
#include <fstream>
#include <string>

void populate_integrals(unsigned int * hold_integrals, string branches_int, unsigned int encode_size);
void decompress_branch(unsigned int * hold_integrals, vector<unsigned int> my_set_of_ids, Bipartition &B, string branches_frac);
unsigned int get_bitstring_length(string bitstring);
unsigned int get_ntrees(string str);
unsigned int get_unique(string str, unsigned int ntrees);
void parse_and_get(string str, string check, unsigned int & var);
unsigned int decode(string encoded, unsigned int * found);
void decode_bitstring(string bitstring, boost::dynamic_bitset<> &bs, unsigned int maxLength);
void load_data_from_trz_file(string file, BipartitionTable &Tab);

string compute_tree(
    LabelMap lm,
    vector< boost::dynamic_bitset<> > my_bs,
    vector< float > my_branches,
    unsigned id,
    bool branch);

#endif
