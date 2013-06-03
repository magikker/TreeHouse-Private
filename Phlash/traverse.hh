#ifndef TRAVERSE_HH_
#define TRAVERSE_HH_

#include <iostream>
#include <string>
#include "parsing.h"
#include "hash.hh"

using namespace std;
extern bool ROOTED;
extern bool FULLMATRIX;
extern bool QOPT;
extern bool TRIVIAL;
extern unsigned int NUM_TREES2;
extern unsigned int TOTAL_TREES;
extern unsigned int * n_biparts;
//for weighted Phlash
struct weighted_element{
  unsigned int t_id;
  float weight;
};
extern struct weighted_element ** windex;
extern struct weighted_element ** weighted_hashtable;


void mc_procedure(unsigned long long M1, unsigned long long M2, unsigned int start, unsigned int end, FILE * fp, LabelMap &lm, bool rooted);
void lv_procedure(unsigned long long M1, unsigned long long M2, unsigned int start, unsigned int end, FILE * fp, LabelMap & lm, bool rooted);
void handle_newick_error(int err);
void calculate_unique_bipartitions(unsigned int threshold, unsigned int & unique_bipartitions);
void process_trees(string infilename, string qfilename, unsigned int & hashtable_length, unsigned int * row);
void process_trees_weighted(string infilename,  string qfilename, unsigned int distance_option, string outfilename, float * ab_array);
void by_row_distance(unsigned int hash_length, unsigned int * row, unsigned int distance_measure, string outfilename);
void weighted_row_distance(float * ab_array, unsigned int distance_measure, string outfilename); 
#endif
