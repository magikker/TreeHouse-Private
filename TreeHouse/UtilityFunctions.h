#ifndef _UTILITY_FUNCTIONS_H
#define _UTILITY_FUNCTIONS_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include <math.h>

//#include "global.h"
#include "THGlobals.h"
#include "pql.h"
#include "AnalysisFunctions.h"
using namespace std;

//I need some work on this one. 
set <unsigned int> duplicates(int treein);

vector<string> to_newick(vector<int> input_from_int);

vector<string> to_newick(set<unsigned int> inval);

string to_newick(int x);

void show_newick(string nwstr, string title, string cssfile, string mode);

void show_newick(vector<string> nwvect, string mode);

void export_newick(string nwstr, string title, string path, string cssfile, string mode);

void export_newick(vector<string> nwvect, string path, string mode);

void show(vector< int> input_from_int, string mode);

void show(set<unsigned int> inval, string mode);

void show(int x, string mode);

void exportt(vector< int> input_from_int, string path, string mode);

void exportt(set<unsigned int> inval, string path, string mode);

void exportt(int x, string path, string mode);

void classification(vector<string> taxa);

void classification(string taxon);

void write_classifications(string filename);

void edit_classification(vector<string> taxa, string rank, string info);

void edit_classification(string taxon, string rank, string info);

bool compare_padded_taxon_names(Taxon* taxon1, Taxon* taxon2);

void load_trait_matrix(string file);

vector<string> taxa_in_group(string group);

vector<string> groups_at_level(string rank);

vector<string> parse_newick(string newick);

bool is_elem_of_vector(string item, vector<string> vect);

string get_clade_of_group(string group, string newick);

void highlight_group_css(string group, vector<string> taxa, string newick, int treeid, map<string, vector <int> > &out_of_place);

int highlight_group_taxa( string group, vector<string> taxa, int tree, string mode, map<string, vector <int> > &out_of_place, bool isolate=false );

int highlight_group(string group, set<unsigned int> treeset, string mode, bool isolate=false);
		      
int highlight_group(string group, vector<int> treevect, string mode, bool isolate=false);

int show_level(string rank, int tree, string mode);

int show_level(string rank, set<unsigned int> treeset, string mode);

int show_level(string rank, vector<int> treevect, string mode);

int show_only_taxa(string group, vector<string> taxa, int tree, string mode);

int show_only(string group, set<unsigned int> treeset, string mode);

int show_only(string group, vector<int> treevect, string mode);

int show_only(vector<string> groupvect, set<unsigned int> treeset, string mode);

int show_only(vector<string> groupvect, vector<int> treevect, string mode);

void taxa_filter (vector<string> taxa, unsigned int tree);

void taxa_filter (vector<string> taxa, set<unsigned int> treeset);

void taxa_filter(vector<string> taxa, vector<int> treevect);

void group_filter(vector<string> groups, unsigned int tree);

void group_filter (vector<string> groups, set<unsigned int> treeset);

void group_filter(vector<string> groups, vector<int> treevect);

void delete_tree(unsigned int tree);

void delete_tree(set<unsigned int> treeset);

void delete_tree(vector<int> treevect);

void write_trz(vector<string> nwvect, string filename);

void write_trz(int tree, string filename);

void write_trz(vector<int> treevect, string filename);

void write_trz(set<unsigned int> treeset, string filename);

int num_labels_in_newick(string newick);

set <unsigned int> unique(set<unsigned int> treesin);

vector <unsigned int> unique(vector<unsigned int> treesin);

vector <int> unique(vector<int> treesin);

int unique_biparts(set< unsigned int > treesin);

vector<string> split(const char *str, char c );

int pANTLR3_COMMON_TOKEN_to_int(pANTLR3_COMMON_TOKEN tok);

string pANTLR3_COMMON_TOKEN_string_lit_to_string(pANTLR3_COMMON_TOKEN tok);

string pANTLR3_COMMON_TOKEN_to_string(pANTLR3_COMMON_TOKEN tok);

double pANTLR3_COMMON_TOKEN_to_double(pANTLR3_COMMON_TOKEN tok);

vector<int> pANTLR3_COMMON_TOKEN_to_intvect(pANTLR3_COMMON_TOKEN tok);

void printBipartTier(int tier);

void printTaxaTierTuples(int tier, int taxa);

void rateTierRogueness(int tier);

void printBipartition(vector<string> leftside, vector<string> rightside);

void printBipartition(int bipartID);

void printSetOfTrees(set<unsigned int> trees);

void printBipartition(vector<unsigned int> leftside, vector<unsigned int> rightside);

void print_bitstring(bool * bitstring, unsigned int length);

void print_vector_of_bs(vector< bool * > bitstrings, unsigned int length_of_bitstrings);

template <class T> //this might not work- needs to be tested
void printVector(vector<T> in){
  cout << endl;
  for(int i = 0; i < in.size(); i++){
	cout << i+1 << ". " << in.at(i) << endl;
	}
}

template <class T> //this might not work- needs to be tested
void printSet(set<T> in){
  for(int i = 0; i < in.size(); i++){
	cout << i+1 << ". " << in.at(i) << endl;
	}
  cout << endl;
}

template <class T> 
void printVectorCompact(vector<T> in){
  cout << endl;
  for(int i = 0; i < in.size(); i++){
	cout << in.at(i) << ", ";}

}
/*
template <class T> 
void printSetCompact(set<T> in){
  cout << endl;
  for(set<typename T>::iterator it = in.begin(); it!=in.end(); it++){
	cout << *it << ", ";}

}*/

void printSetCompactTwo(set<unsigned int> in);

template <class T>
vector<T> copyToVector(T* in, unsigned int size){
  vector<T> returnVec;
  for(int i = 0; i < size; i++){
	returnVec.push_back(in[i]);
	} 
}

void print_list_bs(vector< bool * > list_bs);

void print_hashtable();

void print_set(set<unsigned int> t);

void print_vector_of_strings(vector< string > bitstrings);

//use find instead. 
//bool isInVector(vector<int> toSearch, int x);

unsigned int factorial(int);
#endif
