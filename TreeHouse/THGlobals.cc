#include <vector>
#include <iostream>
#include <set>
#include <algorithm>

#include <stdlib.h>
#include <string.h>

#include "label-map.hh"
#include "BipartitionTable.h"
#include "THGlobals.h"
#include "global.h"

#include "pqlsymbol.h"
#include "newick.h"

using namespace std;

//vector<int> shuffledTaxa;
bool DEBUGMODE;

//bool hetero; // = 0; // If the trees are taxa heterogenious
//unsigned num_trees; // = 0; // number of trees
//unsigned num_taxa; //  = 0; // number of taxa

BipartitionTable biparttable;

//TreeTable treetable;

voidFuncts voidFunctMap;
argFuncts argFunctMap;
std::vector<std::string> functionKeys;

vector< vector<unsigned int> > inverted_index;	
//vector< bool* > taxa_in_trees;  // which taxa are in which trees					NEED

//LabelMap lm;

vector<Taxon *> taxa_info;

set< unsigned int > all_trees;
set< unsigned int > original_trees;
int NUM_TREES_INIT;

std::map<int, vector<int> > tree_dups;
 
vector< pqlsymbol * > query_results;			// the most recent 

std::map<string, pqlsymbol * > symbol_table;
std::map<string, bool > constant_table;

ofstream output_file;

//timing mechanisms
vector<std::clock_t> clocks;

	
unsigned int SetInsertions;
unsigned int SetOps;
double SetTime;
double HetTime;
double AHetTime;
double SearchTime;

bool write_to_output(string input){
  if (output_file.is_open())
  {
    output_file << input;
    return true;
  }
  else cout << "Unable to open file";
  return false;
}

bool change_output_pointer(string input){
	if (output_file.is_open()){
		output_file.close();
	}
	
	output_file.open(input);
	if (output_file.is_open()){
		return true;
	}
	
  else cout << "Unable to open file";
  return false;
}

bool init_output(){
  output_file.open("logs/output.txt", fstream::app);
  
  if (output_file.is_open()){
    output_file << "This is a line.\n";
    output_file << "This is another line.\n";
    return true;
  }
  else cout << "Unable to open file logs/interactive.txt";
  return false;
}
/*
bool is_taxa_homogenious(set<unsigned int> treeset){
	if(!::HETERO){
		return true;
	}
	else{
		for (unsigned int i = 0; i < treeset.size(); i++){
			for (unsigned int j = 0; j < ::NUM_TAXA; j++){
				if (taxa_in_trees[0][j] != taxa_in_trees[i][j]){
					return false;
				}
			}
		}
	}
	return true;
}
*/
/*
bool are_taxa_in_tree(int treeindex, vector<int> setoftaxa){
	cout << "first taxon is: " << setoftaxa[0] << endl;

	for(unsigned int i = 0; i < setoftaxa.size(); i++){
		if(taxa_in_trees[treeindex][setoftaxa[i]] == 0){ 
	//	NOTE- ABOVE LINE USED TO BE:	if(taxa_in_trees[treeindex][i] == 0){ 
			return false;
		}
	}
	return true;
}
*/
/*
int num_taxa_in_tree(int treeindex){
	int count = 0;
	for(unsigned int i = 0; i < ::NUM_TAXA; i ++ ){
		if (taxa_in_trees[treeindex][i] == 1){
			count += 1;
		}
	}
	return count;
}
*/
/*
vector<string> get_taxa_in_tree(unsigned int treeindex){
	set<string> taxaset;
	for(unsigned int i = 0; i < ::NUM_TAXA; i ++ ){
		if (taxa_in_trees[treeindex][i] == 1){
			taxaset.insert(lm.name(i));
		}
	}
	vector<string> taxavect(taxaset.begin(), taxaset.end());
	return taxavect;
}
*/
/*
vector<string> get_all_taxa_vect(){
	vector<string> taxavect;
	for(unsigned int i = 0; i < :lm.size(); i ++ ){
			taxavect.push_back(lm.name(i));
	}
	return taxavect;
}

set<string> get_all_taxa_set(){
	set<string> taxaset;
	for(unsigned int i = 0; i < lm.size(); i ++ ){
			taxaset.insert(lm.name(i));
	}
	return taxaset;
}
*/
/*
void print_taxa_in_trees(){
	cout << "taxa_in_trees" << endl;
	for(unsigned int i = 0; i < taxa_in_trees.size(); i++) {
		for (unsigned int j = 0; j < ::NUM_TAXA; j++){
			cout << taxa_in_trees[i][j]; 
		}
		cout << endl; 
	} 
}
*/
vector<int> get_taxa_with_trait(unsigned int trait_id, int trait_value=1) {
  vector<int> taxa_vect;
  for (unsigned int i=0; i < ::taxa_info.size(); i++) {
    if (::taxa_info[i]->traits[trait_id] == trait_value)
      taxa_vect.push_back(i);
  }
  return taxa_vect;
}

vector<int> get_taxa_without_trait(unsigned int trait_id, int trait_value=1) {
  vector<int> taxa_vect;
  for (unsigned int i=0; i < ::taxa_info.size(); i++) {
    if (::taxa_info[i]->traits[trait_id] != trait_value)
      taxa_vect.push_back(i);
  }
  return taxa_vect;
}
/*
int distance_between_taxa(unsigned int taxon1, unsigned int taxon2, unsigned int tree) {
  vector< bool *> tree_bipartitions = get_tree_bipartitions(tree);
  vector<unsigned int> tree_bs_sizes = get_tree_bs_sizes(tree);
  int distance = 0;
  for (unsigned int i=0; i<tree_bipartitions.size(); i++) {
    if ((taxon1 < tree_bs_sizes[i] && tree_bipartitions[i][taxon1]) != (taxon2 < tree_bs_sizes[i] && tree_bipartitions[i][taxon2]))
      distance++;
  }
  return distance;
}
*/
//editted needs retesting
int distance_between_taxa(unsigned int taxon1, unsigned int taxon2, unsigned int tree) {
  vector< Bipartition > tree_bipartitions = get_tree_bipartitions(tree);
  vector<unsigned int> tree_bs_sizes = get_tree_bs_sizes(tree);
  int distance = 0;
  for (unsigned int i=0; i<tree_bipartitions.size(); i++) {
    if (!tree_bipartitions[i].same_bitstring_value(taxon1, taxon2))
      distance++;
  }
  return distance;
}

int distance_to_common_ancestor(unsigned int taxon1, unsigned int taxon2, unsigned int tree) {
  string taxon1name = ::biparttable.lm.name(taxon1);
  string taxon2name = ::biparttable.lm.name(taxon2);
  string taxastring = taxon1name + " " + taxon2name;
  system(("echo \""+to_newick(tree)+"\" | ../NwUtils/src/nw_clade - "+taxastring+" > temp/clade.txt").c_str());
  ifstream cladefile ("temp/clade.txt");
  string clade;
  getline (cladefile, clade);
  cladefile.close();
  vector<string> clade_parsed = parse_newick(clade);
  int depth = 0;
  int i = 0;
  //count depth from left side of string
  for (i; i<clade_parsed.size(); i++) {
    if (clade_parsed[i] == "(")
      depth++;
    else if (clade_parsed[i] == ")")
      depth--;
    if (clade_parsed[i] == taxon1name)
      return depth;
    else if (clade_parsed[i] == taxon2name)
      break;
  }
  //taxon2 found before taxon1. Locate taxon1.
  for (i; i<clade_parsed.size(); i++) {
    if (clade_parsed[i] == taxon1name)
      break;
  }
  depth = 0;
  //count depth to right side of string
  for (i; i<clade_parsed.size(); i++) {
    if (clade_parsed[i] == "(")
      depth--;
    else if (clade_parsed[i] == ")")
      depth++;
  }
  return depth;
}

vector<int> get_taxa_in_clade(vector<int> taxa, unsigned int tree) {
  //vector<int> taxa_in_clade;
  //vector< bool *> tree_bipartitions = get_tree_bipartitions(tree);
  //vector<unsigned int> tree_bs_sizes = get_tree_bs_sizes(tree);
  //vector< bool *> clade_bipartitions;
  //vector< bool *> clade_bs_sizes;
  //for (unsigned int i=0; i<tree_bipartitions.size(); i++) {
  //  if ((taxon1 < tree_bs_sizes[i] && tree_bipartitions[i][taxon1]) != (taxon2 < tree_bs_sizes[i] && tree_bipartitions[i][taxon2])) {
  //    clade_bipartitions.push_back(tree_bipartitions[i]);
  //    clade_bs_sizes.push_back(tree_bs_sizes[i])
  //}
  string taxastring;
  for (unsigned int i = 0; i < taxa.size(); i++) {
    taxastring += (" " + ::biparttable.lm.name(taxa[i]));
  }
  system(("echo \""+to_newick(tree)+"\" | ../NwUtils/src/nw_clade -"+taxastring+" > temp/clade.txt").c_str());
  ifstream cladefile ("temp/clade.txt");
  string clade;
  getline (cladefile, clade);
  cladefile.close();
  vector<string> clade_labels = parse_newick(clade);
  for (unsigned int i=0; i<clade_labels.size(); i++) {
    if (clade_labels[i] == "(" || clade_labels[i] == ")") {
      clade_labels.erase(clade_labels.begin() + i);
      i--;
    }
  } 
  vector<int> clade_indices = ::biparttable.lm.lookUpLabels(clade_labels);
  //sort indices lowest-to-highest
  std::sort(clade_indices.begin(), clade_indices.end());
  return clade_indices;
}

//we ought to be importing this stuff. 
/*
string to_lower(string str) {
  for (unsigned int i = 0; i < str.size(); i++)
    str[i] = tolower(str[i]);
  return str;
}
*/

vector< Bipartition > get_tree_bipartitions(unsigned int id) {
  vector < Bipartition > tree_bipartitions;
  for (unsigned int i = 0; i < ::biparttable.BipartitionTable.size(); i++) {
    if (::biparttable.treetable[i][id])
      tree_bipartitions.push_back(::biparttable.BipartitionTable[i]);
  }
  return tree_bipartitions;
}

/*
vector< bool *> get_tree_bipartitions(unsigned int id) {
  vector < bool *> tree_bipartitions;
  for (unsigned int i = 0; i < ::biparttable.bipartitions.size(); i++) {
    if (::biparttable.treetable[i][id])
      tree_bipartitions.push_back(::biparttable.bipartitions[i]);
  }
  return tree_bipartitions;
}
*/

vector<unsigned int> get_tree_bs_sizes(unsigned int id) {
  vector <unsigned int> tree_bs_sizes;
  for (unsigned int i = 0; i < ::biparttable.BipartitionTable.size(); i++) {
    if (::biparttable.treetable[i][id])
      tree_bs_sizes.push_back(::biparttable.bitstring_size(i));
  }
  return tree_bs_sizes;
}

/*
vector<unsigned int> get_tree_bs_sizes(unsigned int id) {
  vector <unsigned int> tree_bs_sizes;
  for (unsigned int i = 0; i < ::biparttable.bipartitions.size(); i++) {
    if (::biparttable.treetable[i][id])
      tree_bs_sizes.push_back(::biparttable.length_of_bitstrings[i]);
  }
  return tree_bs_sizes;
}
*/

vector<float> get_tree_branches(unsigned int id) {
  return ::biparttable.tree_branches[id];
}

//This needs to be remaned. Way to general. 
vector<unsigned int> get_tree_data(unsigned int tree_id, vector<bool *>& tree_bipartitions, vector <unsigned int>& tree_bs_sizes, vector<float>& tree_branches) {
	//cout<< "welcome to get_tree_data" <<endl;
	vector <unsigned int> bipart_indices; //which bipartitions we're returning
	for (unsigned int i = 0; i < ::biparttable.BipartitionTable.size(); i++) {
		if (::biparttable.treetable[i][tree_id]) { //if the bipartitin is in the tree_id
			tree_bipartitions.push_back(::biparttable.get_boolarray(i));
			tree_bs_sizes.push_back(::biparttable.bitstring_size(i));
		}
	}
	tree_branches = ::biparttable.tree_branches[tree_id];
	//cout << "returning tree data" <<endl;
	return bipart_indices;
}
/*vector<unsigned int> get_tree_data(unsigned int id, vector < bool *>& tree_bipartitions, vector <unsigned int>& tree_bs_sizes, vector<float>& tree_branches) {
  vector <unsigned int> bipart_indices;
  for (unsigned int i = 0; i < ::biparttable.bipartitions.size(); i++) {
    if (::biparttable.treetable[i][id]) {
      tree_bipartitions.push_back(::biparttable.bipartitions[i]);
      bipart_indices.push_back(i);
      tree_bs_sizes.push_back(::biparttable.length_of_bitstrings[i]);
    }
  }
  tree_branches = ::biparttable.tree_branches[id];
  return bipart_indices;
}*/

string th_compute_tree(BipartitionTable& bpt, unsigned id, bool branch) {
	//cout << "welcome to th_compute_tree in THGLOBALS" <<endl;
	vector< bool * > my_bs;
	vector< float > my_branches;
	vector< unsigned int> bs_sizes;
	get_tree_data(id, my_bs, bs_sizes, my_branches);
	//cout << "got the tree data from get_tree_data" <<endl;
	
	//cout << "my_branches.size() = "<< my_branches.size() << endl;
	
	//for (unsigned int i = 0; i < bs_sizes.size(); i++){
	//	cout << bs_sizes[i] << " : ";
	//	for (unsigned int j = 0; j < bs_sizes[i]; j++){
	//		cout << my_bs[i][j];
	//	}
	//	
	//	cout << " : " << my_branches[i];
	//	
	//	cout << endl;
	//}
	return compute_tree(::biparttable.lm, my_bs, my_branches, id, branch, bs_sizes);
	
}

//delete
void copy_labelmap(LabelMap& lm1, LabelMap& lm2) {
  lm1.clear();
  for (unsigned int i = 0; i < lm1.size(); i++)
    lm2.push(lm1.name(i));
}

//not used
/*void recompute_tree_dups() {
  for (unsigned int i = 0; i < ::NUM_TREES; i++) {
    tree_dups[i].clear();
    for (unsigned int j = 0; j < ::NUM_TREES; j++) {
      bool is_dup = true;
      if (i != j) {
	for (unsigned int k = 0; k < ::biparttable.bipartitions.size(); k++) {
	  is_dup = !(biparttable.treetable[k][i] ^ biparttable.treetable[k][j]);
	  if (!is_dup)
	    break;
	}
	if (is_dup)
	  tree_dups[i].push_back(j);
      }
    }
  }
}
*/


/*
vector<unsigned int> which_trees_double_strict(int bitstringindex, int numTaxaSearched){
	vector<unsigned int> v;
	for(unsigned int i = 0; i < biparttable.searchtable[bitstringindex].size(); i++){ //for each tree in with the biparitition we need to look at which taxa it has... and check that the input taxa are the only taxa marked as 0's in the bipartition. 
		if (num_taxa_in_tree(biparttable.searchtable[bitstringindex][i]) == numTaxaSearched){
			v.push_back(biparttable.searchtable[bitstringindex][i]);
		}
	}
	return v;
}
*/

//new and needs testing
vector<unsigned int> which_trees_double_strict(int bitstringindex, int numTaxaSearched){
	vector<unsigned int> v;
	for(vector<unsigned int >::iterator it = biparttable.BipartitionTable[bitstringindex].trees_begin(); it != biparttable.BipartitionTable[bitstringindex].trees_end(); it++){
		if (biparttable.num_taxa_in_tree(*it) == numTaxaSearched){
			v.push_back(*it);
		}
 	}
	return v;
}


/*
vector<unsigned int> which_trees_single_strict(int bitstringindex, vector<int> positions){
	unsigned int mycount = biparttable.number_of_ones(bitstringindex);
	vector<unsigned int> v;
	if (biparttable.is_one(bitstringindex,positions[0])){
		if(positions.size() == mycount){
			return biparttable.get_trees(bitstringindex); //all
		}
		else{
			return v; // or nothing
		}
	}
	else{
		for(unsigned int i = 0; i < biparttable.searchtable[bitstringindex].size(); i++){ //for each tree in with the biparitition we need to look at which taxa it has... and check that the input taxa are the only taxa marked as 0's in the bipartition. 
			if (num_taxa_in_tree(biparttable.searchtable[bitstringindex][i])- mycount == positions.size()){
				v.push_back(biparttable.searchtable[bitstringindex][i]); 
			}
		}
	}
	return v;
}
*/
vector<unsigned int> which_trees_single_strict(int bitstringindex, vector<int> positions){
	unsigned int mycount = biparttable.number_of_ones(bitstringindex);
	vector<unsigned int> v;
	if (biparttable.is_one(bitstringindex,positions[0])){
		if(positions.size() == mycount){
			return biparttable.get_trees(bitstringindex); //all
		}
		else{
			return v; // or nothing
		}
	}
	else{
		for(vector<unsigned int >::iterator it = biparttable.BipartitionTable[bitstringindex].trees_begin(); it != biparttable.BipartitionTable[bitstringindex].trees_end(); it++){
			if (biparttable.num_taxa_in_tree(*it)- mycount == positions.size()){
				v.push_back(*it); 
			}
		}
	}
	return v;
}

//not sure if it works in all cases. not accounting for the zeros needing to be in the tree. 
/*
vector<unsigned int> which_strict_hetero(int bitstringindex, vector<int> positions1, vector<int> positions2){
	//this function returns the subset of trees found at searchtable[bitstringindex] 
	//where taxa in vector<int> positions1 and vector<int> positions2
	//are the only taxa in the trees 
	vector<unsigned int> v;	
	for(unsigned int i = 0; i < biparttable.searchtable[bitstringindex].size(); i++){ //for each tree in with the biparitition we need to look at which taxa it has... and check that the input taxa are the only taxa marked as 0's in the bipartition. 
		if (num_taxa_in_tree(biparttable.searchtable[bitstringindex][i]) == (positions1.size() + positions2.size()) ){
			v.push_back(biparttable.searchtable[bitstringindex][i]);
		}
	}
	return v;
}
*/
vector<unsigned int> which_strict_hetero(int bitstringindex, vector<int> positions1, vector<int> positions2){
	//this function returns the subset of trees found at searchtable[bitstringindex] 
	//where taxa in vector<int> positions1 and vector<int> positions2
	//are the only taxa in the trees 
	vector<unsigned int> v;	
	for(vector<unsigned int >::iterator it = biparttable.BipartitionTable[bitstringindex].trees_begin(); it != biparttable.BipartitionTable[bitstringindex].trees_end(); it++){
		if (biparttable.num_taxa_in_tree(*it) == (positions1.size() + positions2.size()) ){
			v.push_back(*it);
		}
	}
	return v;
}
/*
vector<unsigned int> which_strict_hetero(int bitstringindex, vector<int> positions){
	unsigned int mycount = biparttable.number_of_ones(bitstringindex);
	vector<unsigned int> v;
	if (biparttable.is_one(bitstringindex,positions[0])){
		if(positions.size() == mycount){
			return biparttable.searchtable[bitstringindex];
		}
		else{
			return v;
		}
	}
	else{
		for(unsigned int i = 0; i < biparttable.searchtable[bitstringindex].size(); i++){ //for each tree in with the biparitition we need to look at which taxa it has... and check that the input taxa are the only taxa marked as 0's in the bipartition. 
			if (num_taxa_in_tree(biparttable.searchtable[bitstringindex][i])-mycount == positions.size()){
				v.push_back(biparttable.searchtable[bitstringindex][i]);
			}
		}
	}
	return v;
}
*/
vector<unsigned int> which_strict_hetero(int bitstringindex, vector<int> positions){
	unsigned int mycount = biparttable.number_of_ones(bitstringindex);
	vector<unsigned int> v;
	if (biparttable.is_one(bitstringindex,positions[0])){
		if(positions.size() == mycount){
			return biparttable.get_trees(bitstringindex);
		}
		else{
			return v;
		}
	}
	else{
		for(vector<unsigned int >::iterator it = biparttable.BipartitionTable[bitstringindex].trees_begin(); it != biparttable.BipartitionTable[bitstringindex].trees_end(); it++){
			if (biparttable.num_taxa_in_tree(*it)-mycount == positions.size()){
				v.push_back(*it);
			}
		}
	}
	return v;
}


// This tells me which trees linked to the bitstringindex have all the passed in taxa. 
// PHASE OUT
/*
vector<unsigned int> which_hetero(int bitstringindex, vector<int> positions){
	vector<unsigned int> v;
	if (biparttable.is_one(bitstringindex,positions[0])){
		return biparttable.searchtable[bitstringindex];
	}
	else{
		for(unsigned int i = 0; i < biparttable.searchtable[bitstringindex].size(); i++){ //for each tree in with the biparitition we need to look at which taxa it has... and check taxa are there. 
			if ( are_taxa_in_tree(biparttable.searchtable[bitstringindex][i], positions) ){
				v.push_back(biparttable.searchtable[bitstringindex][i]);
			}
		}
	}
	return v;
}
*/
vector<unsigned int> which_hetero(int bitstringindex, vector<int> positions){
	vector<unsigned int> v;
	if (biparttable.is_one(bitstringindex,positions[0])){
		return biparttable.get_trees(bitstringindex);
	}
	else{
		for(vector<unsigned int >::iterator it = biparttable.BipartitionTable[bitstringindex].trees_begin(); it != biparttable.BipartitionTable[bitstringindex].trees_end(); it++){
			if ( biparttable.are_taxa_in_tree(*it, positions) ){
				v.push_back(*it);
			}
		}
	}
	return v;
}

bool is_strict_homog(int bitstringindex, vector<int> positions){
	unsigned int mycount = biparttable.number_of_ones(bitstringindex);
		
	if (biparttable.is_one(bitstringindex,positions[0])){
		if(positions.size() == mycount){
			return true;
		}
	}
		
	else{
		if(positions.size() == ::NUM_TAXA-mycount){
			return true;
		}
	}
	return false;		
}

//bool is_compat_ind(int bitstringindex1, int bitstringindex2) {
// bool x1ny1 = true;
//bool x1ny2 = true;
//bool x2ny1 = true;
//bool x2ny2 = true;
//bool fullbitstring1 [NUM_TAXA];
//bool fullbitstring2 [NUM_TAXA];
//for (unsigned int i = 0; i < NUM_TAXA; i++) {
//  if (!biparttable.in_bitstring(bitstringindex1, i) || biparttable.is_zero(bitstringindex1, i))
//    fullbitstring1[i] = 0;
//  else
//    fullbitstring1[i] = 1;
//  if (!biparttable.in_bitstring(bitstringindex2, i) || biparttable.is_zero(bitstringindex2, i))
//    fullbitstring2[i] = 0;
//  else
//    fullbitstring2[i] = 1;
//}
//for (unsigned int i = 0; i < NUM_TAXA; i++) {
//  if (!x1ny1 && !x1ny2 && !x2ny1 && !x2ny2)
//    break;
//  else if (fullbitstring1[i] == 1 && fullbitstring2[i] == 1)
//    x1ny1 = false;
//  else if (fullbitstring1[i] == 1 && fullbitstring2[i] == 0)
//    x1ny2 = false;
//  else if (fullbitstring1[i] == 0 && fullbitstring2[i] == 1)
//    x2ny1 = false;
//  else if (fullbitstring1[i] == 0 && fullbitstring2[i] == 0)
//    x2ny2 = false;
//}
//return (x1ny1 || x1ny2 || x2ny1 || x2ny2);
//}
/*
bool is_compat(bool *bitstring1, int length1, bool *bitstring2, int length2) {
  bool x1ny1 = true;
  bool x1ny2 = true;
  bool x2ny1 = true;
  bool x2ny2 = true;
  bool fullbitstring1 [NUM_TAXA];
  bool fullbitstring2 [NUM_TAXA];
  for (unsigned int i = 0; i < NUM_TAXA; i++) {
    if (i >= length1 || bitstring1[i] == 0)
      fullbitstring1[i] = 0;
    else
      fullbitstring1[i] = 1;
    if (i >= length2 || bitstring2[i] == 0)
      fullbitstring2[i] = 0;
    else
      fullbitstring2[i] = 1;
  }
  for (unsigned int i = 0; i < NUM_TAXA; i++) {
    if (!x1ny1 && !x1ny2 && !x2ny1 && !x2ny2)
      break;
    else if (fullbitstring1[i] == 1 && fullbitstring2[i] == 1)
      x1ny1 = false;
    else if (fullbitstring1[i] == 1 && fullbitstring2[i] == 0)
      x1ny2 = false;
    else if (fullbitstring1[i] == 0 && fullbitstring2[i] == 1)
      x2ny1 = false;
    else if (fullbitstring1[i] == 0 && fullbitstring2[i] == 0)
      x2ny2 = false;
  }
  return (x1ny1 || x1ny2 || x2ny1 || x2ny2);
}
*/
//for bipartitions to be compatible one has to be a subset of the other. 
bool is_compat(boost::dynamic_bitset<>  bitstring1, boost::dynamic_bitset<>  bitstring2) {
  bool x1ny1 = true;
  bool x1ny2 = true;
  bool x2ny1 = true;
  bool x2ny2 = true;
  bool fullbitstring1 [NUM_TAXA];
  bool fullbitstring2 [NUM_TAXA];
  //GRB need a get non-truncated_bitstring function? 
  for (unsigned int i = 0; i < NUM_TAXA; i++) {
    if (i >= bitstring1.size() || bitstring1[i] == 0)
      fullbitstring1[i] = 0;
    else
      fullbitstring1[i] = 1;
    if (i >= bitstring2.size() || bitstring2[i] == 0)
      fullbitstring2[i] = 0;
    else
      fullbitstring2[i] = 1;
  }
  //checking if one is the subset of the other. 
  for (unsigned int i = 0; i < NUM_TAXA; i++) {
    if (!x1ny1 && !x1ny2 && !x2ny1 && !x2ny2)
      break;
    else if (fullbitstring1[i] == 1 && fullbitstring2[i] == 1)
      x1ny1 = false;
    else if (fullbitstring1[i] == 1 && fullbitstring2[i] == 0)
      x1ny2 = false;
    else if (fullbitstring1[i] == 0 && fullbitstring2[i] == 1)
      x2ny1 = false;
    else if (fullbitstring1[i] == 0 && fullbitstring2[i] == 0)
      x2ny2 = false;
  }
  return (x1ny1 || x1ny2 || x2ny1 || x2ny2);
}

bool is_compat(int bitstringindex1, int bitstringindex2) {
  return is_compat(biparttable.get_bitstring(bitstringindex1), biparttable.get_bitstring(bitstringindex2));
    }

bool is_compat(boost::dynamic_bitset<>  bitstring1, int bitstringindex2) {
  return is_compat(bitstring1, biparttable.get_bitstring(bitstringindex2)); 
    }
/*
vector<unsigned int> find_incompat_old(bool *bitstring, int length, set<unsigned int> inputtrees) {
  vector<unsigned int> incompat_indices;
  unsigned int i = 0;
 top:
  for (i; i < biparttable.bipartitions.size(); i++) {
    for (unsigned int j = 0; j < biparttable.searchtable[i].size(); j++) {
      for (set<unsigned int>::const_iterator pos = inputtrees.begin(); pos != inputtrees.end(); pos++) {
	if (biparttable.searchtable[i][j] == *pos) {
	  if (!is_compat(bitstring, length, i)) {
	    incompat_indices.push_back(i);
	    i++;
	    goto top;
	  }
	  break;
	}
      }
    }
  }
  return incompat_indices;
}*/

/*
vector<unsigned int> find_incompat(bool *bitstring, int length, set<unsigned int> inputtrees) {
  vector<unsigned int> incompat_indices;
  for (unsigned int i = 0; i < biparttable.biparttable_size(); i++) {
    set<unsigned int> trees_with_bipart(biparttable.BipartitionTable[i].trees_begin(), biparttable.BipartitionTable[i].trees_end());
    set<unsigned int> isect;
    std::set_intersection(inputtrees.begin(), inputtrees.end(), trees_with_bipart.begin(), trees_with_bipart.end(), std::inserter(isect, isect.begin()));
    if (!isect.empty() && !is_compat(bitstring, length, i)) {
      incompat_indices.push_back(i);
    }
  }
  return incompat_indices;
}
* */
vector<unsigned int> find_incompat(boost::dynamic_bitset<> bitstring, set<unsigned int> inputtrees) {
  vector<unsigned int> incompat_indices;
  for (unsigned int i = 0; i < biparttable.biparttable_size(); i++) {
    set<unsigned int> trees_with_bipart(biparttable.BipartitionTable[i].trees_begin(), biparttable.BipartitionTable[i].trees_end());
    set<unsigned int> isect;
    std::set_intersection(inputtrees.begin(), inputtrees.end(), trees_with_bipart.begin(), trees_with_bipart.end(), std::inserter(isect, isect.begin()));
    if (!isect.empty() && !is_compat(bitstring, i)) {
      incompat_indices.push_back(i);
    }
  }
  return incompat_indices;
}

void debugstatement(string input){
	if(DEBUGMODE){
		cout << input << endl;
	}
}

string get_time_stamp(){
	//gets current data and time
	time_t now = time(0);
	tm* localtm = localtime(&now);
	return asctime(localtm);
}

void start_clock(){
	//starts a clock and pushes it onto the back of the clock stack
	std::clock_t start = std::clock();
	clocks.push_back(start);
}
 
double stop_clockbp(){
	//returns time since clock on back of stack. Pops that clock
	double secs = ( ( std::clock() - clocks.back() ) / (double)CLOCKS_PER_SEC );
	clocks.pop_back();
	return secs;
}

double stop_clockb(){
	//returns time since clock on back of stack.
	double secs = ( ( std::clock() - clocks.back() ) / (double)CLOCKS_PER_SEC );
	return secs;
}

double stop_clockf(){
	//returns time since clock on front of stack.
	double secs = ( ( std::clock() - clocks.front() ) / (double)CLOCKS_PER_SEC );
	return secs;
}
