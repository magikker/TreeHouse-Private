#ifndef _BIPARTITIONTABLE_H
#define _BIPARTITIONTABLE_H

#include <map>
#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
#include "global.h"
#include "Bipartition.h"
using namespace std;

class BipartitionTable {

public:

	vector<Bipartition> BipartitionTable;
	LabelMap lm;
	
	//vector< bool* > taxa_in_trees;  // which taxa are in which trees					NEED
	vector< boost::dynamic_bitset<> > taxa_in_trees; // which taxa are in which trees					NEED
	
	//by bipartitions
	/**///vector < vector < unsigned int > > searchtable; //number of biparts long by which trees with the bipart wide
	vector < bool * > treetable; //number of biparts long by number of trees wide
	/**///vector < bool * > bipartitions; //number of biparts long by number of taxa wide
	/**///vector< unsigned int> length_of_bitstrings; //number of biparts long. value correlates to bipartitions
	
	vector< vector <float> > tree_branches;
	
	//for consensus trees
	vector<unsigned int> global_indices;
	vector< float > branches; //seems to only be used for compute tree which is called by the consensus functions
	
	
void print_taxa_in_trees(){
	cout << "taxa_in_trees" << endl;
	for(unsigned int i = 0; i < taxa_in_trees.size(); i++) {
		for (unsigned int j = 0; j < ::NUM_TAXA; j++){
			cout << taxa_in_trees[i][j]; 
		}
		cout << endl; 
	} 
}

bool is_taxa_homogenious(set<unsigned int> treeset){
	if(!::HETERO){
		return true;
	}
	else{
		for (unsigned int i = 0; i < treeset.size(); i++){
			for (unsigned int j = 0; j < lm.size(); j++){
				if (taxa_in_trees[0][j] != taxa_in_trees[i][j]){
					return false;
				}
			}
		}
	}
	return true;
}
vector<string> get_taxa_in_tree(unsigned int treeindex){
	set<string> taxaset;
	for(unsigned int i = 0; i < lm.size(); i ++ ){
		if (taxa_in_trees[treeindex][i] == 1){
			taxaset.insert(lm.name(i));
		}
	}
	vector<string> taxavect(taxaset.begin(), taxaset.end());
	return taxavect;
}

bool are_taxa_in_tree(int treeindex, vector<int> setoftaxa){
	for(unsigned int i = 0; i < setoftaxa.size(); i++){
		if(taxa_in_trees[treeindex][i] == 0){
			return false;
		}
	}
	return true;
}


int num_taxa_in_tree(int treeindex){
	int count = 0;
	for(unsigned int i = 0; i < lm.size(); i ++ ){
		if (taxa_in_trees[treeindex][i] == 1){
			count += 1;
		}
	}
	return count;
}


	//test needed
	vector<string> get_common_taxa(){
		set<string> taxaset = lm.get_all_taxa_set();
		for(unsigned int i = 0; i < taxa_in_trees.size(); i++ ){
			for(unsigned int j = 0; j < ::NUM_TAXA; j++ ){
				if (taxa_in_trees[i][j] == 0){
					taxaset.erase(lm.name(j));
				}
			}
		}
		vector<string> taxavect(taxaset.begin(), taxaset.end());
		return taxavect;
	}
	//test needed
	vector<string> get_mask_for_homog(){
		set<string> taxaset;
		for(unsigned int i = 0; i < taxa_in_trees.size(); i++ ){
			for(unsigned int j = 0; j < ::NUM_TAXA; j++ ){
				if (taxa_in_trees[i][j] == 0){
					taxaset.insert(lm.name(j));
				}
			}
		}
		vector<string> taxavect(taxaset.begin(), taxaset.end());
		return taxavect;
	}


	
	std::vector<unsigned int>::iterator trees_begin(int bitstringindex) { 
		return BipartitionTable[bitstringindex].trees_begin(); 
	}
	std::vector<unsigned int>::iterator trees_end(int bitstringindex) { 
		return BipartitionTable[bitstringindex].trees_end(); 
	}
	std::vector<unsigned int>::reverse_iterator trees_rbegin(int bitstringindex) { 
		return BipartitionTable[bitstringindex].trees_rbegin(); 
	}
	std::vector<unsigned int>::reverse_iterator trees_rend(int bitstringindex) { 
		return BipartitionTable[bitstringindex].trees_rend(); 
	}
	std::vector<unsigned int>::const_iterator trees_cbegin(int bitstringindex) const{
		return BipartitionTable[bitstringindex].trees_cbegin();
	}
	std::vector<unsigned int>::const_iterator trees_cend(int bitstringindex) const{
		return BipartitionTable[bitstringindex].trees_cend();
	}
	std::vector<unsigned int>::const_reverse_iterator trees_crbegin(int bitstringindex) const{
		return BipartitionTable[bitstringindex].trees_crbegin();
	}
	std::vector<unsigned int>::const_reverse_iterator trees_crend(int bitstringindex) const{
		return BipartitionTable[bitstringindex].trees_crend();
	}

	unsigned int bitstring_size(int bitstringindex){
		return BipartitionTable[bitstringindex].bitstring_size();
	}

	bool get_bit(int bitstringindex, int bitindex){
		return BipartitionTable[bitstringindex].get_bit(bitindex);
	}
	
	bool* get_boolarray (int bitstringindex){
		bool *boolarray = new bool[bitstring_size(bitstringindex)];
		for (unsigned int i = 0; i < bitstring_size(bitstringindex); i++){
			boolarray[i] = get_bit(bitstringindex, i);
		}
		return boolarray;
	}
	
	float get_ave_branchlength(int bitstringindex){
		return BipartitionTable[bitstringindex].get_ave_branchlength();
	}
	vector<float> get_branchlengths(int index){
		return BipartitionTable[index].get_branchlengths();
	}
	//I need to return some specific stuff for the compute_tree function
	vector<bool*> get_compute_tree_bipartitions(){
		vector<bool*> bipartitions;
		for (unsigned int i = 0; i < BipartitionTable.size(); i++){
			bipartitions.push_back(get_boolarray(i));
		}
		return bipartitions;
	}
	
	vector<unsigned int> get_compute_tree_bipartitions_bitlens(){
		vector<unsigned int> bipartitionlens;
		for (unsigned int i = 0; i < BipartitionTable.size(); i++){
			bipartitionlens.push_back(bitstring_size(i));
		}
		return bipartitionlens;
	}
	
	vector<float> get_compute_tree_bipartitions_branchlens(){
		vector<float> avebranchlens;
		for (unsigned int i = 0; i < BipartitionTable.size(); i++){
			avebranchlens.push_back(get_ave_branchlength(i));
		}
		return avebranchlens;
	}
	
	vector< unsigned int> get_trees(int bitstringindex){
		return BipartitionTable[bitstringindex].get_trees();
	}
	
	boost::dynamic_bitset<> get_bitstring(int bitstringindex){
		return BipartitionTable[bitstringindex].get_bitstring();
	}
	
	void add_tree(int bitstringindex, int treeindex, float branchlength){
		BipartitionTable[bitstringindex].add_tree(treeindex, branchlength);
	}
	
	unsigned int biparttable_size(){
		return BipartitionTable.size();
	}
	
	unsigned int trees_size(int index){
		return BipartitionTable[index].trees_size();
	}
	unsigned int get_tree(int bitstringindex, int treeindex){
		return BipartitionTable[bitstringindex].get_tree(treeindex);
	}

	bool in_bitstring(int bitstringindex, int place){
		return BipartitionTable[bitstringindex].in_bitstring(place);
	}

	int number_of_ones(int bitstringindex){
		return BipartitionTable[bitstringindex].number_of_ones();
	}
	
	int number_of_zeros(int bitstringindex){
		return BipartitionTable[bitstringindex].number_of_zeros();
	}

	bool is_zero(int bitstringindex, int position){
		return BipartitionTable[bitstringindex].is_zero(position);
	}

	bool same_bitstring_value(int bitstringindex, int position1, int position2){
		return BipartitionTable[bitstringindex].same_bitstring_value(position1, position2);
	}

	bool same_bitstring_value(int bitstringindex, vector<int> positions){
		
		if (is_zero(bitstringindex,positions[0])){
			for (unsigned int i = 0; i < positions.size(); i++){
				//cout << " in is zero i = " << i << endl;
				if(is_one(bitstringindex, positions[i])){
					//cout << "failed on if(is_one(bitstringindex, positions[i])), when i = " << i << endl; 
					//cout << "positions[i] = " << positions[i] <<endl;
					return false;
				}
			}
			return true;
		}
		
		for (unsigned int i = 0; i < positions.size(); i++){
			//cout << " in is one i = " << i << endl;
			if(is_zero(bitstringindex, positions[i])){
				//cout << "failed on if(is_zero(bitstringindex, positions[i])), when i = " << i << endl;  
				return false;
			}
		}
		return true;
	}

	
	bool is_one(int bitstringindex, int position){
		return BipartitionTable[bitstringindex].is_one(position);
	}

	unsigned int get_bipart_table_size(){
		return BipartitionTable.size();
	}

	


	void print_length_of_bitstrings(){
		cout << "We show the hash table below, with the corresponding bitstring reps of the bipartitions:" << endl << endl;
		//keep in mind that hashtable and hash_lengths are global variables
		for (unsigned int i = 0; i < BipartitionTable.size(); i++){
			cout << BipartitionTable[i].bitstring_size() << " ";
		}
		cout << endl;
	}
	
	void print_biparttable(){
		for (unsigned int i = 0; i < BipartitionTable.size(); i++){
			BipartitionTable[i].print_line();
		}
	}
	
	/*
	void print_hashtable(){
		//TODO: Make columns line up 
		cout << "We show the hash table below, with the corresponding bitstring reps of the bipartitions:" << endl << endl;
		//int nbpt = 0;
		int index = 0;
		//keep in mind that hashtable and hash_lengths are global variables
		for (unsigned int i = 0; i < searchtable.size(); i++){ //biparittions
			cout << index << ". ";			
			for (unsigned int j = 0; j < length_of_bitstrings[i]; j++){ //each bit
			  cout << bipartitions[i][j];
			}
			cout << " --> [ ";
			for (unsigned int j = 0; j < searchtable[i].size(); j++){ 
			  cout << searchtable[i][j] << " ";
			}
			cout << "]";
			if (treetable.size() > 0) {
			  cout << " = [";
			  for (unsigned int j = 0; j < ::NUM_TREES; j++) {
			    cout << treetable[i][j];
			  }
			  if (treetable[i][0])
			    {//nbpt++;
				}
			  cout << "]" << endl;
			  index++;
			}
			else {
			  cout << endl;
			}
			//cout << nbpt << endl;
		}
		//print_length_of_bitstrings();
	}
	*/
	
	void print_bitstring(int index){
		BipartitionTable[index].print_bitstring(true);
	}

	boost::dynamic_bitset<> non_trunc_bitstring(int bitstringindex){ //gives the full bitstring (i.e. with 0s at the end)
		boost::dynamic_bitset<> returnVal(lm.size());
		for(unsigned int i = 0; i < BipartitionTable[bitstringindex].bitstring_size(); i++){
			returnVal[i] = BipartitionTable[bitstringindex].get_bit(i);
		}
		return returnVal;
	}
	/*
	bool* full_bitstring(int index){ //gives the full bitstring (i.e. with 0s at the end)
		bool returnVal[::NUM_TAXA];		
		for(int i = 0; i < length_of_bitstrings[index]; i++){
			returnVal[i] = bipartitions[index][i];
			}
		for(int j = length_of_bitstrings[index]; j < ::NUM_TAXA; j++){
			returnVal[j] = 0;
			}
		return returnVal;
		}
	*/
	void print_bitstring_and_trees(int index){
		BipartitionTable[index].print_line();
	}

	void print_treetable() {
	  for (unsigned int i = 0; i < treetable.size(); i++) {
	    cout << "[";
	    for (unsigned int j = 0; j < ::NUM_TREES; j++) {
	      cout << treetable[i][j];
	    }
	    cout << "]";
	  }
	}
	
	
//	void homogenize(){	
//	}
	
	void create_random_bt(){
		
		unsigned int numbiparts = 20;
		unsigned int lengthofbitstring = ::NUM_TAXA;
		unsigned int numberoftrees = 40;
		
		
		
		//BipartitionTable.push_back(b);

		
		for(unsigned int i = 0; i < numbiparts; i++ ){
			Bipartition b(lengthofbitstring);
			//bool *bs = new bool[lengthofbitstring];
			for (unsigned int j = 0; j < lengthofbitstring; j++){
				if (rand() % 2 == 0){
					//bs[j] = (bool)1;
					b.set(j,1); //
				}
				//else{
				//	bs[j] = (bool)0;
				//}
			}
			
			//bipartitions.push_back(bs);		
			
			
			//bool *treebs = new bool[numberoftrees];
			//vector<unsigned int> trees;
			//vector<float> treebranchs;
			for(unsigned int j = 0; j < numberoftrees; j++ ){
				if (rand() % 4 == 0){
					//treebranchs.push_back(1.0);
					//trees.push_back(j);
					//treebs[j] = 1;
					b.add_tree(j,1.0); //
				}
				//else{
				//	treebs[j] = 0;
				//}
			}
			BipartitionTable.push_back(b); //
			
			//tree_branches.push_back(treebranchs);
			//treetable.push_back(treebs);
			//searchtable.push_back(trees);
			//length_of_bitstrings.push_back(lengthofbitstring);
			//branches.push_back(1.0);
			
			b.print_line(); //
		}
		
	}
	
	void apply_taxa_mask(vector<unsigned int> taxa_mask){
		//vector < vector < unsigned int > > searchtable;
		//vector < bool * > treetable;
		//vector < bool * > bipartitions;
		//vector< unsigned int> length_of_bitstrings;
		//vector< float > branches;
		//vector< vector <float> > tree_branches;

	/*	for (unsigned int i = 0; i < bipartitions.size(); i++){ //for each bipartition
			bool *bs = new bool[length_of_bitstrings[i]];
			for(unsigned int j = 0; j < length_of_bitstrings[i]; j++){ //for each bit
				unsigned int place_in_bs = 0;
				if (std::find(taxa_mask.begin(), taxa_mask.end(), j)!=taxa_mask.end()){ //if it's in mask
					continue;
				}
				else{
					bs[place_in_bs] = bipartitions[i][j];
					place_in_bs++;
				}
			}
			bipartitions[i] = bs;
		}
	}*/
	
	/*bool* apply_taxa_mask(bool* bitstring, unsigned int bitstringSize, vector<unsigned int> taxa_mask){
	bool *bs = new bool[::NUM_TAXA];
	
	unsigned int j = 0;
	for (unsigned int i = 0; i < bitstringSize(); i++){
		if (std::find(taxa_mask.begin(), taxa_mask.end(), i)!=taxa_mask.end()){
			bs[j] = bitstring[i];
			j++;
		}
		else{
			continue;
		}
	}
	
	return bs;
}

BipartitionTable make_taxa_masked_bipart_table(vector<unsigned int> taxa_mask){
	
	std::sort(taxa_mask.begin(), taxa_mask.end(), std::greater<unsigned int>());
	
	BipartitionTable masked_bt;
	for (unsigned int i = 0; i < ::biparttable.get_size(); i++){ //for each bipartition
		bool *tempBitstring = ::biparttable.bipartitions[i];
		
		
		//for ( unsigned int j = 0; j < taxa_mask.size(); j++ ){
		//	tempBitstring.erase(tempBitstring.begin() + taxa_mask[j]);
		//} 
		
		masked_bt.bipartitions.push_back(apply_taxa_mask(::biparttable.bipartitions[i],masked_bt.length_of_bitstrings[i], taxa_mask));
		masked_bt.branches.push_back(1.0);
		//masked_bt.length_of_bitstrings.push_back(tempBitstring.size());
		
	}
	masked_bt.print_hashtable();
	return masked_bt;
	* */
}
	
	
	
};

#endif
