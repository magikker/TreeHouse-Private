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

	//by bipartitions
	vector < vector < unsigned int > > searchtable; //number of biparts long by which trees with the bipart wide
	vector < bool * > treetable; //number of biparts long by number of trees wide
	vector < bool * > bipartitions; //number of biparts long by number of taxa wide
	vector< unsigned int> length_of_bitstrings; //number of biparts long. value correlates to bipartitions
	
	vector< vector <float> > tree_branches;

	//for consensus trees
	vector<unsigned int> global_indices;
	vector< float > branches; //seems to only be used for compute tree which is called by the consensus functions
	
	bool in_bitstring(int bitstring, int place){
		//cout << "length_of_bitstrings[bitstring] " << length_of_bitstrings[bitstring] << ", Place = "  << place << endl;
		if(length_of_bitstrings[bitstring] > place){
			//cout << "it's in the bitstring" << endl;
			return true;
		}
		//cout << "it's NOT in the bitstring" << endl;
		return false;
	}
	
	int number_of_ones(int bitstringindex){
		int count = 0;
		for(unsigned int i = 0; i < length_of_bitstrings[bitstringindex]; i ++ ){
			if (bipartitions[bitstringindex][i] == 1){
				count += 1;
			}
		}
		return count;
	}
	
	int number_of_zeros(int bitstringindex){ //NOTE- this only counts the number of zeros of the literal bitstring, not included omitted zeros at the end
		int count = 0;
		for(unsigned int i = 0; i < length_of_bitstrings[bitstringindex]; i ++ ){
			if (bipartitions[bitstringindex][i] == 0){
				count += 1;
			}
		}
		return count;
	}

	bool is_zero(int bitstringindex, int position){
		if( (! in_bitstring(bitstringindex, position)) || bipartitions[bitstringindex][position] == false){
			return true;
		}
		else 
			return false;
	}

	bool same_bitstring_value(int bitstringindex, int position1, int position2){
		bool value;
		
		if(in_bitstring(bitstringindex, position1))
			value = bipartitions[bitstringindex][position1];
		else
			value = false;
		
		if(in_bitstring(bitstringindex, position2))
			return value == bipartitions[bitstringindex][position2];
		else
			return value == false;
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
		//cout << "bitstringindex = " << bitstringindex << ", position = " << position <<endl;
		if(in_bitstring(bitstringindex, position) && bipartitions[bitstringindex][position] == true){
			//cout << "bipartitions[bitstringindex][position] = " << bipartitions[bitstringindex][position] << endl;
			//cout << "It's a one" << endl;
			return true;
		}
		else 
			//cout << "It's a Zero" << endl;
			return false;
	}


	bool is_one(int bitstringindex, int position, unsigned int &bitcomps){
		bitcomps += 1;
		if(in_bitstring(bitstringindex, position) ){
			bitcomps += 1;
			if( bipartitions[bitstringindex][position] == true){
				return true;
			}
		}
		return false;
	}

	bool is_zero(int bitstringindex, int position, unsigned int &bitcomps){
		
		bitcomps += 1;
		if( (! in_bitstring(bitstringindex, position)) ){
			return true;
		}
		bitcomps += 1;
		if ( bipartitions[bitstringindex][position] == false){
			return true;
		}
		return false;
	}
	
	bool same_bitstring_value(int bitstringindex, vector<int> positions, unsigned int &bitcomps){
		
		if (is_zero(bitstringindex,positions[0], bitcomps)){
			for (unsigned int i = 0; i < positions.size(); i++){
				if(is_one(bitstringindex, positions[i], bitcomps)){
					return false;
				}
			}
			return true;
		}
		else{
			for (unsigned int i = 0; i < positions.size(); i++){
				if(is_zero(bitstringindex, positions[i], bitcomps)){
					return false;
				}
			}
			return true;
		}
		return true;
	}

	unsigned int get_size(){
		return searchtable.size();
	}


	void print_length_of_bitstrings(){
		cout << "We show the hash table below, with the corresponding bitstring reps of the bipartitions:" << endl << endl;
		//keep in mind that hashtable and hash_lengths are global variables
		for (unsigned int i = 0; i < length_of_bitstrings.size(); i++){
			cout << length_of_bitstrings[i] << " ";
		}
		cout << endl;
	}
	
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
	
	void print_bitstring(int index){
		for (unsigned int j = 0; j < length_of_bitstrings[index]; j++){
			  cout << bipartitions[index][j];
		}
		cout << endl;
	}

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
	
	void print_bitstring_and_trees(int index){
		for (unsigned int j = 0; j < length_of_bitstrings[index]; j++){
			  cout << bipartitions[index][j] ;
		}
		cout << " --> [ ";
			for (unsigned int j = 0; j < searchtable[index].size(); j++){
			  cout << searchtable[index][j] << " ";
			}
			cout << "]" << endl;
	}

	void print_searchtable(){
		for (unsigned int i = 0; i < searchtable.size(); i++){
			for (unsigned int j = 0; j < searchtable[i].size(); j++){
				cout << searchtable[i][j] << ", ";
			}
			cout << endl;
		}
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
			bool *bs = new bool[lengthofbitstring];
			for (unsigned int j = 0; j < lengthofbitstring; j++){
				if (rand() % 2 == 0){
					bs[j] = (bool)1;
					b.set(j,1); //
				}
				else{
					bs[j] = (bool)0;
				}
			}
			
			bipartitions.push_back(bs);		
			
			
			bool *treebs = new bool[numberoftrees];
			vector<unsigned int> trees;
			vector<float> treebranchs;
			for(unsigned int j = 0; j < numberoftrees; j++ ){
				if (rand() % 4 == 0){
					treebranchs.push_back(1.0);
					trees.push_back(j);
					treebs[j] = 1;
					b.add_tree(j,1.0); //
				}
				else{
					treebs[j] = 0;
				}
			}
			BipartitionTable.push_back(b); //
			
			tree_branches.push_back(treebranchs);
			treetable.push_back(treebs);
			searchtable.push_back(trees);
			length_of_bitstrings.push_back(lengthofbitstring);
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

		for (unsigned int i = 0; i < bipartitions.size(); i++){ //for each bipartition
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
	}

  
	
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
}
	*/
	
	
};

#endif
