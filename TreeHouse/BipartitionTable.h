#ifndef _BIPARTITIONTABLE_H
#define _BIPARTITIONTABLE_H

#include <map>
#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>
#include "Bipartition.h"
#include "TreeSet.h"

using namespace std;

class BipartitionTable {

private:
	
struct classcomp {
	bool operator() (const boost::dynamic_bitset<>& lhs, const boost::dynamic_bitset<>& rhs) const{ 
		if(lhs.count() ==  rhs.count()){
			if(lhs.size() ==  rhs.size()){
				return lhs<rhs;
			}
			else{return lhs.size()<rhs.size();}
		}
		else{return lhs.count()<rhs.count();}
	}
};

	
public:

	std::map< boost::dynamic_bitset<>, TreeSet, classcomp > CladeMap;
	std::vector<std::map< boost::dynamic_bitset<>, TreeSet >::iterator > MapBenchMarks;
	


	vector<Bipartition> BipartTable;
	LabelMap lm;

	//by bipartitions
	vector < boost::dynamic_bitset<> > treetable; //number of biparts long by number of trees wide
	//this treetable is an overlap with the inverted index. We shouldn't need both. 
	
	set<unsigned int> trivial_bipartitions; //keep track of which bipartitions are trivial
	//not sure why this is being kept. 
	

	//by trees
	vector< boost::dynamic_bitset<> > taxa_in_trees; // which taxa are in which trees					NEED
	//needs to be moved to the treetable. 
	
	vector< vector <float> > tree_branches;
	//needs to be moved to the tree table


	//this stuff needs to be phased out. 
	//for consensus trees
	vector<unsigned int> global_indices;
	vector< float > branches; //seems to only be used for compute tree which is called by the consensus functions
	vector< vector< unsigned int> > inverted_index;
	std::map<int, vector<int> > tree_dups;


	//NumTaxa ought to be pulled from LM.size()
	unsigned int NumTrees; //NumTrees as a replacement to NUM_TREES for now... Not sure if it's actually needed. 
	bool weighted;
	bool hetero;
	
	set< unsigned int > get_trees(	std::map< boost::dynamic_bitset<>, TreeSet >::iterator it){
		if (it->second.is_inverse == true){
			return decompress(it);
		}
		return it->second.tree_ids;
	}
	
	boost::dynamic_bitset<> get_bitstring(std::map< boost::dynamic_bitset<>, TreeSet >::iterator it){
		return it->first;
	}
	
	
	set< unsigned int > decompress(std::map< boost::dynamic_bitset<>, TreeSet >::iterator iter){
		std::set<unsigned int> result;
	
		std::set<unsigned int>::iterator invit=iter->second.tree_ids.begin(); 
		unsigned int i = 0;
		while (i!= NumTrees && invit!=iter->second.tree_ids.end()){
			if(i < *invit){
				//std::cout << i << ' ' ;
				result.insert(i);
				++i;
			}
			else if(*invit < i){
				++invit;
			} 
			else{ 
				++invit; 
				++i;
			}
		}
		while (i!= NumTrees){
			result.insert(i);
			//std::cout << i << ' ' ;
			++i;
		}	
		return result;
	}
	
	bool contains_tree(std::map< boost::dynamic_bitset<>, TreeSet >::iterator it, unsigned int treeid){
		
		if (it->second.is_inverse == true){
			if(it->second.tree_ids.find(treeid) != it->second.tree_ids.end()){
				return false;
			}
			else{
				return true;
			}
		}
		else{
			if(it->second.tree_ids.find(treeid) != it->second.tree_ids.end()){
				return true;
			}
		}
		return false;
	}
	
	
	
	
	
//I need some work on this one. 
set <unsigned int> duplicates(int treein){
	//Returns the unique trees in a tree set. 
	set<unsigned int> duplicatelist;
	
	if (tree_dups[treein].size() > 0 ){
		if(tree_dups[treein][0] < treein){
			vector<int> tempvect = tree_dups[tree_dups[treein][0]];
			copy(tempvect.begin(), tempvect.end(), inserter(duplicatelist, duplicatelist.end()));
			//duplicatelist = ::dups[dups[treein][0]];
			duplicatelist.insert(tree_dups[treein][0]);
		}
		
		else{
			vector<int> tempvect = tree_dups[treein];
			copy(tempvect.begin(), tempvect.end(), inserter(duplicatelist, duplicatelist.end()));
			//duplicatelist = ::dups[treein];
		}
	}
	return duplicatelist;
}

void print_taxa_in_trees(){
	cout << "taxa_in_trees" << endl;
	for(unsigned int i = 0; i < taxa_in_trees.size(); i++) {
		for (unsigned int j = 0; j < taxa_in_trees[i].size(); j++){
			cout << taxa_in_trees[i][j]; 
		}
		cout << endl; 
	} 
}

bool is_taxa_homogenious(set<unsigned int> treeset){
	//if(!::HETERO){
	//	return true;
	//}
	//else{
		for (unsigned int i = 0; i < treeset.size(); i++){
			for (unsigned int j = 0; j < lm.size(); j++){
				if (taxa_in_trees[0][j] != taxa_in_trees[i][j]){
					return false;
				}
			}
		}
	//}
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
			for(unsigned int j = 0; j < taxa_in_trees[i].size(); j++ ){
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
			for(unsigned int j = 0; j < taxa_in_trees[i].size(); j++ ){
				if (taxa_in_trees[i][j] == 0){
					taxaset.insert(lm.name(j));
				}
			}
		}
		vector<string> taxavect(taxaset.begin(), taxaset.end());
		return taxavect;
	}


	
	std::vector<unsigned int>::iterator trees_begin(int bitstringindex) { 
		return BipartTable[bitstringindex].trees_begin(); 
	}
	std::vector<unsigned int>::iterator trees_end(int bitstringindex) { 
		return BipartTable[bitstringindex].trees_end(); 
	}
	std::vector<unsigned int>::reverse_iterator trees_rbegin(int bitstringindex) { 
		return BipartTable[bitstringindex].trees_rbegin(); 
	}
	std::vector<unsigned int>::reverse_iterator trees_rend(int bitstringindex) { 
		return BipartTable[bitstringindex].trees_rend(); 
	}
	std::vector<unsigned int>::const_iterator trees_cbegin(int bitstringindex) const{
		return BipartTable[bitstringindex].trees_cbegin();
	}
	std::vector<unsigned int>::const_iterator trees_cend(int bitstringindex) const{
		return BipartTable[bitstringindex].trees_cend();
	}
	std::vector<unsigned int>::const_reverse_iterator trees_crbegin(int bitstringindex) const{
		return BipartTable[bitstringindex].trees_crbegin();
	}
	std::vector<unsigned int>::const_reverse_iterator trees_crend(int bitstringindex) const{
		return BipartTable[bitstringindex].trees_crend();
	}
	
	std::vector<float>::iterator branchlengths_begin(int bitstringindex) { 
		return BipartTable[bitstringindex].branchlengths_begin(); 
	}
	std::vector<float>::iterator branchlengths_end(int bitstringindex) { 
		return BipartTable[bitstringindex].branchlengths_end(); 
	}
	std::vector<float>::reverse_iterator branchlengths_rbegin(int bitstringindex) { 
		return BipartTable[bitstringindex].branchlengths_rbegin(); 
	}
	std::vector<float>::reverse_iterator branchlengths_rend(int bitstringindex) { 
		return BipartTable[bitstringindex].branchlengths_rend(); 
	}
	std::vector<float>::const_iterator branchlengths_cbegin(int bitstringindex) const{
		return BipartTable[bitstringindex].branchlengths_cbegin();
	}
	std::vector<float>::const_iterator branchlengths_cend(int bitstringindex) const{
		return BipartTable[bitstringindex].branchlengths_cend();
	}
	std::vector<float>::const_reverse_iterator branchlengths_crbegin(int bitstringindex) const{
		return BipartTable[bitstringindex].branchlengths_crbegin();
	}
	std::vector<float>::const_reverse_iterator branchlengths_crend(int bitstringindex) const{
		return BipartTable[bitstringindex].branchlengths_crend();
	}	
	
	

	unsigned int bitstring_size(int bitstringindex){
		return BipartTable[bitstringindex].bitstring_size();
	}

	bool get_bit(int bitstringindex, int bitindex){
		return BipartTable[bitstringindex].get_bit(bitindex);
	}
	
	bool* get_boolarray (int bitstringindex){
		bool *boolarray = new bool[bitstring_size(bitstringindex)];
		for (unsigned int i = 0; i < bitstring_size(bitstringindex); i++){
			boolarray[i] = get_bit(bitstringindex, i);
		}
		return boolarray;
	}
	
	float get_ave_branchlength(int bitstringindex){
		return BipartTable[bitstringindex].get_ave_branchlength();
	}
	
	vector<float> get_branchlengths(int index){
		return BipartTable[index].get_branchlengths();
	}

	vector<unsigned int> get_ones_indices(int index){
		return BipartTable[index].get_ones_indices();
	}
	
	
	std::set<unsigned int> getOnes(int index){
		return BipartTable[index].getOnes();
	}
	
	std::set<unsigned int> getZeros(int index){
		return BipartTable[index].getZeros();
	}
	
	
	//I need to return some specific stuff for the compute_tree function
	vector<boost::dynamic_bitset<> > get_compute_tree_bipartitions(){
		vector<boost::dynamic_bitset<> > bipartitions;
		for (unsigned int i = 0; i < BipartTable.size(); i++){
			bipartitions.push_back(get_bitstring(i));
		}
		return bipartitions;
	}
	/*
	vector<unsigned int> get_compute_tree_bipartitions_bitlens(){
		vector<unsigned int> bipartitionlens;
		for (unsigned int i = 0; i < BipartTable.size(); i++){
			bipartitionlens.push_back(bitstring_size(i));
		}
		return bipartitionlens;
	}
	*/
	vector<float> get_compute_tree_bipartitions_branchlens(){
		vector<float> avebranchlens;
		for (unsigned int i = 0; i < BipartTable.size(); i++){
			avebranchlens.push_back(get_ave_branchlength(i));
		}
		return avebranchlens;
	}
		
	vector< unsigned int> get_trees(int bitstringindex){
		return BipartTable[bitstringindex].get_trees();
	}
	
	boost::dynamic_bitset<> get_bitstring(int bitstringindex){
		return BipartTable[bitstringindex].get_bitstring();
	}
	
	void add_tree(int bitstringindex, int treeindex, float branchlength){
		BipartTable[bitstringindex].add_tree(treeindex, branchlength);
	}
	
	unsigned int biparttable_size(){
		return BipartTable.size();
	}
	
	unsigned int trees_size(int index){
		return BipartTable[index].trees_size();
	}
	unsigned int get_tree(int bitstringindex, int treeindex){
		return BipartTable[bitstringindex].get_tree(treeindex);
	}

	bool in_bitstring(int bitstringindex, int place){
		return BipartTable[bitstringindex].in_bitstring(place);
	}

	int number_of_ones(int bitstringindex){
		return BipartTable[bitstringindex].number_of_ones();
	}
	
	//number of zeros might not be meaningful as bitstring len varies. 
	int number_of_zeros(int bitstringindex){		
		return BipartTable[bitstringindex].number_of_zeros();
	}

	bool is_zero(int bitstringindex, int position){
		return BipartTable[bitstringindex].is_zero(position);
	}

	bool are_zeros(int bitstringindex, vector<int> positions){
		for (unsigned int i = 0; i < positions.size(); i++){
			if(is_one(bitstringindex, positions[i])){
				return false;
			}
		}
		return true;
	}

	bool same_bitstring_value(int bitstringindex, int position1, int position2){
		return BipartTable[bitstringindex].same_bitstring_value(position1, position2);
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
		return BipartTable[bitstringindex].is_one(position);
	}

	bool are_ones(int bitstringindex, vector<int> positions){
		for (unsigned int i = 0; i < positions.size(); i++){
			if(is_zero(bitstringindex, positions[i])){
				return false;
			}
		}
		return true;
	}

	unsigned int get_bipart_table_size(){
		return BipartTable.size();
	}

	

	void print_length_of_bitstrings(){
		cout << "We show the hash table below, with the corresponding bitstring reps of the bipartitions:" << endl << endl;
		//keep in mind that hashtable and hash_lengths are global variables
		for (unsigned int i = 0; i < BipartTable.size(); i++){
			cout << BipartTable[i].bitstring_size() << " ";
		}
		cout << endl;
	}
	
	void print_biparttable(){
		for (unsigned int i = 0; i < BipartTable.size(); i++){
			cout << i << ". ";
			BipartTable[i].print_line();
		}
	}

	void print_inverted_index(){
		cout << "Trivial Bipartitions:\n";
		for(set<unsigned int>::iterator it = trivial_bipartitions.begin(); it!=trivial_bipartitions.end(); it++){
			cout << *it << "  ";
			}
			cout << endl << endl;
		for(unsigned int i = 0; i < inverted_index.size(); i++){
			cout << i << ". ";
			for(unsigned int j = 0; j < inverted_index.at(i).size(); j++){
				unsigned int bipart = inverted_index.at(i).at(j);
				if(trivial_bipartitions.find(bipart)==trivial_bipartitions.end()){ //if the bipart is non-trivial
					cout << inverted_index.at(i).at(j) << ", ";
					}
				}
				cout << endl;	
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
	
	void print_clade(int index){
		for (unsigned int i = 0; i < BipartTable[index].bitstring_size(); i++){
			if (BipartTable[index].is_one(i)){
				cout << lm.name(i) << ' ' ;
			}
		}
		cout << endl;
	}
	
	void print_clade_table(){
		for (unsigned int i = 0; i < BipartTable.size(); i++){
			print_clade(i);
			cout << ": ";
			BipartTable[i].print_trees(true);
			}
	}
	
	void print_bipartitions(int index){
		cout << "Tree " << index << " contains:" <<endl; 
		for (unsigned int i = 0; i < treetable.size(); i++ ){
			if(treetable[i][index] == 1){
				print_bitstring(i);
			}
		}
	}

	void print_clades(int index){
		cout << "Tree " << index << " contains:" <<endl; 
		for (unsigned int i = 0; i < treetable.size(); i++ ){
			if(treetable[i][index] == 1){
				print_clade(i);
			}
		}
	}
	
	
	/*
	void print_CladeMap_trees(){
		typedef std::map<std::string, std::map<std::string, std::string>>::iterator it_type;
		for(it_type iterator = m.begin(); iterator != m.end(); iterator++) {
			// iterator->first = key
			// iterator->second = value
			// Repeat if you also want to iterate through the second map.
		} 
		
		if(CladeMap.inverse_ids == false){
			// std::cout << "myset contains:";
			for (std::set<unsigned int>::iterator it=CladeMap.tree_ids.begin(); it!=CladeMap.tree_ids.end(); ++it){
				std::cout << ' ' << *it;
			}
			std::cout << '\n';
		}
	
	else{
		
	}
}
* 
* 
* 	while (first1!=last1 && first2!=last2)
  {
    if (*first1<*first2) { *result = *first1; ++result; ++first1; }
    else if (*first2<*first1) ++first2;
    else { ++first1; ++first2; }
  }
  return std::copy(first1,last1,result);
}
		
	
* 
*/	
	
	void print_data_info(){
	
		unsigned int unique_clades = 0;
		unsigned int ttids = 0;
		unsigned int storedtid = 0;
		unsigned int compressedlines = 0;
		vector<unsigned int> clade_distrobution;
		typedef std::map< boost::dynamic_bitset<>, TreeSet >::iterator clade_it_type;

		for(unsigned int i = 0; i < lm.size(); i++){
			clade_distrobution.push_back(0);
		}
		
		for(clade_it_type iterator = CladeMap.begin(); iterator != CladeMap.end(); iterator++) {
			clade_distrobution[iterator->first.count()] += 1;
			unique_clades += 1;
			storedtid += iterator->second.tree_ids.size();
			if  (iterator->second.is_inverse == true){
				compressedlines += 1;
				ttids += (NumTrees - iterator->second.tree_ids.size());
			}
			else{
				ttids += iterator->second.tree_ids.size();
			}
		}
	
		cout << "unique_clades = " << unique_clades << endl;
		cout << "ttids = " << ttids << endl;
		cout << "storedtid = " << storedtid << endl;
		cout << "compressedlines = " << compressedlines << endl;
		cout << "numtaxa = " << lm.size() << endl;

		cout << "clade distrobution" << endl;
		
		cout << "[";
		for(unsigned int i = 0; i < lm.size(); i++){
			cout << i << ",";
		}
		cout << "]" << endl;
		
		
		cout << "[";
		for(unsigned int i = 0; i < lm.size(); i++){
			cout << clade_distrobution[i] << ",";
		}
		cout << "]"<< endl;
		

	}

	void print_CladeMap(){
		
		typedef std::map< boost::dynamic_bitset<>, TreeSet >::iterator clade_it_type;

		set <unsigned int> all_trees;
		for(unsigned int i = 0; i < NumTrees; i++){
			all_trees.insert(i);
		}

		
		//typedef std::map<std::string, std::map<std::string, std::string>>::iterator it_type;
		for(clade_it_type iterator = CladeMap.begin(); iterator != CladeMap.end(); iterator++) {
			cout << iterator->second.tree_ids.size();
			if  (iterator->second.is_inverse == false){
				cout << " + "; 
			}
			if (iterator->second.is_inverse == true){
				cout << " - "; 
			}
			
			cout << iterator->first;
			cout << " : [";
			if(iterator->second.is_inverse == true){
				//cout << "in a -";
				std::set<unsigned int>::iterator invit=iterator->second.tree_ids.begin(); 
				std::set<unsigned int>::iterator allit = all_trees.begin(); 
		
				while (allit!= all_trees.end() && invit!=iterator->second.tree_ids.end()){
					//cout << "in the loop";
					if(*allit < *invit){
						std::cout << *allit << ' ' ;
						++allit;
					}
					else if(*invit < *allit){
						++invit;
					} 
					else{ 
						++invit; 
						++allit;
					}
				}
				while (allit!= all_trees.end()){
					std::cout << *allit << ' ' ;
					++allit;
				}	
			}
			else{
				//cout << "IN a +";
				for (std::set<unsigned int>::iterator it=iterator->second.tree_ids.begin(); it!=iterator->second.tree_ids.end(); ++it){
					std::cout << *it << ' ' ;
				}				
			}
			cout << "]" << endl; 
		// iterator->second = value
		// Repeat if you also want to iterate through the second map.
		} 
	print_data_info();
	}
	
	
	
	
	void print_bitstring(int index){
		BipartTable[index].print_bitstring(true);
	}



	boost::dynamic_bitset<> non_trunc_bitstring(int bitstringindex){ //gives the full bitstring (i.e. with 0s at the end)
		boost::dynamic_bitset<> returnVal(lm.size());
		for(unsigned int i = 0; i < BipartTable[bitstringindex].bitstring_size(); i++){
			returnVal[i] = BipartTable[bitstringindex].get_bit(i);
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
		BipartTable[index].print_line();
	}
	

	
	
	void print_treetable() {
	  for (unsigned int i = 0; i < treetable.size(); i++) {
	    cout << "[";
	    for (unsigned int j = 0; j < treetable[i].size(); j++) {
	      cout << treetable[i][j];
	    }
	    cout << "]";
	  }
	}
	
	bool is_trivial(int index){
		cout << "Warning, function called that's not heterogenous safe" << endl;
		if(number_of_ones(index)==1 || number_of_ones(index)==(lm.size()-1)){
			return true;
		}
		return false;
	}

	void calculate_trivial_bipartitions(){
		cout<< "calculate_trivial_bipartitions: warning function call not hetero safe"<<endl;
		for(unsigned int i = 0; i < BipartTable.size(); i++){
			if(number_of_ones(i)<2 || ((number_of_zeros(i)+lm.size()-BipartTable.at(i).bitstring_size()) < 2)){
				trivial_bipartitions.insert(i);
			}
		}
	}
	
//	void homogenize(){	
//	}
	
	void create_random_bt(){
		
		unsigned int numbiparts = 20;
		unsigned int lengthofbitstring = lm.size();
		unsigned int numberoftrees = 40;
		
		
		
		//BipartTable.push_back(b);

		
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
			BipartTable.push_back(b); //
			
			//tree_branches.push_back(treebranchs);
			//treetable.push_back(treebs);
			//searchtable.push_back(trees);
			//length_of_bitstrings.push_back(lengthofbitstring);
			//branches.push_back(1.0);
			
			b.print_line(); //
		}
		
	}
	
	vector<unsigned int> homog_taxa(){
		vector<unsigned int> taxa_indices; 
		
		for(unsigned int i = 0; i < lm.size(); i++ ){
			bool flag = true;
			for(unsigned int j = 0; j < taxa_in_trees.size(); j++ ){
				if (taxa_in_trees[j][i] == 0){ //a tree doesn't have it
					flag = false;
					break;
				}
			}
			if (flag == false){
				taxa_indices.push_back(i);
			}
		}
		
		return taxa_indices;
	  //cout << "the taxa not contained by all trees are "<< endl;
		
	  //for (unsigned int i = 0; i < taxa_indices.size(); i++){
      //      cout << taxa_indices[i] << " ";
      //}
	
		
	}
	
	
	void apply_taxa_mask(vector<unsigned int> taxa_mask){
		//remove taxa from bipartitions
		for (unsigned int i = 0; i < BipartTable.size(); i++){ //for each bipartition
			for (unsigned int j = 0; j < taxa_mask.size(); j++){
				BipartTable[i].set(taxa_mask[j],0);
			}
		}
		
		//removed now duplicated bipartitions
		for (unsigned int i = 0; i < BipartTable.size(); i++){ //for each bipartition
			vector<unsigned int> matches;
			for (unsigned int j = i+1; j < BipartTable.size(); j++){
				if(BipartTable[i].getOnes() == BipartTable[j].getOnes()){
					cout << "found a match between" << endl;
					BipartTable[i].print_line();
					BipartTable[j].print_line();
					matches.push_back(j);
					for(unsigned int k = 0; k < BipartTable[j].trees_size(); k++){
						BipartTable[i].add_tree(BipartTable[j].get_tree(k),BipartTable[j].get_branchlength(k));
					}
					
				}
				if(BipartTable[j].number_of_ones() == 0){
					matches.push_back(j);
				}
				
			}
			
			while (!matches.empty()){
				unsigned int match = matches.back();
				//for(unsigned int k = 0; k < BipartTable[match].trees_size(); k++){
				//	BipartTable[i].add_tree(BipartTable[match].get_tree(k),BipartTable[match].get_branchlength(k));
				//}
				matches.pop_back();
				BipartTable.erase(BipartTable.begin()+match);
			}
	
		}
		
		//Fix taxa in trees
		for(unsigned int i = 0; i < taxa_in_trees.size(); i++){
			for(unsigned int j = 0; j <  taxa_mask.size(); j++){
				taxa_in_trees[i].set(taxa_mask[j],0);
			}
			
		}

		//Fix treetable
		treetable.clear();
		for (unsigned int i = 0; i < BipartTable.size(); i++){ //for each bipartition
			boost::dynamic_bitset<> bits(NumTrees); 
			vector<unsigned int> trees = BipartTable[i].get_trees();
			for (unsigned int j = 0; j < trees.size(); j++){
				bits.set(trees[j],1);
			}
			treetable.push_back(bits);
		}
		hetero = false;
	}
};

#endif
