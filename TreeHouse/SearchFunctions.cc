#include "SearchFunctions.h"
using namespace std;

//GRB NEW
set<unsigned int> search_clade(vector<string> RequiredTaxa) {
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return search_clade(required);
}

//GRB NEW // Direct look up in log(n) time. Het Safe
set<unsigned int> search_clade(vector<int> required){
	set<unsigned int> result_trees;
	
	if (required.size() == 0){ 
		return result_trees;
	}
	
	boost::dynamic_bitset<> search_bs(biparttable.lm.size());
	for (vector<int>::iterator iter = required.begin(); iter != required.end(); iter++){
		search_bs.set(*iter, true);
	}
	
	map<boost::dynamic_bitset<>,TreeSet>::iterator it=biparttable.CladeMap.find(search_bs);
	
	if (it != biparttable.CladeMap.end()){
		::SetOps += 1;
		start_clock();
		result_trees = biparttable.get_trees(it);
		::SetTime += stop_clockbp();
		::SetInsertions += result_trees.size();
		}
	
	return result_trees;
	
}

//GRB NEW
set<unsigned int> search_bl(vector<string> RequiredTaxa, string GreaterOrLessThan, float bl) {
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return search_bl(required, GreaterOrLessThan, bl);
}

//GRB NEW // Direct look up in log(n) time. Het Safe
set<unsigned int> search_bl(vector<int> required, string GreaterOrLessThan, float bl){
	set<unsigned int> result_trees;
	vector<float> branchlengths;
	
	if (required.size() == 0){ 
		return result_trees;
	}
	
	boost::dynamic_bitset<> search_bs(biparttable.lm.size());
	for (vector<int>::iterator iter = required.begin(); iter != required.end(); iter++){
		search_bs.set(*iter, true);
	}
	
	map<boost::dynamic_bitset<>,TreeSet>::iterator it=biparttable.CladeMap.find(search_bs);
	
	if (it != biparttable.CladeMap.end()){
		result_trees = biparttable.get_trees(it);
		branchlengths = biparttable.get_branchlengths(it);
	}

	if (result_trees.size() != branchlengths.size() ){
		cout << "something is really wrong" << endl;
	}
	
	vector<unsigned int> temptrees(result_trees.begin(), result_trees.end());
	result_trees.clear();

	for (unsigned int i =0; i < temptrees.size(); i++){
		if(GreaterOrLessThan == ">"){
			//cout << branchlengths[i] << " > " << bl << endl;
			if (branchlengths[i] > bl){
				result_trees.insert(temptrees[i]);
			}
		}
		else if(GreaterOrLessThan == "<"){
			//cout << branchlengths[i] << " < " << bl << endl;
			if (branchlengths[i] < bl){
				result_trees.insert(temptrees[i]);
			}	
		}
	}
	return result_trees;
}

//GRB NEW // x direct lookups. x log(n) time. Het Safe
//Executes the lookup work for the subtree search. 
set<unsigned int> search_by_multiple_clades(vector<vector<int>> required){
	
	//cout << "required" << required.size() << endl;
	
	vector<vector<int>>::iterator iter =required.begin(); 
	set<unsigned int> result_trees = search_clade(*iter);
	iter++;
	//cout << "Out" << result_trees.size() << endl;
	for(;iter != required.end(); iter++ ){
		set<unsigned int> found_trees;
		set<unsigned int> temp;
		found_trees = search_clade(*iter);
		std::set_intersection(result_trees.begin(), result_trees.end(), found_trees.begin(), found_trees.end(),  std::inserter(temp, temp.end()));
		result_trees.swap(temp);
		//cout << "In" << result_trees.size() << endl;
	}

	return result_trees;
	
}

//GRB NEW
set<unsigned int> search_bl_outliers(vector<string> RequiredTaxa, string GreaterOrLessThan, float bl) {
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return search_bl_outliers(required, GreaterOrLessThan, bl);
}

//GRB NEW // Direct look up in log(n) time. Het Safe
set<unsigned int> search_bl_outliers(vector<int> required, string GreaterOrLessThan, float bl){
	set<unsigned int> result_trees;
	vector<float> branchlengths;
	
	if (required.size() == 0){ 
		return result_trees;
	}
	
	boost::dynamic_bitset<> search_bs(biparttable.lm.size());
	for (vector<int>::iterator iter = required.begin(); iter != required.end(); iter++){
		search_bs.set(*iter, true);
	}
	
	map<boost::dynamic_bitset<>,TreeSet>::iterator it=biparttable.CladeMap.find(search_bs);
	
	if (it != biparttable.CladeMap.end()){
		result_trees = biparttable.get_trees(it);
		branchlengths = biparttable.get_branchlengths(it);
	}


	double sum = std::accumulate(branchlengths.begin(), branchlengths.end(), 0.0);
	double mean = sum / branchlengths.size();
	double sq_sum = std::inner_product(branchlengths.begin(), branchlengths.end(), branchlengths.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / branchlengths.size() - mean * mean);		

	if (result_trees.size() != branchlengths.size() ){
		cout << "something is really wrong" << endl;
	}
	
	vector<unsigned int> temptrees(result_trees.begin(), result_trees.end());
	result_trees.clear();

	for (unsigned int i = 0; i < temptrees.size(); i++){
		if(GreaterOrLessThan == ">"){
			//cout << branchlengths[i] << " > " << bl << endl;
			if (branchlengths[i] > (bl*stdev)){
				result_trees.insert(temptrees[i]);
			}
		}
		else if(GreaterOrLessThan == "<"){
			//cout << branchlengths[i] << " < " << bl << endl;
			if (branchlengths[i] < (bl*stdev)){
				result_trees.insert(temptrees[i]);
			}	
		}
	}
	return result_trees;
}


//GRB NEW // x direct lookups of log(n) time. X is number of clades in subtree. // het safe
// Only works with a subtree match, not a full tree. This is because we don't treat a whole tree as a clade, but we capture the "whole tree clade here."
set<unsigned int> get_trees_by_subtree(string subtree)
{
	//cout << "We show the string that repersent the subtree:" << endl;
	//cout << subtree << endl;
	
	set<unsigned int> trees;

	NEWICKTREE *newickTree;
	int err;
	char * cs = strdup(subtree.c_str());
	newickTree = stringloadnewicktree(cs, &err);

	vector< vector < int > > subtrees; 

	if (!newickTree) {
		switch (err) {
			case -1:
			cout << "Out of memory" << endl;
			break;
		case -2:
			cout << "parse error" << endl;
			break;
		case -3:
			cout << "Can't load file" << endl;
			break;
		default:
			cout << "Error " << err << endl;
		exit(0);
		}
    }

    else {
		//unsigned int numBitstr=0;
		bool * bs = dfs_compute_bitstrings(newickTree->root, NULL, subtrees);
		delete[] bs;
		
  		//cout << "We show bitstrings that repersent the subtree:" << endl;
  		//for (unsigned int i = 0; i < subtrees.size(); i++)
		//{
    	//	for (unsigned int j = 0; j < subtrees[i].size(); j++)
		//	{
      	//		cout << subtrees[i][j] << " " ;
    	//	}
    	//	cout << endl;
  		//}
		
		vector< vector < int > >::iterator it = subtrees.begin() ;
		while(it != subtrees.end() ){
			if (it->size() == 1){
				it = subtrees.erase(it);
			}
			else{
				++it;
			}
		}
		
  		//cout << "We show bitstrings that repersent the subtree:" << endl;
  		//for (unsigned int i = 0; i < subtrees.size(); i++)
		//{
    	//	for (unsigned int j = 0; j < subtrees[i].size(); j++)
		//	{
      	//		cout << subtrees[i][j] << " " ;
    	//	}
    	//	cout << endl;
  		//}
		
		trees = search_by_multiple_clades(subtrees);
		//print_vector_of_bs(bitstrings, ::NUM_TAXA);
		//if (MAXVAL == 0)
		//{
		//	MAXVAL = numBitstr -1;
		//}		
		killnewicktree(newickTree);
    }
    free(cs);
	return trees;
}

set<unsigned int> search_ktet(vector<string> leftside, vector<string> rightside){
	vector<int> leftside_required = ::biparttable.lm.lookUpLabels(leftside);
	vector<int> rightside_required = ::biparttable.lm.lookUpLabels(rightside);

	return search_ktet(leftside_required, rightside_required);
}

//GRB NEW
set<unsigned int> search_subclade(vector<string> RequiredTaxa) {
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return search_subclade(required);
}

//GRB NEW
set<unsigned int> search_subclade(vector<int> required){

	set<unsigned int> trees;
	set<unsigned int> matchingtrees;
	typedef std::map< boost::dynamic_bitset<>, TreeSet >::iterator clade_it_type;

	
	if (required.size() > biparttable.lm.size() ){
		cout << "error searching for more taxa than are total in in the set" <<endl;
		return trees;
	}
	
	boost::dynamic_bitset<> required_search_bs(biparttable.lm.size());
	for (vector<int>::iterator iter = required.begin(); iter != required.end(); iter++){
		//cout << *iter << " ";
		required_search_bs.set(*iter, true);
	}
	
	boost::dynamic_bitset<> temp_bs(biparttable.lm.size());
	
	if(required.size() > 0){
		for(clade_it_type iter = biparttable.MapBenchMarks[required_search_bs.count()]; iter != biparttable.MapBenchMarks[biparttable.lm.size()]; iter++) {
		//for(clade_it_type iter = biparttable.MapBenchMarks[0]; iter != biparttable.MapBenchMarks[biparttable.lm.size()]; iter++) {
		
		//cout  << "LEFt iter->first.count() = " << iter->first.count() << endl;
			if (trees.size() == biparttable.NumTrees){
				break;
			}
			if (required_search_bs.is_subset_of(iter->first)){
				matchingtrees = ::biparttable.get_trees(iter);
				::SetOps += 1;
				::SetInsertions += matchingtrees.size();
				start_clock();
				trees.insert(matchingtrees.begin(), matchingtrees.end());
				::SetTime += stop_clockbp();
			}
		}
	}	
	return trees;
}

//GRB NEW
set<unsigned int> search_ktet(vector<int> leftside, vector<int> rightside){
	//cout << "We show the bipartitions found by the search and corresponding trees" << endl;
	//keep in mind that hashtable and hash_lengths are global variables
	//cout << "hetero = " << ::biparttable.hetero << endl;
	
	//cout << "leftside.size() = " << leftside.size() << endl;
	//for (int i = 0; i < leftside.size(); i++)
	//	cout << leftside[i] << " ";
	//cout << endl;

	//cout << "rightside.size() = " << rightside.size() << endl;
	//for (int i = 0; i < rightside.size(); i++)
	//	cout << rightside[i] << " ";
	//cout << endl;
	
	//cout << "biparttable.lm.size() = " << biparttable.lm.size() << endl;
	
	
	set<unsigned int> trees;
	set<unsigned int> matchingtrees;
	typedef std::map< boost::dynamic_bitset<>, TreeSet >::iterator clade_it_type;

	/*
	set <unsigned int> all_trees;
	for(unsigned int i = 0; i < NumTrees; i++){
		all_trees.insert(i);
	}
	*/
	
	
	if (leftside.size() + rightside.size() > biparttable.lm.size() ){
		cout << "error searching for more taxa than are total in in the set" <<endl;
		return trees;
	}
	
	boost::dynamic_bitset<> left_search_bs(biparttable.lm.size());
	for (vector<int>::iterator iter = leftside.begin(); iter != leftside.end(); iter++){
		//cout << *iter << " ";
		left_search_bs.set(*iter, true);
	}
	
	boost::dynamic_bitset<> right_search_bs(biparttable.lm.size());
	for (vector<int>::iterator iter = rightside.begin(); iter != rightside.end(); iter++){
		//cout << *iter << " ";
		right_search_bs.set(*iter, true);
	}

	boost::dynamic_bitset<> temp_bs(biparttable.lm.size());

	//cout <<"leftside.size() " <<leftside.size() << "rightside.size() "<< rightside.size()  << endl;
	
	if(leftside.size() > 0){
		for(clade_it_type iter = biparttable.MapBenchMarks[left_search_bs.count()]; iter != biparttable.MapBenchMarks[biparttable.lm.size() - right_search_bs.count() +1  ]; iter++) {
			//cout  << "LEFt iter->first.count() = " << iter->first.count() << endl;
			if (trees.size() == biparttable.NumTrees){
				break;
			}
			if (left_search_bs.is_subset_of(iter->first)){
				temp_bs = right_search_bs;
				temp_bs &= iter->first;
				if( !temp_bs.any() ){
					//cout << "yep " << iter->first << endl;
					matchingtrees = ::biparttable.get_trees(iter);
					::SetOps += 1;
					::SetInsertions += matchingtrees.size();
					start_clock();
					trees.insert(matchingtrees.begin(), matchingtrees.end());
					::SetTime += stop_clockbp();

				}
			}
		}
	}
	
	if(rightside.size() > 0){
		for(clade_it_type iter = biparttable.MapBenchMarks[right_search_bs.count()]; iter != biparttable.MapBenchMarks[biparttable.lm.size() - left_search_bs.count() +1 ]; iter++) {
			//cout  << "Right iter->first.count() = " << iter->first.count() << endl;
			if (trees.size() == biparttable.NumTrees){
				break;
			}
			if (right_search_bs.is_subset_of(iter->first)){
				temp_bs = left_search_bs;
				temp_bs &= iter->first;
				if( !temp_bs.any() ){
					//cout << "yep " << iter->first << endl;
					matchingtrees = ::biparttable.get_trees(iter);
					::SetOps += 1;
					::SetInsertions += matchingtrees.size();
					start_clock();
					trees.insert(matchingtrees.begin(), matchingtrees.end());
					::SetTime += stop_clockbp();

				}
			}
		}
	}	
	
	if (::biparttable.hetero == true){ // get the set of trees that contain all searched taxa and do a set difference with matches. 
		start_clock();
		vector<int> searched_taxa; 
		searched_taxa.reserve( leftside.size() + rightside.size() ); // preallocate memory
		searched_taxa.insert( searched_taxa.end(), leftside.begin(), leftside.end() );
		searched_taxa.insert( searched_taxa.end(), rightside.begin(), rightside.end() );
		set<unsigned int> trees_wo_taxa = get_trees_without_taxa(searched_taxa);		
	
		std::set<unsigned int> s; 

		std::set_difference(trees.begin(), trees.end(), trees_wo_taxa.begin(), trees_wo_taxa.end(),  std::inserter(s, s.end()));
 	
		//trees.swap(s);
		trees = s;
		::HetTime += stop_clockbp();
	}
	
	return trees;

}













set<unsigned int> get_subset_trees(int tree){ //return trees which are subsets of the given tree, 
						//i.e. trees whose only bipartitions are contained within the given tree
  set<unsigned int> retSet;
  
  vector<unsigned int> b = ::biparttable.inverted_index.at(tree);
  set<unsigned int> treeBiparts(b.begin(), b.end());
  for(unsigned int i = 0; i < ::biparttable.NumTrees; i++){
	if(i!=tree){ //don't include itself 
		vector<unsigned int> biparts = ::biparttable.inverted_index.at(i);
	 	set<unsigned int> bipartSet(biparts.begin(), biparts.end());
		if(includes(treeBiparts.begin(), treeBiparts.end(), bipartSet.begin(), bipartSet.end())){
			retSet.insert(i);
			} 
		}
	}	

  return retSet;
}

set<unsigned int> get_superset_trees(int tree){ //return trees which have all of the bipartitions of the given tree
  set<unsigned int> retSet;
  vector<unsigned int> b = ::biparttable.inverted_index.at(tree);
  set<unsigned int> treeBiparts(b.begin(), b.end());
  for(unsigned int i = 0; i < ::biparttable.NumTrees; i++){
	if(i!=tree){ //don't include itself 
		vector<unsigned int> biparts = ::biparttable.inverted_index.at(i);
	 	set<unsigned int> bipartSet(biparts.begin(), biparts.end());
		if(includes(bipartSet.begin(), bipartSet.end(), treeBiparts.begin(), treeBiparts.end())){
			retSet.insert(i);
			} 
		}
	}	

  return retSet;
}

/*
set<unsigned int> clade_search(vector<int> required, bool strict){
	set<unsigned int> trees;
	vector<unsigned int> matchingtrees;
	
	
	for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //Each bipartition in the treeset		
		if (required.size() == 0){ 
			continue;
		}
		if (::biparttable.are_ones(i, required) ) { 
			if(strict){
				if (::biparttable.number_of_ones(i) != required.size()){
					continue;
				}
			}
			matchingtrees = ::biparttable.get_trees(i);
		}
		trees.insert(matchingtrees.begin(), matchingtrees.end());
	}
	return trees;
}

set<unsigned int> clade_search(vector<string> RequiredTaxa, bool strict) {
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return clade_search(required, strict);
}
*/
set<unsigned int> clade_size_search(vector<int> required, int size){
	//out << "entering clade_size_search" << endl;
	set<unsigned int> trees;
	vector<unsigned int> matchingtrees;
	
	for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //Each bipartition in the treeset		
		if (required.size() == 0){ 
			continue;
		}
		if (::biparttable.are_ones(i, required) ) { 
			if (::biparttable.number_of_ones(i) != size){
				continue;
			}
			matchingtrees = ::biparttable.get_trees(i);
		}
		trees.insert(matchingtrees.begin(), matchingtrees.end());
	}
	return trees;
}

set<unsigned int> clade_size_search(vector<string> RequiredTaxa, int size) {
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return clade_size_search(required, size);
}


set<unsigned int> smallest_clade_search(vector<int> required){
	set<unsigned int> trees;
	int smallest = ::biparttable.lm.size();
	vector<int> matchingbiparts;
	
	vector<unsigned int> matchingtrees;
	
	for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //Each bipartition in the treeset		
		if (required.size() == 0){ 
			continue;
		}
		if (::biparttable.are_ones(i, required) ) { 
			if (::biparttable.number_of_ones(i) < smallest){
				smallest = ::biparttable.number_of_ones(i);
				matchingbiparts.clear();
				matchingbiparts.push_back(i);
			}
			if (::biparttable.number_of_ones(i) == smallest){
				matchingbiparts.push_back(i);
			}
		}
	}
	
	for (unsigned int i = 0; i < matchingbiparts.size(); i++){ //Each bipartition which defines a min sized clade for the taxa	
		matchingtrees = ::biparttable.get_trees(matchingbiparts[i]);
		trees.insert(matchingtrees.begin(), matchingtrees.end());
	}
	cout << "The smallest clade containing the input taxa is of size: " << smallest << endl;
	return trees;
}

set<unsigned int> smallest_clade_search(vector<string> RequiredTaxa) {
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return smallest_clade_search(required);
}



set<unsigned int> halfbi_size_search(vector<int> required, int size)
{
	set<unsigned int> retSet;
	struct goodBiparts {
		vector<int> goodBipartitions;
		vector<bool> heteroMode;
		vector<int> size;	
	}; 
	goodBiparts B;


	for(unsigned int b = 0; b < ::biparttable.biparttable_size(); b++){ //for each bipartition
		
		//cout << "Bipartition # " << b << ", ";
		//for(int k = 0; k < biparttable.length_of_bitstrings[b]; k++){
		//	cout << biparttable.bipartitions[b][k];
		//} cout << endl;
		
		
		//create a partial bitstring out of the bipartition that represents all of the places of required	
		//bool partialBS[required.size()];
		boost::dynamic_bitset<> partialBS(required.size());	
			
		for(unsigned int i = 0; i < required.size(); i++){ //fill the partial bitstring
			//check that required is in the bitstring
			if(required.at(i) > ::biparttable.lm.size()-1){
					cerr << "Error: element " << required.at(i) << " is not a valid taxon!" << endl;
					return retSet;
				}				
			if(biparttable.bitstring_size(b) <= i){ //we have a 0
					partialBS[i] = 0;
				}			
			else{
					partialBS[i] = biparttable.get_bit(b,required.at(i));
				}
			}
		
		//cout << "Partial BS is: "; 
		for(unsigned int i = 0; i < partialBS.size(); i++) { 
			cout << partialBS[i]; 
		} 
		cout << endl; 		

		//now, check that the partial bitstring is either all 1s or all 0s (meaning its in the clade)
		if(partialBS.count() == partialBS.size() or partialBS.count() == 0){
			//now we know that all of the taxa in required are in a clade. Find out if the clade is the 1s or the 0s
			bool cladeBit = partialBS[0];
			//now, count the number of members on that side of the bipartition

			/*note- in the case of hetero data, if the bipart clade contains 0s, this presents problems
			This is because the trees which contain that bipartition may or may not have the taxa which are 0 in the bipartition
			If cladeBit is 1, on the other hand, we know that all trees containing that bipartition have all taxa in the clade 
			*/
			int bSize;			
			if(cladeBit){
				bSize = biparttable.number_of_ones(b);
				}
			else{
				//we're counting 0s, but trailing zeros are cut off, so we have to account for that since number_of_zeros doesn't
				bSize =  ::biparttable.number_of_zeros(b) + (::biparttable.lm.size() - biparttable.bitstring_size(b)); 
				}	
			//cout << "Clade found, bSize is: " << bSize << ", Bipart index is " << b << endl;						
			if(::biparttable.hetero && !cladeBit){
				//if the data is hetero and the clade is all 0s, then we don't truly know the size without looking at trees
				//the real size could be smaller since some 0s don't exist in individual trees
				if(bSize >= size){
					B.goodBipartitions.push_back(b);
					B.heteroMode.push_back(true);
					B.size.push_back(bSize);
					//cout << "Good bipartition added for hetero data! Bipart index is " << b <<", bSize is " << bSize << endl;
				}
			}
			else{
				if(bSize == size){
					B.goodBipartitions.push_back(b);
					B.heteroMode.push_back(false);
					B.size.push_back(bSize);
					//cout << "Good bipartition added! Bipart index is " << b <<", bSize is " << bSize << endl;
					}
			}
			}		
		}	


	//now we have a list of bipartitions which satisfy the desired clade. Add all trees that contain these to the return set
	for(int i = 0; i < B.goodBipartitions.size(); i++){ //for each bipartition whose trees we want
		//cout << "Good Bipartition: " << B.goodBipartitions[i] << endl;
		//cout << "Search table is: "; for(int q = 0; q < biparttable.searchtable[B.goodBipartitions[i]].size(); q++) { cout << biparttable.searchtable[B.goodBipartitions[i]][q];} cout << endl;
		
		for(int j = 0; j < biparttable.trees_size(B.goodBipartitions[i]); j++){ //for each tree in the searchtable
			if(B.heteroMode.at(i)){ //this means the clade we found contained all 0s and the data set is hetero
				//we need to find the real bSize for each tree.
				int bSize = B.size.at(i); //the size of the clade for this bipartition. 
				int tree = biparttable.get_tree(B.goodBipartitions[i],j);
				bSize = bSize + ::biparttable.num_taxa_in_tree(tree) - ::biparttable.lm.size(); //since the clade was all 0s, n of those 0s don't actually exist in the tree, NUM_TAXA - num_taxa_in_tree(tree)
				//cout << "Tree " << tree << ", adjusted bSize is " << bSize << endl;
				if(bSize == size && biparttable.are_taxa_in_tree(tree, required)){
					retSet.insert(biparttable.get_trees(B.goodBipartitions[i])[j]);
					}
				}
			else{			
				retSet.insert(biparttable.get_tree(B.goodBipartitions[i],j));
				//cout << biparttable.searchtable[B.goodBipartitions[i]][j] << " inserted!" << endl;
				}			
			}
		}
	return retSet;
}

//wrapper function for if we're given a set of strinsg
set<unsigned int> halfbi_size_search(vector<string> RequiredTaxa, int size) {
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return halfbi_size_search(required, size);
}

set<unsigned int> smallest_halfbi(vector<int> required){


	set<unsigned int> retSet;
	int smallest = biparttable.biparttable_size();

	struct goodBiparts {
		vector<int> goodBipartitions;
		vector<bool> allZeros;
		vector<int> size;	
	}; 
	goodBiparts B;	
	
	int index = 0; //an index to where in the vector is. This vector will act in a circular way. Every time we find a smaller clade, we reset the index to 0

	for(int b = 0; b < biparttable.biparttable_size(); b++){ //for each bipartition
		boost::dynamic_bitset<> partialBS(required.size());
		for(int i = 0; i < required.size(); i++){ //fill the partial bitstring
			if(required.at(i) > ::biparttable.lm.size()-1){
				cerr << "Error: element " << required.at(i) << " is not a valid taxon!" << endl;
				return retSet;
				}
			if(biparttable.bitstring_size(b) <= i){ //we have a 0
				partialBS[i] = 0;
				}			
			else{
				partialBS[i] = biparttable.get_bit(b,required.at(i));
				}
			}
		if(partialBS.count() == partialBS.size() or partialBS.count() == 0){
			bool cladeBit = partialBS[0];
			int bSize;			
			if(cladeBit){
				bSize = ::biparttable.number_of_ones(b);
				}
			else{
				//we're counting 0s, but trailing zeros are cut off, so we have to account for that since number_of_zeros doesn't
				bSize =  ::biparttable.number_of_zeros(b) + (::biparttable.lm.size() - biparttable.bitstring_size(b)); 
				}		

			if(::biparttable.hetero){
				B.goodBipartitions.push_back(b);
				B.allZeros.push_back(!cladeBit);
				B.size.push_back(bSize);
				}
			else{
				if(bSize < smallest){
					//cout << "new smallest clade size found of size " << bSize << "!" << endl;
					smallest = bSize; 
					index = 0;
					if(B.goodBipartitions.size()<=index){
						B.goodBipartitions.push_back(b);
						}
					else{ 
						B.goodBipartitions.at(index) = b;}
					index++;
					}			
				else if(bSize==smallest){
					cout << "Bipartition matches smallest clade, inserting at index " << index << endl;
					if(B.goodBipartitions.size()<=index){
						B.goodBipartitions.push_back(b);
						}
					else{ 
						B.goodBipartitions.at(index) = b;}
					index++;
					}
			}
			}		
		}	
	if(!::biparttable.hetero) {
		cout << endl << "The smallest clade found is of size " << smallest << ", found in trees:";
		for(int i = 0; i < index; i++){ //for each bipartition whose trees we want
			for(int j = 0; j < biparttable.trees_size(B.goodBipartitions[i]); j++){ //for each tree in the searchtable
				retSet.insert(biparttable.get_tree(B.goodBipartitions[i],j));
				}
			}
		}
	else{
		//we don't actually know what the smallest clade is, we only know all of the bipartitions which match the clade
		vector<int> returnTrees; //the vector of trees that we are going to return with the smallest size
		int returnTreesIndex = 0;
		//smallest is still set to the number of taxa. We'll use it to keep track of where in returnTrees we are
		for(int i = 0; i < B.goodBipartitions.size(); i++){ //for each good bipartition
			for(int j = 0; j < biparttable.trees_size(B.goodBipartitions[i]); j++){ //for each tree with that bipartition
				int tree = biparttable.get_tree(B.goodBipartitions[i],j);				
				int bSize = B.size.at(i);
				if(B.allZeros.at(i)){ //if the clade is all 0s, we have to make sure all the required biparts are in it and adjust the size
					bSize = bSize + ::biparttable.num_taxa_in_tree(tree) - ::biparttable.lm.size();
					//also check that all the required taxa exist in this tree
					if(!::biparttable.are_taxa_in_tree(tree, required)){cout << "Taxa not in tree " << tree << endl; bSize = ::biparttable.lm.size() + 1;} //make it impossible for bSize to be smallest 
					}
				//now we know the clade size (i.e. bSize) for this tree. We have to make sure bSize isn't 0 or -1, which would mean the required taxa aren't in the tree
				if(bSize==0 || bSize==-1) {bSize = ::biparttable.lm.size() + 1;}
				if(bSize == smallest){
					if(returnTrees.size() <= returnTreesIndex){
						returnTrees.push_back(tree);
						}
					else{
						returnTrees.at(returnTreesIndex) = tree;
						}
					returnTreesIndex++;
					}
				else if(bSize < smallest){
					smallest = bSize;
					returnTreesIndex = 0;
					if(returnTrees.size() <= returnTreesIndex){
						returnTrees.push_back(tree);
						}
					else{
						returnTrees.at(returnTreesIndex) = tree;
						}
					returnTreesIndex++;
					}
				}
			}
		//now, we have all of the trees we need in returnTrees, and an index of how many are the smallest
		cout << "Hetero Data: Smallest clade size is " << smallest << endl;
		for(int x = 0; x < returnTreesIndex; x++){
			retSet.insert(returnTrees.at(x));
			}
		}
			
	return retSet;

}


set<unsigned int> smallest_halfbi(vector<string> RequiredTaxa){
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return smallest_halfbi(required);
}

set<unsigned int> get_trees_with_taxa(vector<int> required){
    set<unsigned int> trees;
	//set<int>::iterator it;

    //for (unsigned int j = 0; j < required.size(); j++) {
	//	cout << required[j] << endl;
	//}
	
    for (unsigned int i = 0; i < ::biparttable.taxa_in_trees.size(); i++) { //for each tree
		int flag = 0;

        for (unsigned int j = 0; j < required.size(); j++) { //for each search taxa
			if(::biparttable.taxa_in_trees[i][required[j]] == 0) {
				flag = 1;
				break;
			}
		} 
		
		if (flag == 0) {
			trees.insert(i);
		}
    }

    return trees;
}

//Takes an input tree and returns the trees(s) that are the most similar to the input tree based on number of shared bipartitions
set<unsigned int> similarity_search(string inputtree){
	set<unsigned int> simTrees; // Return set
	unsigned int highest_sim = 0;// Holds the largest number of shared bipartitions found
	vector < vector < int > > bipartitions;
	::map<unsigned int, unsigned int> temp;

	bipartitions = compute_bitstrings_h(inputtree); // Converts the input tree from a string to bitstring representation

	for (unsigned int i = 0; i < bipartitions.size(); i++){//For each bipartition
		for (unsigned int j = 0; j <bipartitions[i].size(); j++) {//For each tree
			unsigned int ctree = bipartitions[i][j];//The tree at this current bipartition
			//cout << "The trees are: " << ctree << endl;
			temp[ctree]++; //Increments the value of the map corresponding to the tree in question
			if(temp[ctree] > highest_sim){
				highest_sim = temp[ctree];
			}
		}
	}
	//Finds the keys from the map which have values equal to the highest similarity value
	for (map<unsigned int,unsigned int>::iterator j = temp.begin(); j != temp.end(); j++){
		if(j->second == highest_sim){
		simTrees.insert(j->first);	
		}
	}
	float num_bipartitions = bipartitions.size();
	float match_percent = highest_sim /  num_bipartitions; //The percentage of similarity

	cout << "The maximum similarity is " << match_percent << "%" << endl;
	cout << "The trees with this similarity are : ";
	
	//Prints the trees
	//for(set<unsigned int>::iterator it = simTrees.begin(); it != simTrees.end(); it++){
	//		cout << *it << " " ;
	//	}
	return simTrees;
}

set<unsigned int> get_trees_without_taxa(vector<int> excluded){
    set<unsigned int> trees;
    //set< int>::iterator it;
	
    for ( unsigned int i = 0; i < ::biparttable.NumTrees; i++) {
		int flag = 0;
		
		for ( unsigned int j = 0; j < excluded.size(); j++) {
			if(::biparttable.taxa_in_trees[i][excluded[j]] == 1) {
				flag = 1;
				break;
			}
		} 

		if (flag == 0) {
			trees.insert(i);
		}
    }
    //it = unique (trees.begin(), trees.end()); 
	//trees.resize( it - trees.begin() );  
    return trees;
}

set<unsigned int> get_trees_without_taxa(vector<string> ExcludedTaxa){
	vector<int> excluded = ::biparttable.lm.lookUpLabels(ExcludedTaxa);
	return get_trees_without_taxa(excluded); 
}

set<unsigned int> get_trees_with_taxa(vector<string> RequiredTaxa) {
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return get_trees_with_taxa(required);
}

set<unsigned int> get_trees_by_taxa(vector<string> RequiredTaxa, vector<string> ExcludedTaxa) {
    set<unsigned int> trees;

	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	vector<int> excluded = ::biparttable.lm.lookUpLabels(ExcludedTaxa);
	
    for (unsigned int i = 0; i < ::biparttable.NumTrees; i++) {
		int flag = 0;

        for (unsigned int j = 0; j < required.size(); j++) {
			if(::biparttable.taxa_in_trees[i][required[j]] == 0) {
				flag = 1;
				break;
			}
		} 
		
		for (unsigned int j = 0; j < excluded.size(); j++) {
			if(::biparttable.taxa_in_trees[i][excluded[j]] == 1) {
				flag = 1;
				break;
			}
		} 

		if (flag == 0){
			trees.insert(i);
		}
    }

    return trees;
}


//A helper function to clean up redundant code in calling dfs_compute_bitstrings
//simply takes as input the string to be converted to bitstring representation
vector< vector <int > > compute_bitstrings_h(string inputstring) {
	NEWICKTREE *newickTree;
	int err;
	char * cs = strdup(inputstring.c_str());
	newickTree = stringloadnewicktree(cs, &err);
	vector< vector < int > > treeout; //Return variable

	//Error handling for stringloadnewicktree
	if (!newickTree) {
		switch (err) {
			case -1:
				cout << "Out of memory" << endl;
				break;
			case -2:
				cout << "Parse error" << endl;
				break;
			case -3:
				cout << "Can't load file" << endl;
				break;
			default:
				cout << "Error " << err << endl;
			exit(0);
			}
	}
	else {
		//Compute the bitstring representation and save it to treeout
		bool * bs = dfs_compute_bitstrings(newickTree->root, NULL, treeout);
		delete[] bs;
	}

	free(cs); 

	return treeout;
}

bool * dfs_compute_bitstrings(NEWICKNODE* startNode, NEWICKNODE* parent, vector< vector < int > > &solution ){
  //~ if (HETERO && !HCHECK){
      //~ cout << "if (HETERO && !HCHECK)" << endl;
      //~ return NULL;
  //~ }

  vector <int> subtree;

  // If the node is leaf node, just set the place of the taxon name in the bit string to '1' 
  //bitstrings here are not used to compute hash values -- however, we do collect them.
  startNode->myparent = parent;
  ///fprintf(stderr, "AT NODE: ");
  if (startNode->Nchildren == 0) { // leaf node
    string temp(startNode->label);
    unsigned int idx = ::biparttable.lm.position(temp);
    bool * bs = new bool[::biparttable.lm.size()];
	
    for (unsigned int i = 0; i < ::biparttable.lm.size(); i++){
      bs[i] = 0;
	}
    bs[idx] = true;
	
	//Error Checking
    int numOnes = 0;
    for (unsigned int i = 0; i < ::biparttable.lm.size(); i++){
      numOnes+=bs[i];
	  if (bs[i] ==  true){
        subtree.push_back(i);
      }
	}
    if (numOnes > 1 || numOnes == 0) {
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring creation error\n";
      exit(0);
    }
	//End Error Checking

	//print_bitstring(bs, ::NUM_TAXA);
	solution.push_back(subtree);
	//cout << "Here we go" << endl;
	//cout << solution.size() << endl;
	//print_vector_of_bs(solution, ::NUM_TAXA);
    return bs;
  }

  else {
    bool * ebs[startNode->Nchildren];
    //fprintf(stderr, "I am an internal node!\n");
    for (int i=0; i<startNode->Nchildren; ++i) {
      ebs[i] = dfs_compute_bitstrings(startNode->child[i], startNode, solution);
    }
    
    bool * bs = new bool[::biparttable.lm.size()];
    for (unsigned int i  = 0; i < ::biparttable.lm.size(); i++)
      bs[i] = 0;

    //fprintf(stderr, "the bitstrings I'm going to be comparing are:\n");
    /*for (int i=0; i<startNode->Nchildren; ++i) {
      fprintf(stderr, "%s\n", (startNode->child[i])->label);
      }
      fprintf(stderr, "\n");*/

    for (int i=0; i<startNode->Nchildren; ++i) {
      //cerr << "bitstring: " << *(ebs[i]) << endl;
	  if (ebs[i]) {
	    for (unsigned int j = 0; j < ::biparttable.lm.size(); j++)
	      bs[j] |= ebs[i][j];
	    delete [] ebs[i];
	    ebs[i] = NULL;
	  }
	  else {
	    if (!::biparttable.hetero) {
	      cout << "ERROR: null bitstring\n";
	      exit(0);
        }
	  }
    }

   int numOnes = 0;
   for (unsigned int i= 0; i < ::biparttable.lm.size(); i++) {
     numOnes+=bs[i];
	 if (bs[i] ==  true) {
       subtree.push_back(i);
     }
   }
    if (numOnes < 1) {
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring OR error\n";
      exit(0);
    }

	//print_bitstring(bs, ::NUM_TAXA);
	solution.push_back(subtree);
	//cout << "Here we go" << endl;
	//cout << solution.size() << endl;
	//print_vector_of_bs(solution, ::NUM_TAXA);
    return bs;
  } //end else
}



