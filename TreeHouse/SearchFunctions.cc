#include "SearchFunctions.h"



set<unsigned int> clade_size_search(vector<int> required, int size)
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
		bool partialBS[required.size()];
				
		for(unsigned int i = 0; i < required.size(); i++){ //fill the partial bitstring
			//check that required is in the bitstring
			if(required.at(i) > ::NUM_TAXA-1){
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
		
		cout << "Partial BS is: "; 
		for(unsigned int i = 0; i < required.size(); i++) { 
			cout << partialBS[i]; 
		} 
		cout << endl; 		

		//now, check that the partial bitstring is either all 1s or all 0s (meaning its in the clade)
		if(areBitsSame(partialBS, required.size())){
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
				bSize =  ::biparttable.number_of_zeros(b) + (::NUM_TAXA - biparttable.bitstring_size(b)); 
				}	
			//cout << "Clade found, bSize is: " << bSize << ", Bipart index is " << b << endl;						
			if(::HETERO && !cladeBit){
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
		
		for(int j = 0; j < biparttable.num_trees(B.goodBipartitions[i]); j++){ //for each tree in the searchtable
			if(B.heteroMode.at(i)){ //this means the clade we found contained all 0s and the data set is hetero
				//we need to find the real bSize for each tree.
				int bSize = B.size.at(i); //the size of the clade for this bipartition. 
				int tree = biparttable.get_tree(B.goodBipartitions[i],j);
				bSize = bSize + ::biparttable.num_taxa_in_tree(tree) - ::NUM_TAXA; //since the clade was all 0s, n of those 0s don't actually exist in the tree, NUM_TAXA - num_taxa_in_tree(tree)
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
set<unsigned int> clade_size_search(vector<string> RequiredTaxa, int size) {
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return clade_size_search(required, size);
}

set<unsigned int> smallest_clade(vector<int> required){


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
		bool partialBS[required.size()];
		for(int i = 0; i < required.size(); i++){ //fill the partial bitstring
			if(required.at(i) > ::NUM_TAXA-1){
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
		if(areBitsSame(partialBS, required.size())){
			bool cladeBit = partialBS[0];
			int bSize;			
			if(cladeBit){
				bSize = ::biparttable.number_of_ones(b);
				}
			else{
				//we're counting 0s, but trailing zeros are cut off, so we have to account for that since number_of_zeros doesn't
				bSize =  ::biparttable.number_of_zeros(b) + (::NUM_TAXA - biparttable.bitstring_size(b)); 
				}		

			if(::HETERO){
				B.goodBipartitions.push_back(b);
				B.allZeros.push_back(!cladeBit);
				B.size.push_back(bSize);
				}
			else{
				if(bSize < smallest){
					cout << "new smallest clade size found of size " << bSize << "!" << endl;
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
	if(!::HETERO) {
		cout << endl << "The smallest clade found is of size " << smallest << ", found in trees:";
		for(int i = 0; i < index; i++){ //for each bipartition whose trees we want
			for(int j = 0; j < biparttable.num_trees(B.goodBipartitions[i]); j++){ //for each tree in the searchtable
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
			for(int j = 0; j < biparttable.num_trees(B.goodBipartitions[i]); j++){ //for each tree with that bipartition
				int tree = biparttable.get_tree(B.goodBipartitions[i],j);				
				int bSize = B.size.at(i);
				if(B.allZeros.at(i)){ //if the clade is all 0s, we have to make sure all the required biparts are in it and adjust the size
					bSize = bSize + ::biparttable.num_taxa_in_tree(tree) - ::NUM_TAXA;
					//also check that all the required taxa exist in this tree
					if(!::biparttable.are_taxa_in_tree(tree, required)){cout << "Taxa not in tree " << tree << endl; bSize = ::NUM_TAXA + 1;} //make it impossible for bSize to be smallest 
					}
				//now we know the clade size (i.e. bSize) for this tree. We have to make sure bSize isn't 0 or -1, which would mean the required taxa aren't in the tree
				if(bSize==0 || bSize==-1) {bSize = ::NUM_TAXA + 1;}
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


set<unsigned int> smallest_clade(vector<string> RequiredTaxa){
	vector<int> required = ::biparttable.lm.lookUpLabels(RequiredTaxa);
	return smallest_clade(required);
}

set<unsigned int> get_trees_with_taxa(vector<int> required){
    set<unsigned int> trees;
	//set<int>::iterator it;

    //for (unsigned int j = 0; j < required.size(); j++) {
	//	cout << required[j] << endl;
	//}
	
    for (unsigned int i = 0; i < ::NUM_TREES; i++) { //for each tree
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
	
    for ( unsigned int i = 0; i < ::NUM_TREES; i++) {
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
	
    for (unsigned int i = 0; i < ::NUM_TREES; i++) {
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

//GRB at min needs a rename
set<unsigned int> search_hashtable_ktets(vector < vector < int > > subtrees){ 
	set<unsigned int> total_trees( ::all_trees.begin(), ::all_trees.end()); //assume all trees have all biparts. 

	for(unsigned int i = 0; i < subtrees.size(); i++) { // for each subtree
		set<unsigned int> trees;
		set<unsigned int> temp;

	for (unsigned int j = 0; j < ::biparttable.biparttable_size(); j++) { // for each bipartition in the hashtable
		int numones = ::biparttable.number_of_ones(j);
		bool foundFlag = true;
		for(unsigned int k = 0; k < subtrees[i].size(); k++){ // for each piece of the subtree
			if ( numones != (int)subtrees[i].size() ){	 
				foundFlag= false;
				break;
			}
		 
		 //cout << "subtrees[i][k] = " << subtrees[i][k] << endl;
         if(subtrees[i][k] >= (int)::biparttable.bitstring_size(j) || ::biparttable.get_bit(j, subtrees[i][k]) == 0 ){
		   //cout << "That taxa relationship wasn't in the tree" << endl;
           foundFlag= false;
		   break;
		 }
       }

       if (foundFlag == true){
	     trees.insert(::biparttable.trees_begin(j), ::biparttable.trees_end(j));
       }
     }

	std::set_intersection(total_trees.begin(), total_trees.end(), trees.begin(), trees.end(),  std::inserter(temp, temp.begin()) );
		
	swap(total_trees, temp);

  }

  return total_trees;
}

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
		trees = search_hashtable_ktets(subtrees);
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

/*
set<unsigned int> search_hashtable_strict_old(vector<int> leftside, vector<int> rightside, int side){
  //cout << "We show the bipartitions found by the search and corresponding trees" << endl;
  //keep in mind that hashtable and hash_lengths are global variables
  bool foundFlag = true;

  set<unsigned int> trees;
  bool lvalue = 0;
  bool rvalue = 0;
  
  for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //Each bipartition in the treeset

    foundFlag = true; //assume true and check each case for a contridiction to set it as false. If it passes the tests then it is true. 
	
	// if the first taxa on the left and the first taxa on the right side of the search bipart are on the same side of the bipart
	// we know we can bail out before checking anything else. But first we need to make sure there are taxa on each side of the | or we
	// get a segfault. 

	if (leftside.size() > 0){ //if there are taxa on the left side
		if (! ::biparttable.in_bitstring(i,leftside[0])) { //if the first taxa on the left side has been truncated it's treated as a 0
			lvalue = 0;
		}
		else{
			lvalue = ::biparttable.get_bit(i, leftside[0]); 
		}

		for (unsigned int j = 0; j < leftside.size(); j++){  // Each taxa on the left side 
			//cout << length_of_bitstrings[i] <<" " << leftside[j] << " " << lvalue << endl;

			if(! ::biparttable.in_bitstring(i,leftside[0])){	
				if(lvalue == 1){	
					//cout<< "if(length_of_bitstrings[i] < leftside[j]  && lvalue == 1)" << endl;
					foundFlag = false;
			  		break; //Break at any mistakes. 
				}
			}
			else if (::biparttable.get_bit(i,leftside[j]) != lvalue){
				//cout<< "broke on if (list_bs[i][leftside[j]] != lvalue)" << endl;
				foundFlag = false;
		  		break; //Break at any mistakes. 
			}
		}
	}

	if (rightside.size() > 0){
		if (::biparttable.get_bit(i,rightside[0]) > ::biparttable.bitstring_size(i)) {
			rvalue = 0;
		}
		else{
			rvalue = ::biparttable.get_bit(i,rightside[0]);
		}

		for (unsigned int j = 0; j < rightside.size(); j++){ // Each taxa on the left side 
			//cout << length_of_bitstrings[i] <<" " << rightside[j] << " " << rvalue << endl;

			if(::biparttable.bitstring_size(i) <= rightside[j]){	
				if(rvalue == 1){	
				//	cout<< "if(length_of_bitstrings[i] < rightside[j]  && rvalue == 1)" << endl;
					foundFlag = false;
			  		break; //Break at any mistakes. 
				}
			}
			else if (::biparttable.get_bit(i,rightside[j]) != rvalue){
				//cout<< "broke on if (list_bs[i][rightside[j]] != rvalue)" << endl;
				foundFlag = false;
		  		break; //Break at any mistakes. 
			}
		}
	}

	if (leftside.size() > 0 && rightside.size() > 0 && lvalue == rvalue){
		//cout<< "broke on if (leftside.size() > 0 && rightside.size() > 0 && lvalue == rvalue)" << endl;
		foundFlag = false;
		continue;	
	}

	if (side > 0){
		//cout << "side is > 0" << endl;
		//unsigned int mycount = count (::list_bs[i], ::list_bs[i]+::length_of_bitstrings[i], 1);
		unsigned int mycount = ::biparttable.number_of_ones(i); 
		
		if (side == 1 || side == 3){
			//cout << foundFlag << " "  << side << " " <<  lvalue << " " << leftside.size() << " " << (NUM_TAXA-mycount) << " " << mycount << endl;
			if(lvalue == 1 && leftside.size() != mycount){
				//cout << "Hit lvalue == 1 && leftside.size() != mycount)" << endl;
				foundFlag = false;
	  			continue;
			}
			else if(lvalue == 0 && leftside.size() != (NUM_TAXA-mycount)){
				//cout << "(lvalue == 0 && leftside.size() != (NUM_TAXA-mycount))" << endl;
				foundFlag = false;
	  			continue;
			}
		}
		
		if (side == 2 || side == 3){
			if(rvalue == 1 && rightside.size() != mycount){
				//cout << "(rvalue == 1 && rightside.size() != mycount)" << endl;
				foundFlag = false;
	  			continue;
			}
			else if(rvalue == 0 && rightside.size() != (NUM_TAXA-mycount)){
				//cout << "(rvalue == 0 && rightside.size() != (NUM_TAXA-mycount))" << endl;
				foundFlag = false;
	  			continue;
			}
		}
	}
	//If it is true I don't know if 0's mean that it's on the otherside of the bipartition or not in the tree at all. So check. 

	//if we got here and it's true it passed all the checks. 
	if (foundFlag == true){
      for (unsigned int l = 0; l < ::biparttable.trees_size(i); l++){
		if (HETERO == true){
			bool taxaFlag = true;
			for (unsigned int j = 0; j < leftside.size(); j++){ // Each taxa on the left side
				if (::taxa_in_trees[::biparttable.get_tree(i,l)][leftside[j]] == 0){
					taxaFlag = false;
					break;
				}
			}
			for (unsigned int j = 0; j < rightside.size(); j++){ // Each taxa on the left side 
				if (::taxa_in_trees[::biparttable.get_tree(i,l)][rightside[j]] == 0){
					taxaFlag = false;
					break;
				}
			}
			if (taxaFlag == true){
				trees.insert(::biparttable.get_tree(i,l));
			}
		}
		else{		
			trees.insert(::biparttable.get_tree(i,l));
		}
      }
    }
	}
  return trees;
}
*/

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

//A very similar version of this function is found in HashTableSearch.cc (differences are some additional code is commented out here
//and &solution is a vector< vector <unsigned int> > in HashTableSearch.cc
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
    unsigned int idx = ::biparttable.lm[temp];
    bool * bs = new bool[::NUM_TAXA];
	
    for (unsigned int i = 0; i < ::NUM_TAXA; i++){
      bs[i] = 0;
	}
    bs[idx] = true;
	
	//Error Checking
    int numOnes = 0;
    for (unsigned int i = 0; i < ::NUM_TAXA; i++){
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
    
    bool * bs = new bool[::NUM_TAXA];
    for (unsigned int i  = 0; i < ::NUM_TAXA; i++)
      bs[i] = 0;

    //fprintf(stderr, "the bitstrings I'm going to be comparing are:\n");
    /*for (int i=0; i<startNode->Nchildren; ++i) {
      fprintf(stderr, "%s\n", (startNode->child[i])->label);
      }
      fprintf(stderr, "\n");*/

    for (int i=0; i<startNode->Nchildren; ++i) {
      //cerr << "bitstring: " << *(ebs[i]) << endl;
	  if (ebs[i]) {
	    for (unsigned int j = 0; j < ::NUM_TAXA; j++)
	      bs[j] |= ebs[i][j];
	    delete [] ebs[i];
	    ebs[i] = NULL;
	  }
	  else {
	    if (!HETERO) {
	      cout << "ERROR: null bitstring\n";
	      exit(0);
        }
	  }
    }

   int numOnes = 0;
   for (unsigned int i= 0; i < ::NUM_TAXA; i++) {
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



set<unsigned int> search_hashtable_strict(vector<int> leftside, vector<int> rightside, int side){
	//cout << "We show the bipartitions found by the search and corresponding trees" << endl;
	//keep in mind that hashtable and hash_lengths are global variables
	//cout << "hetero = " << hetero << endl;
	
	set<unsigned int> trees;
	vector<unsigned int> matchingtrees;
	
	//Testing features for 0 by, 1 by, and missnamed taxa.

	if(HETERO){
		if(leftside.size() == 0 && rightside.size() == 0){
			return ::all_trees;
		}
		if(leftside.size() <= 1){
			return get_trees_with_taxa(rightside);
		}
		if(rightside.size() <= 1){
			return get_trees_with_taxa(leftside);
		}
	}
	
	
	for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //Each bipartition in the treeset
		//foundFlag = true; //assume true and check each case for a contridiction to set it as false. If it passes the tests then it is true. 
		// if the first taxa on the left and the first taxa on the right side of the search bipart are on the same side of the bipart
		// we know we can bail out before checking anything else. But first we need to make sure there are taxa on each side of the | or we
		// get a segfault. 
		
		//cout << "Looking at bipartition " << i << endl;
		
		if (leftside.size() > 0){
			if (! ::biparttable.same_bitstring_value(i, leftside) ) { // taxa on the left side of bipartition are all on 1's or 0's
				//cout<< "broke on if (! same_bitstring_value(i, leftside) )" << endl;
				continue;
			}
		}
	
		if (rightside.size() > 0){
			if (! ::biparttable.same_bitstring_value(i, rightside) ) { // taxa on the right side of bipartition are all on 1's or 0's
				//cout<< "broke on if (! same_bitstring_value(i, rightside) )" << endl;
				continue;
			}
		}
	
		if (leftside.size() > 0 && rightside.size() > 0 && (::biparttable.is_one(i,leftside[0]) == ::biparttable.is_one(i,rightside[0]) ) ){
			//cout<< "broke on if (leftside.size() > 0 && rightside.size() > 0 && lvalue == rvalue)" << endl;
			//cout << "leftside[0] = " << leftside[0] << " which is a " << is_one(i,leftside[0]) << endl;
			//cout << "rightside[0] = " << rightside[0] << " which is a " << is_one(i,rightside[0]) << endl;
				continue;	
		}
		
		//we passed the test to match a bipartition
		//However we might only want a subset of those trees thanks to
		//taxa Heterogenious trees / bipartition truncation and Strict bipartition matching
		
		if (::HETERO == true){//if the trees are taxa heterogenious
			//cout << "if the trees are taxa heterogenious" << endl;
			if (side > 0){//if at least one side is strict
				//cout << "side is > 0" << endl;
				if (side == 1){//left side strict
					matchingtrees = ::which_strict_hetero(i, leftside);
				}
				if (side == 2){//right side strict
					matchingtrees = ::which_strict_hetero(i, rightside);		
				}
				if (side == 3){//both sides strict
					matchingtrees = ::which_strict_hetero(i, leftside, rightside);	
				}
			}
			else{//neither is strict
				if(::biparttable.is_one(i,leftside[0])){
					matchingtrees = ::which_hetero(i, rightside);
				}
				else{
					matchingtrees = ::which_hetero(i, leftside);
				}
			}			
		}
		else{ // taxa homogenious trees, all trees have all taxa
			//cout << "taxa homogenious trees, all trees have all taxa" << endl;
			if (side > 0){//if at least one side is strict
				//cout << "side is > 0" << endl;
				if (side == 1 && ::is_strict_homog(i, leftside) ){//left side strict
					matchingtrees = ::biparttable.get_trees(i);
				}
				if (side == 2 && ::is_strict_homog(i, rightside) ){//right side strict
					matchingtrees = ::biparttable.get_trees(i);
				}
				if (side == 3 && ::is_strict_homog(i, leftside) && ::is_strict_homog(i, rightside) ){//both sides strict
					matchingtrees = ::biparttable.get_trees(i);
				}
			}
			else{//neither is strict
				matchingtrees = ::biparttable.get_trees(i);
			}			
		}
		trees.insert(matchingtrees.begin(), matchingtrees.end());
		//cout << "trees.size() = " << trees.size() << endl;
	}
return trees;
}

//not safe, Infinite loop when too many taxa can be requested
/*vector<int> random_nums_generated(unsigned int number){
	set<int> randoms;

	while(randoms.size() != number){
		randoms.insert(rand() % ::NUM_TAXA);
	}
	vector<int> v(randoms.begin(), randoms.end());
	return v;
}

//not safe, Infinite loop when too many taxa can be requested. 
vector<int> random_nums_generated2(unsigned int number, std::uniform_int_distribution<> d,  std::minstd_rand g){
	set<int> randoms;
	
	while(randoms.size() != number){
		randoms.insert(d(g));
	}
	vector<int> v(randoms.begin(), randoms.end());
	return v;
}
*/
int random_search2(int left, int right, int side, int iterations){
	cout<< "entering rsearch2" <<endl;
	vector<int>::iterator it;
	set<unsigned int> trees;		
	
	::SetInsertions = 0;	
	::SetOps = 0;
	::SetTime = 0;
	::HetTime = 0;
	::AHetTime = 0;	
	int TreesFound = 0;
	int SucessSearch = 0;
	
	start_clock();
	srand ( time(NULL) );
	//std::uniform_int_distribution<> d(0, ::NUM_TAXA-1);
	//std::minstd_rand g(time(NULL));
	vector<int> CulRandVect;
	vector<unsigned int> CulTrees;

	vector<int> randomVect;
	
	for (unsigned int i = 0; i < ::NUM_TAXA; ++i){
		randomVect.push_back(i);
	}
	
	cout << "rsearch2 initialization complete" << endl;
	
	
	for (int i = 0; i < iterations; i++){
		//vector<int> randomVect;

		set<unsigned int> trees;

		random_shuffle ( randomVect.begin(), randomVect.end() );
		
		cout << "suffled" << endl;
		
		for (unsigned int j = 0; j < randomVect.size(); j++){
			CulRandVect.push_back(randomVect[j]);
		}
		
		
		std::vector<int> leftside(randomVect.begin(), randomVect.begin() + left);
		std::vector<int> rightside(randomVect.begin() + left, randomVect.begin() + left+right);
		
		cout << "left and right side created" << endl;
		
		vector <int>::iterator it;
		for (it=leftside.begin(); it!=leftside.end(); it++)
			cout << *it << " ";
		cout << "| ";
		for (it=rightside.begin(); it!=rightside.end(); it++)
			cout << *it << " ";
		cout << endl;    

		
		//trees = search_hashtable_strict_and_timed(leftside, rightside, side);
		cout << "about to run a search"<<endl;
		trees = search_hashtable_strict_and_timed(leftside, rightside, side);
		TreesFound += trees.size();
		
		
		set<unsigned int>::iterator myIterator; 
		for(myIterator = trees.begin(); myIterator != trees.end(); ++myIterator){
			CulTrees.push_back(*myIterator); 
		}
		
		if (trees.size() > 0){
			SucessSearch +=1;
		}
		
	}
	
	double sum = std::accumulate(CulRandVect.begin(), CulRandVect.end(), 0.0);
	double mean = sum / CulRandVect.size();
	double sq_sum = std::inner_product(CulRandVect.begin(), CulRandVect.end(), CulRandVect.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / CulRandVect.size() - mean * mean);
	
	cout << "Minimum element in random is: " << *( std::min_element( CulRandVect.begin(), CulRandVect.end() ) ) << endl;
	cout << "Maximum element in random is: " << *( std::max_element( CulRandVect.begin(), CulRandVect.end() ) ) << endl;
	cout << "Sum of random is: " << sum << endl;
	cout << "mean of random is: " << mean << endl;
	cout << "sq_sum of random is: " << sq_sum << endl;
	cout << "stdev of random is: " << stdev << endl;


	if (CulTrees.size() > 0){
		sum = std::accumulate(CulTrees.begin(), CulTrees.end(), 0.0);
		mean = sum / CulTrees.size();

		sq_sum = std::inner_product(CulTrees.begin(), CulTrees.end(), CulTrees.begin(), 0.0);
		stdev = std::sqrt(sq_sum / CulTrees.size() - mean * mean);
		
		cout << "Minimum element in CulTrees is: " << *( std::min_element( CulTrees.begin(), CulTrees.end() ) ) << endl;
		cout << "Maximum element in CulTrees is: " << *( std::max_element( CulTrees.begin(), CulTrees.end() ) ) << endl;
		cout << "Sum of CulTrees is: " << sum << endl;
		cout << "mean of CulTrees is: " << mean << endl;
		cout << "sq_sum of CulTrees is: " << sq_sum << endl;
		cout << "stdev of CulTrees is: " << stdev << endl;
	}

	
	cout << "ExeTime = " << stop_clockbp() << endl;
	cout << "SetTime = " << ::SetTime << endl;
	cout << "HetTime = " << ::HetTime << endl;
	cout << "AHetTime = " << ::AHetTime << endl;
	cout << "SetOps = " << ::SetOps << endl;
	cout << "SetInsertions = " << ::SetInsertions << endl;
	cout << "TreesFound = " << TreesFound << endl;
	cout << "SucessfulSearches = " << SucessSearch << endl;
	cout << "ding" << endl;
	
	return SucessSearch;
}

int random_search(int left, int right, int side, int iterations){

	vector<int> treeFreq(::NUM_TREES, 0);
	
	vector<int>::iterator it;

	::SetInsertions = 0;	
	::SetOps = 0;
	::SetTime = 0;
	::HetTime = 0;
	::AHetTime = 0;	
	int TreesFound = 0;
	int SucessSearch = 0;
	srand ( unsigned ( time (NULL) ) );
	start_clock();
	vector<int> CulRandVect;
	vector<unsigned int> CulTrees;
	
	vector<int> randomVect;
	
	for (unsigned int i = 0; i < ::NUM_TAXA; ++i){
		randomVect.push_back(i);
	}
	
	for (int i = 0; i < iterations; i++){

		set<unsigned int> trees;

		random_shuffle ( randomVect.begin(), randomVect.end() );

		for (unsigned int j = 0; j < randomVect.size(); j++){
			CulRandVect.push_back(randomVect[j]);
		}		
		
		std::vector<int> leftside(randomVect.begin(), randomVect.begin() + left);
		std::vector<int> rightside(randomVect.begin() + left, randomVect.begin() + left+right);
		
		//vector <int>::iterator it;
		//for (it=leftside.begin(); it!=leftside.end(); it++)
			//cout << *it << " ";
		//cout << "| ";
		//for (it=rightside.begin(); it!=rightside.end(); it++)
			//cout << *it << " ";
		//cout << endl;    
		
		trees = search_hashtable_strict(leftside, rightside, side);
		
		TreesFound += trees.size();
		
		
 
		for(set<unsigned int>::iterator myIterator = trees.begin(); myIterator != trees.end(); ++myIterator){
			CulTrees.push_back(*myIterator); 
		}
		
		if (trees.size() > 0){
			SucessSearch +=1;
		}
		
	}
	
	double sum = std::accumulate(CulRandVect.begin(), CulRandVect.end(), 0.0);
	double mean = sum / CulRandVect.size();
	double sq_sum = std::inner_product(CulRandVect.begin(), CulRandVect.end(), CulRandVect.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / CulRandVect.size() - mean * mean);
	
	cout << "Minimum element in random is: " << *( std::min_element( CulRandVect.begin(), CulRandVect.end() ) ) << endl;
	cout << "Maximum element in random is: " << *( std::max_element( CulRandVect.begin(), CulRandVect.end() ) ) << endl;
	cout << "Sum of random is: " << sum << endl;
	cout << "mean of random is: " << mean << endl;
	cout << "sq_sum of random is: " << sq_sum << endl;
	cout << "stdev of random is: " << stdev << endl;


	if (CulTrees.size() > 0){
		sum = std::accumulate(CulTrees.begin(), CulTrees.end(), 0.0);
		mean = sum / CulTrees.size();

		sq_sum = std::inner_product(CulTrees.begin(), CulTrees.end(), CulTrees.begin(), 0.0);
		stdev = std::sqrt(sq_sum / CulTrees.size() - mean * mean);
		
		cout << "Minimum element in CulTrees is: " << *( std::min_element( CulTrees.begin(), CulTrees.end() ) ) << endl;
		cout << "Maximum element in CulTrees is: " << *( std::max_element( CulTrees.begin(), CulTrees.end() ) ) << endl;
		cout << "Sum of CulTrees is: " << sum << endl;
		cout << "mean of CulTrees is: " << mean << endl;
		cout << "sq_sum of CulTrees is: " << sq_sum << endl;
		cout << "stdev of CulTrees is: " << stdev << endl;
	}

	
	cout << "ExeTime = " << stop_clockbp() << endl;
	cout << "SetTime = " << ::SetTime << endl;
	cout << "HetTime = " << ::HetTime << endl;
	cout << "AHetTime = " << ::AHetTime << endl;
	cout << "SetOps = " << ::SetOps << endl;
	cout << "SetInsertions = " << ::SetInsertions << endl;
	cout << "TreesFound = " << TreesFound << endl;
	cout << "SucessfulSearches = " << SucessSearch << endl;
	cout << "ding" << endl;
	
	return SucessSearch;
}




set<unsigned int> search_hashtable_strict_and_timed(vector<int> leftside, vector<int> rightside, int side){
	cout << "We show the bipartitions found by the search and corresponding trees" << endl;
	//keep in mind that hashtable and hash_lengths are global variables
	cout << "hetero = " << ::HETERO << endl;
	
	set<unsigned int> trees;
	vector<unsigned int> matchingtrees;
	
	for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //Each bipartition in the treeset
		//foundFlag = true; //assume true and check each case for a contridiction to set it as false. If it passes the tests then it is true. 
		// if the first taxa on the left and the first taxa on the right side of the search bipart are on the same side of the bipart
		// we know we can bail out before checking anything else. But first we need to make sure there are taxa on each side of the | or we
		// get a segfault. 
		
		//cout << "Looking at bipartition " << i << endl;
		
		if (leftside.size() > 0){
			if (! ::biparttable.same_bitstring_value(i, leftside) ) { // taxa on the left side of bipartition are all on 1's or 0's
				//cout<< "broke on if (! same_bitstring_value(i, leftside) )" << endl;
				continue;
			}
		}
		
		if (rightside.size() > 0){
			if (! ::biparttable.same_bitstring_value(i, rightside) ) { // taxa on the right side of bipartition are all on 1's or 0's
				//cout<< "broke on if (! same_bitstring_value(i, rightside) )" << endl;
				continue;
			}
		}
		
		if (leftside.size() > 0 ){
			if(rightside.size() > 0 ) {
				if(::biparttable.is_one(i, leftside[0]) == ::biparttable.is_one(i, rightside[0]) ){
			//cout<< "broke on if (leftside.size() > 0 && rightside.size() > 0 && lvalue == rvalue)" << endl;
			//cout << "leftside[0] = " << leftside[0] << " which is a " << is_one(i,leftside[0]) << endl;
			//cout << "rightside[0] = " << rightside[0] << " which is a " << is_one(i,rightside[0]) << endl;	
					continue;	
				}
			}
		}
	
		//we passed the test to match a bipartition
		//However we might only want a subset of those trees thanks to
		//taxa Heterogenious trees / bipartition truncation and Strict bipartition matching
		
		//this block returns the trees with the correct number of taxa. The actual taxa are checked later
		if (::HETERO == true){//if the trees are taxa heterogenious
			start_clock();
			//cout << "if the trees are taxa heterogenious" << endl;
			if (side > 0){//if at least one side is strict
				//cout << "side is > 0" << endl;
				if (side == 1){//left side strict
					matchingtrees = ::which_trees_single_strict(i, leftside);
				}
				if (side == 2){//right side strict
					matchingtrees = ::which_trees_single_strict(i, rightside);
				}
				if (side == 3){//both sides strict 
					matchingtrees = ::which_trees_double_strict(i, leftside.size()+rightside.size());
				}
			}
			else{//neither is strict
				matchingtrees = ::biparttable.get_trees(i);
			}
		::HetTime += stop_clockbp();
		}
		
		
		else{ // taxa homogenious trees, all trees have all taxa
			//cout << "taxa homogenious trees, all trees have all taxa" << endl;
			if (side > 0){//if at least one side is strict
				//cout << "side is > 0" << endl;
				if (side == 1 && ::is_strict_homog(i, leftside) ){//left side strict
					matchingtrees = ::biparttable.get_trees(i);
				}
				if (side == 2 && ::is_strict_homog(i, rightside) ){//right side strict
					matchingtrees = ::biparttable.get_trees(i);
				}
				if (side == 3 && ::is_strict_homog(i, leftside) && ::is_strict_homog(i, rightside) ){//both sides strict
					matchingtrees = ::biparttable.get_trees(i);
				}
			}
			else{//neither is strict
				matchingtrees = ::biparttable.get_trees(i);
			}
		}
				
		::SetOps += 1;
		::SetInsertions += matchingtrees.size();
		start_clock();
		trees.insert(matchingtrees.begin(), matchingtrees.end());
		::SetTime += stop_clockbp();
		//cout << "trees.size() = " << trees.size() << endl;
		
		if (trees.size() == ::NUM_TREES){
			break;
		}
		
	}
	
	if (::HETERO == true){ // get the set of trees that contain all searched taxa and do a set difference with matches. 
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
		::AHetTime += stop_clockbp();
	}
	
	return trees;
}
