#include "SearchFunctions.h"

set<unsigned int> get_trees_with_taxa(vector<int> required){
    set<unsigned int> trees;
	//set<int>::iterator it;

    //for (unsigned int j = 0; j < required.size(); j++) {
	//	cout << required[j] << endl;
	//}
	
    for (unsigned int i = 0; i < ::NUM_TREES; i++) { //for each tree
		int flag = 0;

        for (unsigned int j = 0; j < required.size(); j++) { //for each search taxa
			if(::taxa_in_trees[i][required[j]] == 0) {
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

set<unsigned int> get_trees_without_taxa(vector<int> excluded){
    set<unsigned int> trees;
    //set< int>::iterator it;
	
    for ( unsigned int i = 0; i < ::NUM_TREES; i++) {
		int flag = 0;
		
		for ( unsigned int j = 0; j < excluded.size(); j++) {
			if(::taxa_in_trees[i][excluded[j]] == 1) {
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
	vector<int> excluded = ::lm.lookUpLabels(ExcludedTaxa);
	return get_trees_without_taxa(excluded); 
}

set<unsigned int> get_trees_with_taxa(vector<string> RequiredTaxa) {
	vector<int> required = ::lm.lookUpLabels(RequiredTaxa);
	return get_trees_with_taxa(required);
}

set<unsigned int> get_trees_by_taxa(vector<string> RequiredTaxa, vector<string> ExcludedTaxa) {
    set<unsigned int> trees;

	vector<int> required = ::lm.lookUpLabels(RequiredTaxa);
	vector<int> excluded = ::lm.lookUpLabels(ExcludedTaxa);
	
    for (unsigned int i = 0; i < ::NUM_TREES; i++) {
		int flag = 0;

        for (unsigned int j = 0; j < required.size(); j++) {
			if(::taxa_in_trees[i][required[j]] == 0) {
				flag = 1;
				break;
			}
		} 
		
		for (unsigned int j = 0; j < excluded.size(); j++) {
			if(::taxa_in_trees[i][excluded[j]] == 1) {
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

set<unsigned int> search_hashtable_ktets(vector < vector < int > > subtrees){ 
	set<unsigned int> total_trees( ::all_trees.begin(), ::all_trees.end()); //assume all trees have all biparts. 

	for(unsigned int i = 0; i < subtrees.size(); i++) { // for each subtree
		set<unsigned int> trees;
		set<unsigned int> temp;

	for (unsigned int j = 0; j < ::biparttable.get_size(); j++) { // for each bipartition in the hashtable
		int numones = ::biparttable.number_of_ones(j);
		bool foundFlag = true;
		for(unsigned int k = 0; k < subtrees[i].size(); k++){ // for each piece of the subtree
			if ( numones != (int)subtrees[i].size() ){	 
				foundFlag= false;
				break;
			}
		 
		 //cout << "subtrees[i][k] = " << subtrees[i][k] << endl;
         if(subtrees[i][k] >= (int)::biparttable.length_of_bitstrings[j] || ::biparttable.bipartitions[j][subtrees[i][k]] == 0 ){
		   //cout << "That taxa relationship wasn't in the tree" << endl;
           foundFlag= false;
		   break;
		 }
       }

       if (foundFlag == true){
         for (unsigned int l = 0; l < ::biparttable.searchtable[j].size(); l++){	
	       trees.insert(::biparttable.searchtable[j][l]);
	     }
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


set<unsigned int> search_hashtable_strict_old(vector<int> leftside, vector<int> rightside, int side){
  //cout << "We show the bipartitions found by the search and corresponding trees" << endl;
  //keep in mind that hashtable and hash_lengths are global variables
  bool foundFlag = true;

  set<unsigned int> trees;
  bool lvalue = 0;
  bool rvalue = 0;
  
  for (unsigned int i = 0; i < ::biparttable.get_size(); i++){ //Each bipartition in the treeset

    foundFlag = true; //assume true and check each case for a contridiction to set it as false. If it passes the tests then it is true. 
	
	// if the first taxa on the left and the first taxa on the right side of the search bipart are on the same side of the bipart
	// we know we can bail out before checking anything else. But first we need to make sure there are taxa on each side of the | or we
	// get a segfault. 

	if (leftside.size() > 0){ //if there are taxa on the left side
		if (! ::biparttable.in_bitstring(i,leftside[0])) { //if the first taxa on the left side has been truncated it's treated as a 0
			lvalue = 0;
		}
		else{
			lvalue = ::biparttable.bipartitions[i][leftside[0]]; 
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
			else if (::biparttable.bipartitions[i][leftside[j]] != lvalue){
				//cout<< "broke on if (list_bs[i][leftside[j]] != lvalue)" << endl;
				foundFlag = false;
		  		break; //Break at any mistakes. 
			}
		}
	}

	if (rightside.size() > 0){
		if (::biparttable.bipartitions[i][rightside[0]] > ::biparttable.length_of_bitstrings[i]) {
			rvalue = 0;
		}
		else{
			rvalue = ::biparttable.bipartitions[i][rightside[0]];
		}

		for (unsigned int j = 0; j < rightside.size(); j++){ // Each taxa on the left side 
			//cout << length_of_bitstrings[i] <<" " << rightside[j] << " " << rvalue << endl;

			if(::biparttable.length_of_bitstrings[i] <= rightside[j]){	
				if(rvalue == 1){	
				//	cout<< "if(length_of_bitstrings[i] < rightside[j]  && rvalue == 1)" << endl;
					foundFlag = false;
			  		break; //Break at any mistakes. 
				}
			}
			else if (::biparttable.bipartitions[i][rightside[j]] != rvalue){
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
      for (unsigned int l = 0; l < ::biparttable.searchtable[i].size(); l++){
		if (HETERO == true){
			bool taxaFlag = true;
			for (unsigned int j = 0; j < leftside.size(); j++){ // Each taxa on the left side
				if (::taxa_in_trees[::biparttable.searchtable[i][l]][leftside[j]] == 0){
					taxaFlag = false;
					break;
				}
			}
			for (unsigned int j = 0; j < rightside.size(); j++){ // Each taxa on the left side 
				if (::taxa_in_trees[::biparttable.searchtable[i][l]][rightside[j]] == 0){
					taxaFlag = false;
					break;
				}
			}
			if (taxaFlag == true){
				trees.insert(::biparttable.searchtable[i][l]);
			}
		}
		else{		
			trees.insert(::biparttable.searchtable[i][l]);
		}
      }
    }
	}
  return trees;
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
    unsigned int idx = ::lm[temp];
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
	
	
	for (unsigned int i = 0; i < ::biparttable.get_size(); i++){ //Each bipartition in the treeset
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
					matchingtrees = ::biparttable.searchtable[i];
				}
				if (side == 2 && ::is_strict_homog(i, rightside) ){//right side strict
					matchingtrees = ::biparttable.searchtable[i];
				}
				if (side == 3 && ::is_strict_homog(i, leftside) && ::is_strict_homog(i, rightside) ){//both sides strict
					matchingtrees = ::biparttable.searchtable[i];
				}
			}
			else{//neither is strict
				matchingtrees = ::biparttable.searchtable[i];
			}			
		}
		trees.insert(matchingtrees.begin(), matchingtrees.end());
		//cout << "trees.size() = " << trees.size() << endl;
	}
return trees;
}

//not safe, Infinite loop when too many taxa can be requested
vector<int> random_nums_generated(int number){
	set<int> randoms;

	while(randoms.size() != number){
		randoms.insert(rand() % ::NUM_TAXA);
	}
	vector<int> v(randoms.begin(), randoms.end());
	return v;
}

//not safe, Infinite loop when too many taxa can be requested. 
vector<int> random_nums_generated2(int number, std::uniform_int_distribution<> d,  std::minstd_rand g){
	set<int> randoms;
	
	while(randoms.size() != number){
		randoms.insert(d(g));
	}
	vector<int> v(randoms.begin(), randoms.end());
	return v;
}

int random_search2(int left, int right, int side, int iterations){

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
	std::uniform_int_distribution<> d(0, ::NUM_TAXA-1);
	std::minstd_rand g(time(NULL));
	vector<int> CulRandVect;
	vector<unsigned int> CulTrees;

	vector<int> randomVect;
	
	for (int i = 0; i < ::NUM_TAXA; ++i){
		randomVect.push_back(i);
	}
	
// This version of the random number generator fails when asked for more taxa than in the set. 
	
//	for (unsigned int i = 0; i < iterations; i++){
//		vector<int> randomVect;
//		while(randomVect.size() != left+right){
//			int i = d(g);
//			if(std::find(randomVect.begin(), randomVect.end(), i) != randomVect.end()) {
//				continue;
//			}
//			else{
//				randomVect.push_back(i);
//			}
//		}

	for (unsigned int i = 0; i < iterations; i++){
		vector<int> randomVect;

		set<unsigned int> trees;

		random_shuffle ( randomVect.begin(), randomVect.end() );

		for (int j = 0; j < randomVect.size(); j++){
			CulRandVect.push_back(randomVect[j]);
		}		
		
		std::vector<int> leftside(randomVect.begin(), randomVect.begin() + left);
		std::vector<int> rightside(randomVect.begin() + left, randomVect.begin() + left+right);
				
		//vector <int>::iterator it;
		//for (it=leftside.begin(); it!=leftside.end(); it++)
		//	cout << *it << " ";
		//cout << "| ";
		//for (it=rightside.begin(); it!=rightside.end(); it++)
		//	cout << *it << " ";
		//cout << endl;    

		
		//trees = search_hashtable_strict_and_timed(leftside, rightside, side);
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
	
	for (int i = 0; i < ::NUM_TAXA; ++i){
		randomVect.push_back(i);
	}
	
	for (unsigned int i = 0; i < iterations; i++){

		set<unsigned int> trees;

		random_shuffle ( randomVect.begin(), randomVect.end() );

		for (int j = 0; j < randomVect.size(); j++){
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
	//cout << "We show the bipartitions found by the search and corresponding trees" << endl;
	//keep in mind that hashtable and hash_lengths are global variables
	//cout << "hetero = " << hetero << endl;
	
	set<unsigned int> trees;
	vector<unsigned int> matchingtrees;
	
	for (unsigned int i = 0; i < ::biparttable.get_size(); i++){ //Each bipartition in the treeset
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
				matchingtrees = ::biparttable.searchtable[i];
			}
		::HetTime += stop_clockbp();
		}
		
		
		else{ // taxa homogenious trees, all trees have all taxa
			//cout << "taxa homogenious trees, all trees have all taxa" << endl;
			if (side > 0){//if at least one side is strict
				//cout << "side is > 0" << endl;
				if (side == 1 && ::is_strict_homog(i, leftside) ){//left side strict
					matchingtrees = ::biparttable.searchtable[i];
				}
				if (side == 2 && ::is_strict_homog(i, rightside) ){//right side strict
					matchingtrees = ::biparttable.searchtable[i];
				}
				if (side == 3 && ::is_strict_homog(i, leftside) && ::is_strict_homog(i, rightside) ){//both sides strict
					matchingtrees = ::biparttable.searchtable[i];
				}
			}
			else{//neither is strict
				matchingtrees = ::biparttable.searchtable[i];
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
