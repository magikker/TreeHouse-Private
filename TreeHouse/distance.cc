#include "distance.h"
using namespace std;
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
 string taxon1name = ::biparttable.lm.name(taxon1);
  string taxon2name = ::biparttable.lm.name(taxon2);
	//cout << "taxa names are: " << taxon1name << " and " << taxon2name << endl;

  vector< Bipartition > tree_bipartitions = get_tree_bipartitions(tree);
  vector<unsigned int> tree_bs_sizes = get_tree_bs_sizes(tree);
  int distance = 0;
  for (unsigned int i=0; i<tree_bipartitions.size(); i++) {
    if (!tree_bipartitions[i].same_bitstring_value(taxon1, taxon2))
      distance++;
  }
  return distance;
}

double average_distance_between_taxa(unsigned int taxon1, unsigned int taxon2){
  //only homogeneous?

  if(!::biparttable.hetero){
//  cout << "we have a homogeneous data set!" << endl;
  //takes the average distance between taxa for every tree.
  //first, we need to look at every bipartition and indicate whether the taxa are on a different side of the edge
  //each edge where the taxa differ adds 1 to the distance
  
  unsigned int numBiparts = biparttable.BipartTable.size();
  unsigned long int total = 0;

  for(unsigned int i = 0; i < numBiparts; i++){
	if(!biparttable.BipartTable[i].same_bitstring_value(taxon1, taxon2)){
 	 	//cout << "tree size for bipartition " << i << " is " << biparttable.BipartitionTable[i].trees_size() << endl;
		total+=biparttable.BipartTable[i].trees_size();
		}
	}

/*
 //for testing this function against Mark Adamo's distance_between_taxa function
 unsigned int total2 = 0;
  for(int i = 0; i < ::NUM_TREES; i++)
 {
   total2+=distance_between_taxa(taxon1, taxon2, i);
	}
  cout << "total is: " << total << ", total 2 is " << total2 << endl;
*/

  return total/(double)::biparttable.NumTrees;
  
 }
  else{
  	cout << "We have a heterogeneous data set! Return 0 for now!" << endl;
  	return 0;
	}

return 0;

}

double average_ancestral_distance(unsigned int taxon1, unsigned int taxon2){
  set<unsigned int> allTrees;
  for(unsigned int i = 0; i < ::biparttable.NumTrees; i++){
	allTrees.insert(i);
	}
  return average_ancestral_distance(taxon1, taxon2, allTrees);
}

double average_ancestral_distance(unsigned int taxon1, unsigned int taxon2, set<unsigned int> treeSet){
//note- this uses a bruteforce approach
//because we have to look at rooted trees, we need to get the newick strings for all trees
//for that reason, we can't take any shortcuts by just looking at the bipartition table
  unsigned int total = 0;
  for(unsigned int i = 0; i < treeSet.size(); i++){
	total+=distance_to_common_ancestor(taxon1, taxon2, i);
	}

  return total/(double)treeSet.size();

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
  string nwTree = to_newick(tree);
  //cout << "newick string is: " << nwTree << endl;
  //cout << "printing clade_parsed..." << endl;
  //printVector(clade_parsed);
  int depth = 0;
  unsigned int i = 0;
  //count depth from left side of string
  for (; i<clade_parsed.size(); i++) {
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
  for (; i<clade_parsed.size(); i++) {
    if (clade_parsed[i] == taxon1name)
      break;
  }
  depth = 0;
  //count depth to right side of string
  for (; i<clade_parsed.size(); i++) {
    if (clade_parsed[i] == "(")
      depth--;
    else if (clade_parsed[i] == ")")
      depth++;
  }
  return depth;
}



unsigned int distance_to_root(unsigned int taxon, unsigned int tree){
  string taxonName = ::biparttable.lm.name(taxon);
  cout << "taxon1's name is: " << taxonName << endl;
  string nwTree = to_newick(tree);
  cout << "newick string is: " << nwTree << endl;
  //now, parse the newick string and get the depth of the desired taxon
  unsigned int depth = 0;

  for(unsigned int i = 0; i < nwTree.length(); i++){
	if(nwTree.at(i)=='('){
	++depth;
		}
	else if(nwTree.at(i)==')'){
	--depth;
		}	
	else if(nwTree.at(i)==','){
		}
	else{
		//we found a taxon. Turn it into a string and compare with taxon1name
		unsigned int toJump = 0;
		unsigned int room = nwTree.length()-i-1;
		string t; t.push_back(nwTree.at(i));
		while(room>0){
			char next = nwTree.at(i+1+toJump);
			if(next=='(' || next==',' || next==')'){
				//we're at the end of the taxon name. break out of while loop
				room = 0;
				}
			else{
				t.push_back(next);
				room--;
				toJump++;
				}
			}
		//now, we have the taxon name in string t. Compare it to taxonName
		if(taxonName.compare(t)==0){ //the strings are equal. We found a match!
			return depth;
			}
		else{
			i+=toJump;
			}
		}
	}			

  cout << "distance_to_root error: taxon not found!";
  return 0;

}

double average_depth(unsigned int tree){ //returns the average taxon depth of a tree.

  string nwTree = to_newick(tree);
  //cout << "newick string is: " << nwTree << endl;
  //now, parse the newick string and get the depth of the desired taxon
  unsigned int depth = 0;
  unsigned int totalDepth = 0;


  for(unsigned int i = 0; i < nwTree.length(); i++){
	if(nwTree.at(i)=='('){
	++depth;
		}
	else if(nwTree.at(i)==')'){
	--depth;
		}	
	else if(nwTree.at(i)==',' || nwTree.at(i)==';'){
		}
	else{ 
		unsigned int room = nwTree.length()-i-1;
		unsigned int toJump = 0;		
		while(room>0){
			char next = nwTree.at(i+1+toJump);
			if(next=='(' || next==')' || next==','){
				break;
				}
			else{toJump++; room--;}
			}
		totalDepth+=depth;
		//cout << "Done with taxa, total depth is " << totalDepth << ", string index is " << i << endl;
		}
	}

return totalDepth/(double)biparttable.num_taxa_in_tree(tree);

}

double expected_average_depth(unsigned int n){
  //according to Kirkpatrick and Slatkin (1992), expected average depth is  2 * summation(i = 2 to n) of (1/i)
  double ret = 0;
  if(n<2){
	return 0;
	}
  for(int i = 2; i <= n; i++){
	ret += (double)1/(double)i;
	}
  return 2*ret;

}

double depth_variance(unsigned int tree){

  string nwTree = to_newick(tree);
  unsigned int depth = 0;
  vector<unsigned int> treeDepths;

  for(unsigned int i = 0; i < nwTree.length(); i++){
	if(nwTree.at(i)=='('){
	++depth;
		}
	else if(nwTree.at(i)==')'){
	--depth;
		}	
	else if(nwTree.at(i)==',' || nwTree.at(i)==';'){
		}
	else{ 
		unsigned int room = nwTree.length()-i-1;
		unsigned int toJump = 0;		
		while(room>0){
			char next = nwTree.at(i+1+toJump);
			if(next=='(' || next==')' || next==','){
				break;
				}
			else{toJump++; room--;}
			}
		treeDepths.push_back(depth);
		}
	}

  double total = accumulate(treeDepths.begin(), treeDepths.end(), 0.0);
  double ave = total / (double)treeDepths.size();
  double sumOfSquares = 0;

//now we have all of the depths in a vector. Find the variance of the vector
//remember, variance is the sum of the squared deviations from mean over 
  for(int i = 0; i < treeDepths.size(); i++){
	sumOfSquares+= ( (treeDepths.at(i) - ave) * (treeDepths.at(i) - ave));
	}
  return sumOfSquares/(double)treeDepths.size();
}


vector<unsigned int> checkForDuplicateBitstrings(){
  vector<unsigned int> retVec;
  for(unsigned int x = 0; x < biparttable.BipartTable.size()-1; x++){
	boost::dynamic_bitset<> b = biparttable.BipartTable[x].get_bitstring();
	//pad b with 0s
	while(b.size() < ::biparttable.lm.size()){
		b.push_back(0);
		}
	b.flip();
		for(unsigned int j = x+1; j < biparttable.BipartTable.size(); j++){
			boost::dynamic_bitset<> z = biparttable.BipartTable.at(j).get_bitstring();
				while(z.size() < ::biparttable.lm.size()){
				z.push_back(0);
				}		
			if(z==b){
				retVec.push_back(x);
				cout << "DUPLICATE FOUND AT: " << x << endl;
				}
			}
		}
return retVec;
}


  //NOTE- CURRENTLY ONLY WORKS ON HOMOGENEOUS DATA!
double hamming_distance_total(unsigned int tree1, unsigned int tree2){
 unsigned int total = 0; 
 pair<set<unsigned int>, set<unsigned int>> rfs = rfDistanceSet(tree1, tree2);
 set<unsigned int> rf, rf2;
 rf = rfs.first; rf2 = rfs.second;

 for(set<unsigned int>::iterator i = rf.begin(); i!=rf.end(); i++){ //unique biparts to tree 1 are rows
  	vector<int> h;
	boost::dynamic_bitset<> bi = biparttable.non_trunc_bitstring(*i);
	for(set<unsigned int>::iterator j = rf2.begin(); j!=rf2.end(); j++){ //unique biparts to tree 2 are columns
		int XORdistance = 0;		
		boost::dynamic_bitset<> bj = biparttable.non_trunc_bitstring(*j);
		for(int x = 0; x < bi.size(); x++){
			if(bi[x]^bj[x]) {		
				XORdistance++;
				}
			}
		total+=XORdistance;
		}
	}
  //note- we add to total rfDistance^2 number of times
  // since we want total to reflect the rf distance * average distance between any two given bipartitions,
  // we need to divide total by rf.size
 return (double)total / (double)rf.size();
}



unsigned int hamming_distance_greedy(unsigned int tree1, unsigned int tree2){
  //NOTE- CURRENTLY ONLY WORKS ON HOMOGENEOUS DATA!

//this is average hamming distance * rf distance
//in other words, an all-to-all bipartition hamming distance

 unsigned int total = 0; 

 pair<set<unsigned int>, set<unsigned int>> rfs = rfDistanceSet(tree1, tree2);
 set<unsigned int> rf, rf2;
 rf = rfs.first; rf2 = rfs.second;

 //we need to make sure the rf distance is non-0 or else the program will crash!
 if(rf.size()!=0){
	 //now, calculate the XOR matrix
	 vector< vector <int> > hamming;

	 for(set<unsigned int>::iterator i = rf.begin(); i!=rf.end(); i++){ //unique biparts to tree 1 are rows
	  	vector<int> h;
		boost::dynamic_bitset<> bi = biparttable.non_trunc_bitstring(*i);
		for(set<unsigned int>::iterator j = rf2.begin(); j!=rf2.end(); j++){ //unique biparts to tree 2 are columns
			int XORdistance = 0;		
			boost::dynamic_bitset<> bj = biparttable.non_trunc_bitstring(*j);
			for(int x = 0; x < bi.size(); x++){
				if(bi[x]^bj[x]) {		
					XORdistance++;
					}
				}
			h.push_back(XORdistance);
			}
		hamming.push_back(h); //add new row to the hamming matrix
		}

	  //now we have the matrix of hamming distances. 
	  //keep a sorted list of all the indeces we have to explore

	set<unsigned int> indices;
	  for(unsigned int i = 0; i < hamming.at(0).size(); i++){
		indices.insert(i);
		}
	  //go through the hamming vector for each bipartition
	  //find the min, then remove that index from indices
	  //also add it to the total


	  int min = 0;
	  for(unsigned int i = 0; i < hamming.size(); i++){
		min = hamming.at(i).at(* (indices.begin()) );
		int index = *( indices.begin() ); //holds where we find the min from this row	
		//cout << "min is: " << min << " index is: " << index << endl;	
		//printVectorCompact(hamming.at(i));
	
		for(set<unsigned int>::iterator j = indices.begin(); j!=indices.end(); j++){
			if(hamming.at(i).at(*j) < min) {
				min = hamming.at(i).at(*j);
				index = *j;
				}
			}
		//now we've gotten the min hamming distance from bipartition i. 
		//add its min to to the total and remove index from indices
		total+= min;
		set<unsigned int>::iterator it = indices.find(index);
		if(it==indices.end()){
			cout << "hamming_distance_greedy error: index not found\n";
			cout << "index is: " << index << endl;
			printSetCompactTwo(indices);
			printVectorCompact(hamming.at(i));
			cout << endl << endl;
			}
		else{
			indices.erase(indices.find(index));
			}
		}
	 return total;
  }
  else{
	return 0;
	} 
}


double hamming_distance_average(unsigned int tree1, unsigned int tree2){
  return hamming_distance_total(tree1, tree2) /  (double)rfDistance(tree1, tree2);
}


unsigned int hamming_distance_minimum(unsigned int tree1, unsigned int tree2){

 unsigned int total = 0; 
 pair<set<unsigned int>, set<unsigned int>> rfs = rfDistanceSet(tree1, tree2);
 set<unsigned int> rf, rf2;
 rf = rfs.first; rf2 = rfs.second;


  if(rf.size()!=0){
	 //now, calculate the XOR matrix
	 vector< set <int> > hamming;

	 for(set<unsigned int>::iterator i = rf.begin(); i!=rf.end(); i++){ //unique biparts to tree 1 are rows
	  	set<int> h;
		boost::dynamic_bitset<> bi = biparttable.non_trunc_bitstring(*i);
		for(set<unsigned int>::iterator j = rf2.begin(); j!=rf2.end(); j++){ //unique biparts to tree 2 are columns
			int XORdistance = 0;		
			boost::dynamic_bitset<> bj = biparttable.non_trunc_bitstring(*j);
			for(int x = 0; x < bi.size(); x++){
				if(bi[x]^bj[x]) {		
					XORdistance++;
					}
				}
			h.insert(XORdistance);
			}
		hamming.push_back(h); //add new row to the hamming matrix
		}

	 //now we have the a vector of sets for each bipartition. Take the min from each set (i.e. first element)
	  unsigned int acc = 0;
	  for(vector<set<int>>::iterator i = hamming.begin(); i!=hamming.end(); i++){
		if((*i).size()>0){
			acc+=*((*i).begin());
			}
		}
	  return acc;
	}
  else{
	return 0;
	}
}



unsigned int rfDistance(int tree1, int tree2){

  unsigned long total = 0;
  vector<unsigned int> rf, rf2;
  //rf is the one sided difference from t1 to t2. rf2 is from t2 to t1
  vector<unsigned int> t1 = biparttable.inverted_index.at(tree1);
  vector<unsigned int> t2 = biparttable.inverted_index.at(tree2);
  set_difference(t1.begin(), t1.end(), t2.begin(), t2.end(), inserter(rf, rf.begin()));
  set_difference(t2.begin(), t2.end(), t1.begin(), t1.end(), inserter(rf2, rf2.begin()));
  return (rf.size() + rf2.size()) / 2;

}

pair<set<unsigned int>, set<unsigned int>> rfDistanceSet(int tree1, int tree2){

  set<unsigned int> rf, rf2;
  //rf is the one sided difference from t1 to t2. rf2 is from t2 to t1
  vector<unsigned int> t1 = biparttable.inverted_index.at(tree1);
  vector<unsigned int> t2 = biparttable.inverted_index.at(tree2);
  set_difference(t1.begin(), t1.end(), t2.begin(), t2.end(), inserter(rf, rf.begin()));
  set_difference(t2.begin(), t2.end(), t1.begin(), t1.end(), inserter(rf2, rf2.begin()));
  return make_pair(rf, rf2);
}

void printRFset(int tree1, int tree2){
  pair< set<unsigned int>, set<unsigned int> > rfDistances = rfDistanceSet(tree1, tree2);
  cout << "Bipartitions that are in tree 1 but not tree 2: \n";
  printSet(rfDistances.first);
  cout << "Bipartitions that are in tree 2 but not tree 1: \n";
  printSet(rfDistances.second);

}


void testDistance(){
 
testHamming();
//testNewick();
//testDepthVariance();


}

void testHamming(){
 for(int i = 234; i < 236; i++){
	for(int j = 1000; j < 2000; j++){
		cout << "\nTesting hamming for " << i << " " << j << ":\n";
		cout << "	rf distance = " << rfDistance(i, j);
		cout << ", HDG = " << hamming_distance_greedy(i,j);
	//hamming_distance_greedy(i,j);
	//hamming_distance_total(i,j);
	//hamming_distance_total(i,j);
		cout << ", HDT = " << hamming_distance_total(i,j);
		cout << ", HDA = " << hamming_distance_average(i,j);
		cout << ", HDM = " << hamming_distance_total(i,j);

		}
	cout << "done with " << i << endl;
	}
cout << endl;


}

void testDepthVariance(){
  for(int i = 0; i < 1000; i++){
	double ave;
	cout << "Depth variance is: " << depth_variance(i) << endl;
	}

}

void testAverageDepth(){
double totalAverage = 0;
double expect = expected_average_depth(biparttable.lm.size());
for(int i = 0; i < biparttable.NumTrees; i++){
	double ave;
	cout << "expected average depth for " << i << " is: " << expect;
	cout << ", average depth is: " << (ave = average_depth(i));
	cout << ", difference is " << ave-expect << endl;
	totalAverage += ave;	
	}
cout << "TOTAL AVERAGE OVER " << biparttable.NumTrees << " TREES IS: " << totalAverage/(double)biparttable.NumTrees << endl; 
}

void testNewick(){
 for(int i = 0; i < 10; i++){
    string nwTree = to_newick(i);
    cout << "newick string for tree " << i << " is: " << nwTree << endl;
	}


}
