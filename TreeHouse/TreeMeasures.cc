#include "TreeMeasures.h"
using namespace std;




double calculate_C(unsigned int tree){
  return calculate_C(to_newick(tree));
}

double calculate_C(string nw){
  //C is defined as 2/( n(n-3)+2 ) * (summation from 1 to n-1)(ri - si)
  //where ri and si are the number of taxa to the right and left of each internal node 
  //and ri>si

  //theorem- in a newick string of a bifurcating tree, each opening parenthesis represents an internal node of the tree
  unsigned int total = 0;
  pair<int, int> x = numLeftRight(nw, total);
  cout << "calculate_C: total is: " << total << endl;
  size_t num_taxa = count(nw.begin(), nw.end(), ',') + 1; //count the number of commas, add one, to get NumTaxa
  double multiplier = 2.0 / ( ( (double)num_taxa * (double)(num_taxa-3) ) + 2.0);
  return multiplier * (double)total;  

}

pair<int, int> numLeftRight(string nw, unsigned int &total){
  //helper function for calculate_C, returns the number of taxa to the left and right of an internal node
  
  //remove spaces in string:
  nw.erase(std::remove_if(nw.begin(), nw.end(), ::isspace), nw.end());
//  cout << "CALLING NUMLEFTRIGHT ON " << nw << endl;

  //BASE CASE 1- we have a singleton taxa, as in "A". Return pair(1,0)
  //BASE CASE 2- we have a closing paren before we see another opening pair, as in the case of (A,B). 
  //	return pair(1,1)
  //RECURSIVE CASE- we encounter an opening paren, as is the case of ((A,B),(C,D)). 
	//since the tree is bifurcating, there must be at most TWO substrings, where the first is started
	//by the first opening paren, and the second is started after the closing of the first opening paren
	//for instance, the Newick String  ((A,B),(C,D),E,F) should never be encountered

  //now to detect each base:
  //base case 1- string starts with a single (
  //base case 2- string does not start with a paren
  //recursive case- string starts with ((

  if(nw.length() < 1){
	return make_pair(0,0);
	}
  //now we know the string has at least one character; check what the first one is
  if(nw.at(0)=='('){ //we're in either base case 1 or recursive case
	if(nw.size()<2){
		cout << "numLeftRight error: input string is '(', this shouldn't happen\n";
		return make_pair(0,0);
		}
	if(count(nw.begin()+1, nw.end(), '(') == 0){ //there are no more left parents, BASE CASE 2
		return make_pair(1,1);
		}
	else{ //RECURSIVE CASE
		//now we need to find the two branches of and call numLeftRight on both
		if(nw.at(1)!='('){
		//case one- we have a NwString like (A,(B,C)), i.e. second character in this string 
			//isn't a left paren. 
			int l = 1;
			int r = 0;
			//find the next opening paren
			size_t firstParen = nw.find("(", 1); //don't search first character
			if(firstParen==string::npos){
				cout << "ERROR- could not find opening paren in numLeftRight\n";
				}
			else{
				string right = nw.substr(firstParen, nw.length()-firstParen);
		//		cout << "recursive case 1- calling on " << right << endl;
				r = sumPair(numLeftRight(right, total));
				}
			total += abs (l - r);
		//	cout << abs (l - r) << " added to total for string " << nw << endl;
			return make_pair(l, r);
			}
		else{ //we know string starts with "(("
			//find the split between each branch of the tree
			int numParens = 1; int index = 2;
			while(numParens>0 && index < nw.length()){
				if(nw.at(index)=='('){
					numParens++;
					}
				else if(nw.at(index)==')'){
					numParens--;
					}
				if(numParens!=0){
					++index;
					}
				}
			//now we have the index of the closing paren after the first opening paren
			//there should be a comma after this index
			if(!nw.at(index+1)==','){
				cout << "recursive case error! expected comma to split up left and right!\n";
				}
			string leftString = nw.substr(1, index);
			string rightString = nw.substr(index+2, nw.length()-index-3);
			cout << "recursive call: performing numLeftRight on " << leftString << " and " << rightString << " from original " << nw << endl;
			int left = sumPair( numLeftRight(leftString, total) );
			int right = sumPair( numLeftRight( rightString, total) );
			total += abs (left - right);
	//		cout <<abs (left - right) << " added to total for string " << nw << endl;
			return make_pair(left,right);
			}	
		}
	}
  else{ //BASE CASE 1
	return make_pair(1,0);
	}
}

bool isBifurcating(string nw){
  return isBifurcating(nw, 0);
}

bool isBifurcating(string nw, int depth){
  //cout << "calling isBifurcating on " << nw << endl;
  bool retVal = true;

  //remove spaces in string:
  nw.erase(std::remove_if(nw.begin(), nw.end(), ::isspace), nw.end());
  if(nw.length()<2){
	cout << "isBifurcating error: newick string input size is " << nw.length() << endl;
	return false;
	}	

  //make sure the string starts with a paren
  if(nw.at(0)!='('){
	cout << "isBifurcating error! The newick string doesn't start with a paren! Starts with " << nw.at(0) << endl;
	return false;
	}
  //now we know that the 0th index is an opening paren. 

  int degree = 0;
  int index = 1;
 
  while(nw.length() > index){ //OUTER WHILE LOOP- looks for branches off of root node
	//all we care about is commas, which indicate we just hit a taxon, and parens, which indicate we found nesting
	if(nw.at(index)==','){
		//we just went over a taxon or a paren unit, add to degree
		degree++;
		if(degree>2){
		//cout << "The tree isn't bifurcating! Found at depth " << depth << endl;
		 // if(depth!=1){cout << " found at depth " << depth << endl;}		
		return false;
		}
		index++;
		}
	else if(nw.at(index)=='('){
		//we found a new level of nesting! Recursively call it and skip to when we find a closing paren
		int startIndex = index;	//store the index where we encounter the opening paren	
		depth++;
		int parenNesting = 1;
		//find where the matching closing paren is
		index++;
		while(parenNesting!=0){
			if(nw.at(index)=='('){
				parenNesting++;
				}
			else if(nw.at(index)==')'){
				parenNesting--;
				}
			if(parenNesting!=0){index++;}
			}
		int endIndex = index; //this is where the nested newick branch ends
		//we still have a little work to do- skip to after the next comma, because other this function won't work on ((A,B),(C,D)), for example
	/*	while(index < nw.length()){
			if(nw.at(index)!='('){			
				index++;
				}
			else{
				break;
				}
			}*/
		index++;
		retVal = retVal & isBifurcating(nw.substr(startIndex, endIndex-startIndex+1), depth);
		}
	else if(nw.at(index)==')'){
		//cout << "closing paren found at index " << index << " on string " << nw << endl;		
		degree++; //add to degree for instances like (A,B)
			  //the code should only hit a closing paren if it the Nw string didn't start with two opening parens
		index++;

		}
	else{
		index++;
		}
	}

 // cout << "on NwString " << nw << ", degree is " << degree << endl;
  return (degree==2);

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

double average_depth(set<unsigned int> trees){
  double total = 0;
  for(set<unsigned int>::iterator i = trees.begin(); i!=trees.end(); i++){
	total += average_depth(*i);
	}
  return total/(double)trees.size();

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

void testIsBifurcating(){
  unsigned int numBifurcating = 0;
  for(unsigned int i = 0; i < biparttable.NumTrees; i++){
	bool x = isBifurcating(to_newick(i));
	if(x){ numBifurcating++;
		cout << "Tree " << i << " is bifurcating" << endl;}
	if(!(i%100)) { cout << "DONE WITH " << i << endl;}
	}
  cout << "NUMBER OF TREES BIFURCATING: " << numBifurcating << endl;

/*  string nw = "((A,B),C)";
  cout << isBifurcating(nw) << endl;
  nw = "((A,B),(C,D))";
  cout << isBifurcating(nw) << endl;
  nw = "(A,B)";
  cout << isBifurcating(nw) << endl;
  nw = "(A,(B,C))";
  cout << isBifurcating(nw) << endl;
  //now the following tests should NOT be bifurcating
  nw = "((A,B),(C,D), (E,F))";
  cout << isBifurcating(nw) << endl;
  nw = "(A,B,C,D)";
  cout << isBifurcating(nw) << endl;
  nw = "(((A,B),(C,D)))";
  cout << isBifurcating(nw) << endl;
  nw = "(A)";
  cout << isBifurcating(nw) << endl;
*/
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

double average_distance_between_taxa(unsigned int taxon1, unsigned int taxon2, set<unsigned int> treeSet){
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
		//this means they ARENT on the same edge for that bitstring- add to total
	 	 	//cout << "tree size for bipartition " << i << " is " << biparttable.BipartitionTable[i].trees_size() << endl;
			vector<unsigned int> trees = biparttable.BipartTable[i].get_trees();
			for(int t = 0; t < trees.size(); t++){
				if(treeSet.find(trees.at(t))!=treeSet.end()){
					++total;
					}
				}

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

	  return total/treeSet.size();
	  
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
  for(set<unsigned int>::iterator i = treeSet.begin(); i!= treeSet.end(); i++){
	total+=distance_to_common_ancestor(taxon1, taxon2, *i);
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
