#include "distance.h"

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

double average_distance_between_taxa(unsigned int taxon1, unsigned int taxon2){
  //only homogeneous?

  if(!::HETERO){
  cout << "we have a homogeneous data set!" << endl;
  //takes the average distance between taxa for every tree.
  //first, we need to look at every bipartition and indicate whether the taxa are on a different side of the edge
  //each edge where the taxa differ adds 1 to the distance
  
  unsigned int numBiparts = biparttable.BipartitionTable.size();
  unsigned long int total = 0;

  for(int i = 0; i < numBiparts; i++){
	if(!biparttable.BipartitionTable[i].same_bitstring_value(taxon1, taxon2)){
 	 	//cout << "tree size for bipartition " << i << " is " << biparttable.BipartitionTable[i].trees_size() << endl;
		total+=biparttable.BipartitionTable[i].trees_size();
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

  return total/(double)::NUM_TREES;
  
 }
  else{
  	cout << "We have a heterogeneous data set! Return 0 for now!" << endl;
  	return 0;
	}

return 0;

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

double averageDepth(unsigned int tree){ //returns the average taxon depth of a tree.

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




void testDistance(){

cout << endl << endl;
//int taxon = 3; int tree = 2;
//cout << "Distance to root in taxon " << taxon << " in tree " << tree << " is: " << distance_to_root(taxon,tree) << endl;

//for(int tree = 0; tree < ::NUM_TREES; tree++){
//	cout << "Average depth of taxa in tree " << tree << " is: " << averageDepth(tree) << endl;
//	}

for(int i = 0; i < 500; i++)
{
cout << "average distance between taxa 0 and " << i << " is: " << average_distance_between_taxa(0, i) << endl;
 }

}
