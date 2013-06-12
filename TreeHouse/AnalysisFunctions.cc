#include "AnalysisFunctions.h"


/*
 * Distinguishing biparts and extension to hetero. 
 * In a homog setting a DB is found in all the trees in one set and none
 * of the trees in the other. When in a homog setting trees with 
 * different taxa can't have the "same" bipartitions. We can either 
 * ignore the differing taxa, and use DB to make a comment about the 
 * structure of the trees with those taxa removed, or we can say trees 
 * with varying taxa have no bipartitions in common and calculate the 
 * value. I think we might want to compute both, but would need to 
 * better define DB to mean only one and create a new term for the 
 * other. 
*/

void generate_random_bt(){
	BipartitionTable rand_bt;
	rand_bt.create_random_bt();
	rand_bt.print_hashtable();
	return;
}

//Returns a vector representing the bipartitions which are in all input trees
std::vector<unsigned int> biparts_in_all_trees(set<unsigned int> inputtrees){
	std::vector<unsigned int> bipartsInAll; //vector to return
	for (unsigned int i = 0; i < ::biparttable.get_size(); i++){ //for each bipartition
		vector<unsigned int> temp = ::biparttable.searchtable[i]; //The list of trees for bipart at index i
		set<unsigned int> sinter;
		set_intersection (temp.begin(), temp.end(), inputtrees.begin(), inputtrees.end(), std::inserter(sinter, sinter.begin()));
		//Intersection of the vectors is the size of the input set: bipartition is in all trees
		if (sinter.size() == inputtrees.size() ){
			bipartsInAll.push_back(i);
		}
	}
	
	//for(std::vector<unsigned int>::iterator pos = bipartsInAll.begin(); pos != bipartsInAll.end(); ++pos) {
	//	cout << *pos << " ";
    //}
    //cout << endl;
	return bipartsInAll;
}

//Returns a vector of biparts which are in none of the input trees
std::vector<unsigned int> biparts_in_no_trees(set<unsigned int> inputtrees){
	std::vector<unsigned int> bipartsInNo; // vector to be returned
	for (unsigned int i = 0; i < ::biparttable.get_size(); i++){ //for each bipartition
		vector<unsigned int> temp = ::biparttable.searchtable[i]; //The list of trees for bipart at index i
		set<unsigned int> sinter;
		set_intersection (temp.begin(), temp.end(), inputtrees.begin(), inputtrees.end(), std::inserter(sinter, sinter.begin()));
	        //intersection of the vectors is empty: there are no trees with the bipartition
		if (sinter.size() == 0 ){
			bipartsInNo.push_back(i);
		}
	}
	
	//for(std::vector<unsigned int>::iterator pos = bipartsInNo.begin(); pos != bipartsInNo.end(); ++pos) {
	//	cout << *pos << " ";
    //}
    //cout << endl;
	return bipartsInNo;
}

//Outputs the vector of bipartitions which are in all of the trees in one set, and none of the trees in the other
std::vector<int> distinguishing_bipart(set<unsigned int> inputtrees1, set<unsigned int> inputtrees2){
	vector<int> result; //Vector to be returned
	
	vector<unsigned int> inAll1 = biparts_in_all_trees(inputtrees1);
	vector<unsigned int> inNo1 = biparts_in_no_trees(inputtrees1);
	vector<unsigned int> inAll2 = biparts_in_all_trees(inputtrees2);
	vector<unsigned int> inNo2 = biparts_in_no_trees(inputtrees2);

	vector<unsigned int> inAll;
	std::set_union(inAll1.begin(), inAll1.end(), inAll2.begin(), inAll2.end(), std::inserter(inAll, inAll.end()));

	vector<unsigned int> inNo;
	std::set_union(inNo1.begin(), inNo1.end(), inNo2.begin(), inNo2.end(), std::inserter(inNo, inNo.end()));
	//The intersection will contain all bipartions which are in one set of trees but not he other
	std::set_intersection(inAll.begin(), inAll.end(), inNo.begin(), inNo.end(), std::inserter(result, result.end()));
	//Prints each bipartition
	for (unsigned int i = 0; i < result.size(); i++){ //for each bipartition
		cout << "printing: ";
		printBipartition(i);
	}
	
	return result;
}

//Returns the vector of taxa in are present in all of the input trees
std::vector<string> taxa_in_all_trees(set<unsigned int> inputtrees){
	vector<string> allTaxa = get_all_taxa_vect();// Vector to return
	//For each input tree
	for(std::set<unsigned int>::iterator pos = inputtrees.begin(); pos != inputtrees.end(); ++pos) {
	vector<string> temp = get_taxa_in_tree(*pos);
	std::vector<string> s; //Stores intersection temporarily
	std::set_intersection(allTaxa.begin(), allTaxa.end(), temp.begin(), temp.end(),  std::inserter(s, s.end()));
	allTaxa.swap(s); //Places the value of s in allTaxa and allTaxa in s
	}
	return allTaxa;
}

//Returns a vector of the taxa which are present in none of the input trees
std::vector<string> taxa_in_no_trees(set<unsigned int> inputtrees){
	vector<string> allTaxa = get_all_taxa_vect(); //Vector to return
	//for each input tree	
	for(std::set<unsigned int>::iterator pos = inputtrees.begin(); pos != inputtrees.end(); ++pos) {
	vector<string> temp = get_taxa_in_tree(*pos); //Taxa in current tree
	std::vector<string> s; //Stores the intersection temporarily
	std::set_difference(allTaxa.begin(), allTaxa.end(), temp.begin(), temp.end(),  std::inserter(s, s.end()));
	allTaxa.swap(s); //Places the value of s in allTaxa and allTaxa in s
	}
	return allTaxa;
}

//The taxa that appear in all of one set and none of the other. 
std::vector<string> distinguishing_taxa(set<unsigned int> inputtrees1, set<unsigned int> inputtrees2){
	vector<string> result; //vector to be returned
	
	vector<string> inAll1 = taxa_in_all_trees(inputtrees1);
	vector<string> inNo1 = taxa_in_no_trees(inputtrees1);
	vector<string> inAll2 = taxa_in_all_trees(inputtrees2);
	vector<string> inNo2 = taxa_in_no_trees(inputtrees2);
	
	vector<string> inAll; //Taxa in all of the trees from both sets
	std::set_union(inAll1.begin(), inAll1.end(), inAll2.begin(), inAll2.end(), std::inserter(inAll, inAll.end()));

	vector<string> inNo; //Taxa in none of the trees from either set
	std::set_union(inNo1.begin(), inNo1.end(), inNo2.begin(), inNo2.end(), std::inserter(inNo, inNo.end()));
	//Intersection contains those that are in none of one set and all of the other
	std::set_intersection(inAll.begin(), inAll.end(), inNo.begin(), inNo.end(), std::inserter(result, result.end()));

	return result;
}




/*
//the taxa that only appear in one of the sets. 
std::vector<string> taxa_set_xor(set<unsigned int> inputtrees1, set<unsigned int> inputtrees2){
	set<string> taxaSet1;
	set<string> taxaSet2;
	
    for(std::set<unsigned int>::iterator pos = inputtrees1.begin(); pos != inputtrees1.end(); ++pos) {
		vector<string> temp = get_taxa_in_tree(*pos);
		taxaSet1.insert(temp.begin(), temp.end());
    }
    
	for(std::set<unsigned int>::iterator pos = inputtrees2.begin(); pos != inputtrees2.end(); ++pos) {
		vector<string> temp = get_taxa_in_tree(*pos);
		taxaSet2.insert(temp.begin(), temp.end());
    }
    
    std::vector<string> s; 
	std::set_symmetric_difference(taxaSet1.begin(), taxaSet1.end(), taxaSet2.begin(), taxaSet2.end(),  std::inserter(s, s.end()));
	return s;
}
*/

void bitpartitions_by_frequency(set<unsigned int> inputtrees, float threshold, vector< bool * > &consensus_bs, vector< float > &consensus_branchs, vector< unsigned int> &consensus_bs_sizes){
	for (unsigned int i = 0; i < ::biparttable.get_size(); i++){ //for each bipartition
		vector<unsigned int> temp = ::biparttable.searchtable[i];
		set<unsigned int> sinter;
		set_intersection (temp.begin(), temp.end(), inputtrees.begin(), inputtrees.end(), std::inserter( sinter, sinter.begin() ) );
		
		if (sinter.size() >= threshold){
			consensus_bs.push_back(::biparttable.bipartitions[i]);
			//consensus_branchs.push_back(::list_branches[i]);
			consensus_branchs.push_back(1.0);
			consensus_bs_sizes.push_back(::biparttable.length_of_bitstrings[i]);
		}
	}
}

unsigned int compute_threshold(unsigned int numberofTrees, float threshold){
	float temp = ceil(numberofTrees * (threshold / 100.00) );
	//cout << "temp ceil = " << temp << endl;
	unsigned int returnvalue = (unsigned int)temp;
	if (returnvalue > numberofTrees){
		returnvalue = numberofTrees;
	}
	
	if (returnvalue <= (float)numberofTrees/2.0){
		returnvalue += 1;
	}
	//cout << "Returnvalue = " << returnvalue << endl;
	return returnvalue;
}

string consen(set<unsigned int> inputtrees, float percent){
	int threshold; 
	vector< bool * > consensus_bs;
	vector< float > consensus_branchs;
	vector< unsigned int> consensus_bs_sizes;
	string consensus_tree;
	
	threshold = compute_threshold(inputtrees.size(), percent);

	bitpartitions_by_frequency(inputtrees, threshold, consensus_bs, consensus_branchs, consensus_bs_sizes);

	//cout << "consensus_bs = "<< consensus_bs.size() << endl;
	//cout << "consensus_branchs = "<< consensus_branchs.size() << endl;
	//cout << "consensus_bs_sizes = "<< consensus_bs_sizes.size() << endl;

	//need a function to find return the bitstrings that appear in a certain percent of trees. 
	
	consensus_tree = compute_tree(::lm, consensus_bs, consensus_branchs, 0, false, consensus_bs_sizes); // change false to true to enable branch lengths.

	for (unsigned int i=0; i<consensus_bs.size(); i++) {
	  for (unsigned int j=0; j<consensus_bs_sizes[i]; j++)
	    cout << consensus_bs[i][j];
	  cout << endl;
	}
	
    return consensus_tree;
}


//Returns consensus resolution rate for a set of trees and a given consensus strictness
float consensus_reso_rate(set<unsigned int> inputtrees, float percent){
	string temp = consen(inputtrees, percent);
	float result;
	//Runs the resolution rate function on the output of forming a consensus tree
	result = reso_rate(temp);

	return result;


}

//Returns the resolution rate of a given tree
float reso_rate(string inputtree){

	float result; // Float to be returned
	
	//vector< vector < int > > treei;

	//The bitstring representation of the input tree
	vector< vector <int> > treei = compute_bitstrings_h(inputtree); 	
		//Only valid because homogeneous treesets, same number of taxa in all trees
		//vector<string> temp = get_taxa_in_tree(0);
		float taxa = ::NUM_TAXA;
		//Computes the number of possible bipartitions
		float total_biparts = taxa - 3;
		//Computes the number of non-trivial bipartitions
		float num_tree_biparts = treei.size() - taxa;
			cout << "Number of Bipartitions = " << num_tree_biparts << endl;
		//Computes the resolution value
		result = num_tree_biparts / total_biparts;
			cout << "Result is : " << result << endl;	
	return result;



}


  
BipartitionTable least_conflict_bt(set<unsigned int> inputtrees) {
  BipartitionTable consen_bt = get_consen_bt(inputtrees, 100);
  vector<set <unsigned int> > all_incompat;
  for (unsigned int i = 0; i < ::biparttable.get_size(); i++) { //for each bipartition
      vector<unsigned int> incompat_indices_vect = find_incompat(::biparttable.bipartitions[i], ::biparttable.length_of_bitstrings[i], inputtrees);
      set<unsigned int> incompat_indices(incompat_indices_vect.begin(), incompat_indices_vect.end());
      all_incompat.push_back(incompat_indices);
  }
     
    for (unsigned int threshold = 1; threshold < (::biparttable.get_size()); threshold++) {
      for (unsigned int i = 0; i < ::biparttable.get_size(); i++) {
	if (all_incompat[i].size() == threshold) {
	  set<unsigned int> consen_biparts;
	  for (unsigned int j = 0; j < consen_bt.searchtable.size(); j++) {
	    consen_biparts.insert(consen_bt.global_indices[j]);
	  }
	  set<unsigned int> conflicts;
	  set_intersection(all_incompat[i].begin(), all_incompat[i].end(), consen_biparts.begin(), consen_biparts.end(), std::inserter(conflicts, conflicts.begin()));
	  set<unsigned int> sinter;
	  set_intersection(::biparttable.searchtable[i].begin(),::biparttable.searchtable[i].begin(), inputtrees.begin(), inputtrees.end(), std::inserter(sinter, sinter.begin()));
	  if (conflicts.size() == 0 && sinter.size() >= 0) {
	    int insert_ind = 0;
	    for (unsigned int k = 0; k < consen_bt.global_indices.size(); k++) {
	      if (i < consen_bt.global_indices[k]) {
		insert_ind = k;
		break;
	      }
	    }
	    consen_bt.bipartitions.insert(consen_bt.bipartitions.begin()+insert_ind, ::biparttable.bipartitions[i]);
	    //consensus_branchs.push_back(::list_branches[i]);
	    consen_bt.branches.insert(consen_bt.branches.begin()+insert_ind, 1.0);
	    consen_bt.length_of_bitstrings.insert(consen_bt.length_of_bitstrings.begin()+insert_ind, ::biparttable.length_of_bitstrings[i]);
	    vector<unsigned int> incompat_vect(all_incompat[i].begin(), all_incompat[i].end());
	    consen_bt.searchtable.insert(consen_bt.searchtable.begin()+insert_ind, incompat_vect);
	    consen_bt.global_indices.insert(consen_bt.global_indices.begin()+insert_ind, i);
	  }
	}
      }
    }
    return consen_bt;
}

string least_conflict(set<unsigned int> inputtrees) {
  BipartitionTable consen_bt = least_conflict_bt(inputtrees);
  string consensus_tree = compute_tree(::lm, consen_bt.bipartitions, consen_bt.branches, 0, false, consen_bt.length_of_bitstrings);
  consen_bt.print_hashtable();
  return consensus_tree;
}


BipartitionTable get_consen_bt(set<unsigned int> inputtrees, float percent) {
  BipartitionTable consen_bt;
  unsigned int threshold = compute_threshold(inputtrees.size(), percent);

  for (unsigned int i = 0; i < ::biparttable.get_size(); i++){ //for each bipartition
    vector<unsigned int> temp = ::biparttable.searchtable[i];
    set<unsigned int> sinter;
    set_intersection (temp.begin(), temp.end(), inputtrees.begin(), inputtrees.end(), std::inserter( sinter, sinter.begin() ) );
    if (sinter.size() >= threshold){
      consen_bt.bipartitions.push_back(::biparttable.bipartitions[i]);
      //consensus_branchs.push_back(::list_branches[i]);
      consen_bt.branches.push_back(1.0);
      consen_bt.length_of_bitstrings.push_back(::biparttable.length_of_bitstrings[i]);
      vector<unsigned int> incompat_indices;
      if (percent < 100)
	incompat_indices = find_incompat(::biparttable.bipartitions[i], ::biparttable.length_of_bitstrings[i], inputtrees);
      consen_bt.searchtable.push_back(incompat_indices);
      consen_bt.global_indices.push_back(i);
    }
  }
  return consen_bt;
}

BipartitionTable greedy_consen_bt(set<unsigned int> inputtrees, float percent) {
  BipartitionTable consen_bt = get_consen_bt(inputtrees, percent);

  set<unsigned int> all_incompat;
  for (unsigned int i = 0; i < consen_bt.searchtable.size(); i++) {
    set<unsigned int> incompat_indices(consen_bt.searchtable[i].begin(), consen_bt.searchtable[i].end());
    std::set_union(all_incompat.begin(), all_incompat.end(), incompat_indices.begin(), incompat_indices.end(), std::inserter(all_incompat, all_incompat.begin()));
  }

  set<unsigned int> all_indices;
  for (unsigned int i = 0; i < biparttable.bipartitions.size(); i++) {
    
    set<unsigned int> trees_with_bipart(biparttable.searchtable[i].begin(), biparttable.searchtable[i].end());
    set<unsigned int> isect;
    std::set_intersection(inputtrees.begin(), inputtrees.end(), trees_with_bipart.begin(), trees_with_bipart.end(), std::inserter(isect, isect.begin()));
    if (!isect.empty()) {
      all_indices.insert(i);
    }
  }

  for (unsigned int i = 0; i < consen_bt.global_indices.size(); i++)
    all_indices.erase(consen_bt.global_indices[i]);

  vector<unsigned int> all_compat;
  std::set_difference(all_indices.begin(), all_indices.end(), all_incompat.begin(), all_incompat.end(), back_inserter(all_compat));  

  for (unsigned int base_popularity = ::NUM_TREES; base_popularity >= 1; base_popularity--) { 
    unsigned int i = 0; 
    while (i < all_compat.size()) {
      if (::biparttable.searchtable[all_compat[i]].size() >= base_popularity) {
	int insert_ind;
	for (unsigned int j = 0; j < consen_bt.global_indices.size(); j++) {
	  if (all_compat[i] < consen_bt.global_indices[j]) {
	    insert_ind = j;
	    break;
	  }
	}
	consen_bt.bipartitions.insert(consen_bt.bipartitions.begin()+insert_ind, ::biparttable.bipartitions[all_compat[i]]);
	consen_bt.branches.insert(consen_bt.branches.begin()+insert_ind, 1.0);
	consen_bt.length_of_bitstrings.insert(consen_bt.length_of_bitstrings.begin()+insert_ind, ::biparttable.length_of_bitstrings[all_compat[i]]);
	vector<unsigned int> incompat_indices = find_incompat(::biparttable.bipartitions[all_compat[i]], ::biparttable.length_of_bitstrings[all_compat[i]], inputtrees);
	consen_bt.searchtable.insert(consen_bt.searchtable.begin()+insert_ind, incompat_indices);
	consen_bt.global_indices.insert(consen_bt.global_indices.begin()+insert_ind, all_compat[i]);
	
	set<unsigned int> compat_set(all_compat.begin(), all_compat.end());
	compat_set.erase(all_compat[i]);
	set<unsigned int> incompat_set(incompat_indices.begin(), incompat_indices.end());
	all_compat.clear();
	std::set_difference(compat_set.begin(), compat_set.end(), incompat_set.begin(), incompat_set.end(), back_inserter(all_compat));

	i = 0;
      } 
      else i++;
    }
  }
  return consen_bt;
}

string greedy_consen(set<unsigned int> inputtrees, float percent) {
  BipartitionTable consen_bt = greedy_consen_bt(inputtrees, percent);
  string consensus_tree = compute_tree(::lm, consen_bt.bipartitions, consen_bt.branches, 0, false, consen_bt.length_of_bitstrings);
  consen_bt.print_hashtable();
  return consensus_tree;
}

void test_trait_correlation(int t1ind, int t1val, int t2ind, int t2val, unsigned int tree, int iterations, string folder) {

  vector<int> dists_target;
  vector<int> dists_nontarget;
  vector<int> dep_taxa_init = get_taxa_with_trait(t2ind, t2val);
  vector<int> dep_taxa_next_init = get_taxa_without_trait(t2ind, t2val);
  vector<int> indep_taxa_init = get_taxa_with_trait(t1ind, t1val);
  
  vector<string> parsed_tree = parse_newick(to_newick(tree));
  set<int> taxa_in_tree;

  for (unsigned int i=0; i<parsed_tree.size(); i++) {
    if (parsed_tree[i] != "(" &&
	parsed_tree[i] != ")" &&
	parsed_tree[i] != "," &&
	parsed_tree[i] != ";" &&
	parsed_tree[i] != "") {
      string taxonstring = parsed_tree[i];
      taxa_in_tree.insert(::lm.position(taxonstring));
    }
  }
  
  for (unsigned int i=0; i<iterations; i++) {
    cout << "Iteration " << i << endl;
    int cycles = 0;
    set<int> tempset;
    set<int> dep_taxa (dep_taxa_init.begin(), dep_taxa_init.end());
    set<int> dep_taxa_next (dep_taxa_next_init.begin(), dep_taxa_next_init.end());
    set<int> indep_taxa (indep_taxa_init.begin(), indep_taxa_init.end());
    //reduce sets to only those taxa in specified tree
    std::set_intersection(dep_taxa.begin(), dep_taxa.end(), taxa_in_tree.begin(), taxa_in_tree.end(), std::inserter( tempset, tempset.begin()));
    dep_taxa = tempset;
    tempset.clear();
    std::set_intersection(dep_taxa_next.begin(), dep_taxa_next.end(), taxa_in_tree.begin(), taxa_in_tree.end(), std::inserter( tempset, tempset.begin()));
    dep_taxa_next = tempset;
    tempset.clear();
    std::set_intersection(indep_taxa.begin(), indep_taxa.end(), taxa_in_tree.begin(), taxa_in_tree.end(), std::inserter( tempset, tempset.begin()));
    indep_taxa = tempset;
    tempset.clear();

    while (dep_taxa.size() > 1 && indep_taxa.size() > 1) {
      cout << "     Cycle " << cycles << endl;
      cycles++;
      set<int>::const_iterator iit(dep_taxa.begin());
      advance(iit, rand() % dep_taxa.size());
      int dep_taxon = *iit;
      int closest_with_indep;
      int min_distance = std::numeric_limits<int>::max();
      //find closest taxon with independent trait
      if (indep_taxa.find(dep_taxon) != indep_taxa.end()) {
	closest_with_indep = dep_taxon;
	min_distance = 0;
	cout << "          " << ::lm.name(dep_taxon) << "/" << taxa_info[dep_taxon]->label << " has independent trait: ";
	for (unsigned int i=0; i<taxa_info[dep_taxon]->traits.size(); i++)
	  cout << ::taxa_info[dep_taxon]->traits[i];
	cout << endl;
      }
      else {
	for (set<int>::const_iterator dit(indep_taxa.begin()); dit != indep_taxa.end(); ++dit) {
	  if (dep_taxon == *dit)
	    continue;
	  int new_distance = distance_between_taxa(dep_taxon, *dit, tree);
	  if (new_distance < min_distance) {
	    closest_with_indep = *dit;
	    min_distance = new_distance; 
	  }
	}
      }
      
      vector<int> taxa_in_clade_vect;
      //exclude clade
      if (dep_taxon == closest_with_indep)
	taxa_in_clade_vect.push_back(dep_taxon);
      else {
	vector<int> taxavect;
	taxavect.push_back(dep_taxon);
	taxavect.push_back(closest_with_indep);
	taxa_in_clade_vect = get_taxa_in_clade(taxavect, tree);
      }
      set<int> taxa_in_clade (taxa_in_clade_vect.begin(), taxa_in_clade_vect.end());
      
      std::set_difference(dep_taxa.begin(), dep_taxa.end(), taxa_in_clade.begin(), taxa_in_clade.end(), std::inserter( tempset, tempset.begin()));
      dep_taxa = tempset;
      tempset.clear();
      std::set_difference(dep_taxa_next.begin(), dep_taxa_next.end(), taxa_in_clade.begin(), taxa_in_clade.end(), std::inserter( tempset, tempset.begin()));
      dep_taxa_next = tempset;
      tempset.clear();
      std::set_difference(indep_taxa.begin(), indep_taxa.end(), taxa_in_clade.begin(), taxa_in_clade.end(), std::inserter( tempset, tempset.begin()));
      indep_taxa = tempset;
      tempset.clear();
      
      //swap dep_taxa and dep_taxa_next: switch target and nontarget samples
      tempset = dep_taxa;
      dep_taxa = dep_taxa_next;
      dep_taxa_next = tempset;
      tempset.clear();
      
      //log distance-to-common-ancestor
      //make sure an equal number of target/nontarget distances will be recorded
      if (cycles % 2 == 0) {
	if (dep_taxon == closest_with_indep)
	  dists_nontarget.push_back(0);
	else
	  dists_nontarget.push_back(distance_to_common_ancestor(dep_taxon, closest_with_indep, tree) - 1);
	//cout << "          dist_nontarget=" << dist_nontarget << endl;
      }
      else if (dep_taxa.size() > 1 && indep_taxa.size() > 1) {
	if (dep_taxon == closest_with_indep)
	  dists_target.push_back(0);
	else
	  dists_target.push_back(distance_to_common_ancestor(dep_taxon, closest_with_indep, tree) - 1);
	//cout << "          dist_target=" << dist_target << endl;
      }
    } //end while

  } //end for


  // reporting (make a CSV file)

  system(("mkdir -p "+folder).c_str());
  
  ofstream output (folder+"/"+folder+"-"+to_string(iterations)+".csv");
  
  output << "target_AD,nontarget_AD" << endl;
  
  for(unsigned int i=0; i<dists_target.size()-1; i++) {
    output << dists_target[i] << ',' << dists_nontarget[i] << endl;
  }

  output << dists_target[dists_target.size()-1] << ',' << dists_nontarget[dists_target.size()-1];

  output.close();

  cout << "Total measurements recorded: " << dists_target.size() << endl;

}


  

    
    
    
      
      
	
      
    
    
