// ConsensusFunctions.cc
#include "ConsensusFunctions.h"

//should really phase out
void bitpartitions_by_frequency(set<unsigned int> inputtrees, float threshold, vector< bool * > &consensus_bs, vector< float > &consensus_branchs, vector< unsigned int> &consensus_bs_sizes){
	for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //for each bipartition
		vector<unsigned int> temp = ::biparttable.get_trees(i);
		set<unsigned int> sinter;
		set_intersection (temp.begin(), temp.end(), inputtrees.begin(), inputtrees.end(), std::inserter( sinter, sinter.begin() ) );
		
		if (sinter.size() >= threshold){
			consensus_bs.push_back(::biparttable.get_boolarray(i));
			//consensus_branchs.push_back(::biparttable.get_branchlengths(i));
			consensus_bs_sizes.push_back(::biparttable.bitstring_size(i));
			//GRB Do we need a branch length? 
			consensus_branchs.push_back(1.0);
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
	vector<Bipartition> consensus_bipartitions;
	vector< bool * > consensus_bs;
	vector< float > consensus_branchs;
	vector< unsigned int> consensus_bs_sizes;
	string consensus_tree;

	threshold = compute_threshold(inputtrees.size(), percent);



	//GRB Need to come in and clean up the boolarrays this spawns
	
	bitpartitions_by_frequency(inputtrees, threshold, consensus_bs, consensus_branchs, consensus_bs_sizes);

	//cout << "consensus_bs = "<< consensus_bs.size() << endl;
	//cout << "consensus_branchs = "<< consensus_branchs.size() << endl;
	//cout << "consensus_bs_sizes = "<< consensus_bs_sizes.size() << endl;

	//need a function to find return the bitstrings that appear in a certain percent of trees. 

	
	consensus_tree = compute_tree(::biparttable.lm, consensus_bs, consensus_branchs, 0, false, consensus_bs_sizes); // change false to true to enable branch lengths.

	//consensus_tree = compute_tree(::lm, consensus_bs, consensus_branchs, 0, false, consensus_bs_sizes); // change false to true to enable branch lengths.
	
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
  BipartitionTable consen_bt = get_consen_bt(inputtrees, 100); // start with a strict consensus bipartition table. 
  // these are the bipartition that don't conflict with any other bipartition (threshold 0)
  
  vector<set <unsigned int> > all_incompat;
  for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++) { //for each bipartition find the bipartitions that are incompatible with it. 
      vector<unsigned int> incompat_indices_vect = find_incompat(::biparttable.get_bitstring(i), inputtrees);
      set<unsigned int> incompat_indices(incompat_indices_vect.begin(), incompat_indices_vect.end());
      all_incompat.push_back(incompat_indices);
  }

  for (unsigned int threshold = 1; threshold < (::biparttable.biparttable_size()); threshold++) { //iterate the thresholds, threshold 1 are biparts that conflict with only one other bipartition
    for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++) {
	  if (all_incompat[i].size() == threshold) {
	    set<unsigned int> consen_biparts;
	    for (unsigned int j = 0; j < consen_bt.biparttable_size(); j++) {
	      consen_biparts.insert(consen_bt.global_indices[j]);
	    }
	    
	    set<unsigned int> conflicts;
	    set_intersection(all_incompat[i].begin(), all_incompat[i].end(), consen_biparts.begin(), consen_biparts.end(), std::inserter(conflicts, conflicts.begin()));
	    set<unsigned int> sinter;
	    //BUG ? begin was used twice in a row
	    set_intersection(::biparttable.trees_begin(i),::biparttable.trees_end(i), inputtrees.begin(), inputtrees.end(), std::inserter(sinter, sinter.begin()));
	    if (conflicts.size() == 0 && sinter.size() >= 0) {
	      int insert_ind = 0;
	      for (unsigned int k = 0; k < consen_bt.global_indices.size(); k++) {
	        if (i < consen_bt.global_indices[k]) {
		      insert_ind = k;
		      break;
	        }
	      }
	      //consen_bt.BipartitionTable.insert(consen_bt.BipartitionTable.begin()+insert_ind, ::biparttable.get_bipartition(i));
	      //consensus_branchs.push_back(::list_branches[i]);
	      //consen_bt.branches.insert(consen_bt.branches.begin()+insert_ind, 1.0);
	      //consen_bt.length_of_bitstrings.insert(consen_bt.length_of_bitstrings.begin()+insert_ind, ::biparttable.bitstring_size(i));
	      vector<unsigned int> incompat_vect(all_incompat[i].begin(), all_incompat[i].end());
	      Bipartition tempbipart(::biparttable.get_bitstring(i), incompat_vect, 1.0);
	      consen_bt.BipartitionTable.insert(consen_bt.BipartitionTable.begin()+insert_ind, tempbipart);
	      //consen_bt.searchtable.insert(consen_bt.searchtable.begin()+insert_ind, incompat_vect);
	      consen_bt.global_indices.insert(consen_bt.global_indices.begin()+insert_ind, i);
	    }
	  }
    }
  }
  return consen_bt;
}

string least_conflict(set<unsigned int> inputtrees) {
  BipartitionTable consen_bt = least_conflict_bt(inputtrees);
  vector<bool *> biparts = consen_bt.get_compute_tree_bipartitions();
  string consensus_tree = compute_tree(::biparttable.lm, biparts, consen_bt.get_compute_tree_bipartitions_branchlens(), 0, false, consen_bt.get_compute_tree_bipartitions_bitlens());
  consen_bt.print_biparttable();
  //cleaning up the mess made by using bool *'s
  while(! biparts.empty()){
		delete[] biparts.back();
		biparts.pop_back();
	}
  
  return consensus_tree;
}


BipartitionTable get_consen_bt(set<unsigned int> inputtrees, float percent) {
  BipartitionTable consen_bt;
  unsigned int threshold = compute_threshold(inputtrees.size(), percent);

  for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //for each bipartition
    vector<unsigned int> temp = ::biparttable.get_trees(i);
    set<unsigned int> sinter;
    set_intersection (temp.begin(), temp.end(), inputtrees.begin(), inputtrees.end(), std::inserter( sinter, sinter.begin() ) );
    if (sinter.size() >= threshold){
		
      //consen_bt.bipartitions.push_back(::biparttable.bipartitions[i]);
      //consensus_branchs.push_back(::list_branches[i]);
      //consen_bt.branches.push_back(1.0);
      //consen_bt.length_of_bitstrings.push_back(::biparttable.length_of_bitstrings[i]);
      vector<unsigned int> incompat_indices;
      if (percent < 100){
	    incompat_indices = find_incompat(::biparttable.get_bitstring(i), inputtrees);
      }
      
      //consen_bt.searchtable.push_back(incompat_indices);

      Bipartition tempbipart(::biparttable.get_bitstring(i), incompat_indices, 1.0);
      consen_bt.BipartitionTable.push_back(tempbipart);

      consen_bt.global_indices.push_back(i);
      
    }
  }
  return consen_bt;
}

BipartitionTable greedy_consen_bt(set<unsigned int> inputtrees, float percent) {
  BipartitionTable consen_bt = get_consen_bt(inputtrees, percent);

  set<unsigned int> all_incompat;
  for (unsigned int i = 0; i < consen_bt.biparttable_size(); i++) {
    set<unsigned int> incompat_indices(consen_bt.trees_begin(i), consen_bt.trees_end(i));
    std::set_union(all_incompat.begin(), all_incompat.end(), incompat_indices.begin(), incompat_indices.end(), std::inserter(all_incompat, all_incompat.begin()));
  }

  set<unsigned int> all_indices;
  for (unsigned int i = 0; i < biparttable.biparttable_size(); i++) {
    set<unsigned int> trees_with_bipart(biparttable.trees_begin(i), biparttable.trees_end(i));
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
      if (::biparttable.trees_size(all_compat[i]) >= base_popularity) {
	int insert_ind;
	for (unsigned int j = 0; j < consen_bt.global_indices.size(); j++) {
	  if (all_compat[i] < consen_bt.global_indices[j]) {
	    insert_ind = j;
	    break;
	  }
	}
	
	//consen_bt.bipartitions.insert(consen_bt.bipartitions.begin()+insert_ind, ::biparttable.bipartitions[all_compat[i]]);
	//consen_bt.branches.insert(consen_bt.branches.begin()+insert_ind, 1.0);
	//consen_bt.length_of_bitstrings.insert(consen_bt.length_of_bitstrings.begin()+insert_ind, ::biparttable.length_of_bitstrings[all_compat[i]]);
	vector<unsigned int> incompat_indices = find_incompat(::biparttable.get_bitstring(all_compat[i]), inputtrees);
	//consen_bt.searchtable.insert(consen_bt.searchtable.begin()+insert_ind, incompat_indices);
	Bipartition tempbipart(::biparttable.get_bitstring(all_compat[i]), incompat_indices, 1.0);
    consen_bt.BipartitionTable.insert(consen_bt.BipartitionTable.begin()+insert_ind, tempbipart);
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
  vector<bool *> biparts = consen_bt.get_compute_tree_bipartitions();
  
  string consensus_tree = compute_tree(::biparttable.lm, biparts, consen_bt.get_compute_tree_bipartitions_branchlens(), 0, false, consen_bt.get_compute_tree_bipartitions_bitlens());
  consen_bt.print_biparttable();

  //cleaning up the mess made by using bool *'s
  while(! biparts.empty()){
	delete[] biparts.back();
	biparts.pop_back();
  }
  
  return consensus_tree;
 
}

