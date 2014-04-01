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

 //GRB NEW
 
  void print_branchlength_stats(){

	typedef std::map< boost::dynamic_bitset<>, TreeSet >::iterator clade_it_type;
 
	for(clade_it_type iter = biparttable.MapBenchMarks[1]; iter != biparttable.CladeMap.end(); iter++) {
		cout << biparttable.get_clade_string(iter) << " " << biparttable.get_mean_branchlength(iter) << endl;
	}
 }

 
 
 void print_summary_stats(){

	typedef std::map< boost::dynamic_bitset<>, TreeSet >::iterator clade_it_type;
 
	vector< unsigned int > shared_clades;
	
	vector< unsigned int > whos_responsible;
	vector< unsigned int > taxa_overlap;
	unsigned int num_overlapping_taxa;
	

	unsigned int trees_resposible1;
	unsigned int trees_resposible10;



	for(size_t i = 0; i < biparttable.NumTrees; i++){
		shared_clades.push_back(0);
	}

	for(size_t i = 0; i < biparttable.NumTrees; i++){
		whos_responsible.push_back(0);
	}

	for(size_t i = 0; i < biparttable.lm.size(); i++){
		taxa_overlap.push_back(0);
	}	


	for(clade_it_type iter = biparttable.MapBenchMarks[1]; iter != biparttable.MapBenchMarks[2]; iter++) {
		
		set<unsigned int> temp = biparttable.get_trees(iter);

		if(temp.size() > 1){
			set<unsigned int>::iterator myIterator;
			for(myIterator = temp.begin(); myIterator != temp.end(); myIterator++){
				taxa_overlap[iter->first.find_first()]+=1;
			}
		}
	}
	
	
	

	for(clade_it_type iter = biparttable.MapBenchMarks[2]; iter != biparttable.CladeMap.end(); iter++) {
	
		set<unsigned int> temp = biparttable.get_trees(iter);
	
		shared_clades[temp.size()] += 1;
		
		if(temp.size() > 1){
			set<unsigned int>::iterator myIterator;
			for(myIterator = temp.begin(); myIterator != temp.end(); myIterator++){
				whos_responsible[*myIterator]+=1;
			}
		}
	}
 
	for(size_t i = 0; i < biparttable.NumTrees;i++){
		if(shared_clades[i] > 0){
			cout << shared_clades[i] << " clades shared by " << i << " trees" << endl;  
		}
	}

	for(size_t i = 0; i < biparttable.lm.size();i++){
		if(taxa_overlap[i] > 0){
			num_overlapping_taxa+=1;
			//cout << " taxa " << i << " shared by " << taxa_overlap[i] << " trees" << endl;  
		}
	}
	
	cout << "The taxa that appears in the most trees appears in " << *std::max_element(taxa_overlap.begin(),taxa_overlap.end())<< " trees" << endl;
	cout << "Of the " << biparttable.lm.size() << " taxa in the set, " <<  num_overlapping_taxa << " appear in more than one tree." << endl;
	
 	//for(size_t i = 0; i < biparttable.NumTrees;i++){
	//	if(whos_responsible[i] > 0){
	//		trees_resposible1 +=1;
	//		cout << "Tree " << i << " is responsible for " << whos_responsible[i] << " shared clades" << endl;  
	//	}

	//	if(whos_responsible[i] > 9){
	//		trees_resposible10 +=1;
	//	}			
	//}
	
	//cout << "of the " << biparttable.NumTrees << " trees in the set only " <<  trees_resposible1 << " of them are share at least one clade with another tree " << endl;
	//cout << "of the " << biparttable.NumTrees << " trees in the set only " <<  trees_resposible10 << " of them are share at least 10 clades with another tree " << endl;
 
 
	//for(size_t i = 0; i < biparttable.NumTrees;i++){
	//	cout << "Tree " << i << " has " << get_outgroup_ids(i).size() << " Outgroup Taxa" << endl;
	//}
	
	print_branchlength_stats();
 }
 

 
 
 
 
//GRB NEW 
string psupport(vector < set < unsigned int > > treesets){
	typedef std::map< boost::dynamic_bitset<>, TreeSet >::reverse_iterator clade_rit_type;	
	vector<float> bipartition_percent;
	vector<float> bipartition_psupport;

	vector<boost::dynamic_bitset<> > maj_clades;
	vector<string> clade_psupport;

	
	
	set<unsigned int> trees_we_care_about;
	
	for(unsigned int i = 0; i < treesets.size(); i++){
		trees_we_care_about.insert(treesets[i].begin(), treesets[i].end());
	}

	cout << "trees_we_care_about.size = " << trees_we_care_about.size() << endl;

	int counter = 0;
	
	for(clade_rit_type riter = biparttable.CladeMap.rbegin(); riter != biparttable.CladeMap.rend(); riter++) {
		set<unsigned int> trees_at_bipart = biparttable.get_trees(riter);
		set<unsigned int> res_set;
		std::set_intersection(trees_at_bipart.begin(), trees_at_bipart.end(), trees_we_care_about.begin(), trees_we_care_about.end(),  std::inserter(res_set, res_set.end()));
		
		cout << counter << " ";
		
		float majsup = res_set.size()/float(trees_we_care_about.size());
		bipartition_percent.push_back(majsup);
		
		cout << majsup << " ";		
		
		unsigned int supportingClusts = 0;
		
		for(unsigned int j = 0; j < treesets.size(); j++){
			res_set.clear();
			std::set_intersection(trees_at_bipart.begin(), trees_at_bipart.end(), treesets[j].begin(), treesets[j].end(),  std::inserter(res_set, res_set.end()));
			
			cout << (res_set.size() / (float)treesets[j].size()) << " ";
			
			if(res_set.size() > (treesets[j].size()/2.0)){
				supportingClusts++;
			}
		}
		cout << "supportingClusts = " << supportingClusts << endl;

		bipartition_psupport.push_back(supportingClusts/float(treesets.size()));		
		
		if(majsup > .5){
			maj_clades.push_back(biparttable.get_bitstring(riter));
			
			std::ostringstream ss;
			ss << "b" << counter << "_" << (100 * supportingClusts/float(treesets.size())) << "%";
			//std::string s(ss.str());
			
			
			clade_psupport.push_back(ss.str());
		}
		counter++;
	}
	
	//for(unsigned int i = 0; i < biparttable.biparttable_size(); i++){
	//	cout << i << " " << bipartition_percent[i] << " " << bipartition_psupport[i] << endl;
	//}
	
	
	
	cout << maj_clades.size() << endl;
	
	//string newick = compute_tree_sup(::biparttable.lm, maj_clades, clade_psupport);
	
	string newick = compute_tree_labels(::biparttable.lm, maj_clades, clade_psupport);
	
	
	cout << newick << endl;
	//string newick;
	return newick;
	
	
}


float entropy(vector <set < unsigned int >> treesets){
	float ent_result = 0.0;
	
	unsigned int set_size = 0;
	
	for (unsigned int i = 0; i < treesets.size(); i++ ){
		set_size += treesets[i].size();
	}
	
	for (unsigned int i = 0; i < treesets.size(); i++ ){
		ent_result += (-1) * (treesets[i].size() / float(set_size)) * (log2(treesets[i].size() / float(set_size)) );
	}
	
	return ent_result;
}


float info_gain(vector <set < unsigned int >> treesets, unsigned int bpid){
	
	float infogain_of_s = entropy(treesets);

	vector<unsigned int> trees_at_bipart = biparttable.get_trees(bpid);

	unsigned int set_size = 0;
	for (unsigned int i = 0; i < treesets.size(); i++ ){
		set_size += treesets[i].size();
	}
	
	for (unsigned int i = 0; i < treesets.size(); i++ ){
	
		set<unsigned int> have_set;
		std::set_intersection(trees_at_bipart.begin(), trees_at_bipart.end(), treesets[i].begin(), treesets[i].end(),  std::inserter(have_set, have_set.end()));
	
		set<unsigned int> havenot_set;
		std::set_difference(treesets[i].begin(), treesets[i].end(), have_set.begin(), have_set.end(),  std::inserter(havenot_set, havenot_set.end()));

		vector<set <unsigned int> > haves_and_nots;
		
		haves_and_nots.push_back(have_set);
		haves_and_nots.push_back(havenot_set);
	
		infogain_of_s -= ( treesets[i].size() / float(set_size)) * entropy(haves_and_nots);

	}
	return infogain_of_s;
}


vector<vector<set<unsigned int>>>  split(vector <set < unsigned int >> treesets,  unsigned int bpid){
	vector<unsigned int> trees_at_bipart = biparttable.get_trees(bpid);

	vector<vector<set<unsigned int>>> splits;

	for (unsigned int i = 0; i < treesets.size(); i++ ){
	
		set<unsigned int> have_set;
		std::set_intersection(trees_at_bipart.begin(), trees_at_bipart.end(), treesets[i].begin(), treesets[i].end(),  std::inserter(have_set, have_set.end()));
	
		set<unsigned int> havenot_set;
		std::set_difference(treesets[i].begin(), treesets[i].end(), have_set.begin(), have_set.end(),  std::inserter(havenot_set, havenot_set.end()));

		vector<set <unsigned int> > haves_and_nots;
		
		haves_and_nots.push_back(have_set);
		haves_and_nots.push_back(havenot_set);
	
		splits.push_back(haves_and_nots);
	}
	
	return splits;
}

/*
void ID3wrapper(vector <set < unsigned int >> treesets){
	
	std::map<unsigned int, unsigned int> MyMap;	
	
	for(unsigned int i = 0; i < treesets.size(); i++){				
		set<unsigned int>::iterator myIterator;
		for(myIterator = treesets[i].begin(); myIterator != treesets[i].end(); myIterator++){
			MyMap[*myIterator]=i;
		}
	}
	
	vector <unsigned int> bipartitions;
	for(unsigned int i = 0; i < biparttable.biparttable_size(); i++){				
		bipartitions.push_back(i);
	}
}
* */

/*
void ID3 (std::map<unsigned int, unsigned int> MyMap, vector <unsigned int> bipartitions){

	if(bipartitions.size() == 0){
		cout << "We ran out of bipartitions" << endl;
	}
	
	for (auto it=mymap.begin(); it!=mymap.end(); ++it){
		
	}
	
	for(unsigned int i = 0; i < treesets.size(); i++){		
		if(MyMap.size() ==  0 )
		cout << "got to the end of a path" << endl;
	}		
	
	choose_attribute()
}
*/



/*
void dTree(vector <set < unsigned int >> treesets){

	vector<int> depth;
	int calls_to_split = 0;
	
	int best_depth = 0;
	
	vector <vector <set < unsigned int >>> to_split;

	depth.push_back(0);
	to_split.push_back(treesets);
	
	
	
	while (! to_split.empty()){
		vector<set<unsigned int>> working_trees = to_split.pop_back();

		float best_score = 0.0;
		int best_bip = -1;
		float score = 0.0;
		
		for (unsigned int i = 0; i < biparttable.biparttable_size(); i++){
			score = info_gain(treesets, i);
			if (score > best_score){
				best_score = score;
				best_bip = i;
			}
		}
		vector<vector<set<unsigned int>>> splits = split(working_trees, best_bip);
		
		for (unsigned int i = 0; i < splits.size(); i++){
			for(unsigned int j = 0; j < splits[i].size(); j++){
				if(check(splits[i])){
					
				}
				else{
					
				}
				
			}
		}
	}

	cout << "best score = " << best_score << ", best_bip = " << best_bip << endl;
}
* */

/*
float info_gain(set < unsigned int > treeset, unsigned int bpid){
	float info_gain = 0.0;
	float ent = entropy(treeset, bpid);
	vector<unsigned int> trees_at_bipart = biparttable.get_trees(bpid);

	set<unsigned int> have_set;
	std::set_intersection(trees_at_bipart.begin(), trees_at_bipart.end(), treeset.begin(), treeset.end(),  std::inserter(have_set, have_set.end()));
	float s_have = (have_set.size()/float(treeset.size())) *  entropy(have_set, bpid);

	
	set<unsigned int> have_not_set;
	std::set_intersection(trees_at_bipart.begin(), trees_at_bipart.end(), treeset.begin(), treeset.end(),  std::inserter(have_not_set, have_not_set.end()));
	float s_have_not = (have_not_set.size()/float(treeset.size())) *  entropy(have_not_set, bpid);

	info_gain = ent - ( s_have + s_have_not);
	
	return info_gain;
}


float bipart_entropy(set < unsigned int > treeset, unsigned int bpid){
	float ent_result = 0.0;
	
	vector<unsigned int> trees_at_bipart = biparttable.get_trees(bpid);
	set<unsigned int> res_set;
	std::set_intersection(trees_at_bipart.begin(), trees_at_bipart.end(), treeset.begin(), treeset.end(),  std::inserter(res_set, res_set.end()));

	ent_result = (-1)*(res_set.size()/float(treeset.size())) * log2(res_set.size()/float(treeset.size())) + (-1)*((treeset.size()-res_set.size())/float(treeset.size())) * log2((treeset.size()-res_set.size())/float(treeset.size()));
		
	return ent_result;
}
*/




void generate_random_bt(){
	BipartitionTable rand_bt;
	rand_bt.create_random_bt();
	rand_bt.print_biparttable();
	return;
}

//GRB: Should be moved or deleted... This is not really analysis function
// also there's an inverted index for this exact thing.
//Returns a vector representing the bipartitions which are in the input tree
std::vector<unsigned int> biparts_in_tree(unsigned int inputtree){
	std::vector<unsigned int> biparts; //vector to return
	for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //for each bipartition
		vector<unsigned int> temp = ::biparttable.get_trees(i); //list of trees for bipart
		if (std::find(temp.begin(), temp.end(), inputtree) != temp.end()){//If tree is in biparts list of trees
			biparts.push_back(i);
		}
	}
	return biparts;
}


//Returns a vector representing the bipartitions which are in all input trees
std::vector<unsigned int> biparts_in_all_trees(set<unsigned int> inputtrees){
	std::vector<unsigned int> bipartsInAll;
	for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //for each bipartition
		vector<unsigned int> temp = ::biparttable.get_trees(i);

		std::vector<unsigned int> bipartsInAll; //vector to return
		for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //for each bipartition
			vector<unsigned int> temp = ::biparttable.get_trees(i); //The list of trees for bipart at index i
			set<unsigned int> sinter;
			set_intersection (temp.begin(), temp.end(), inputtrees.begin(), inputtrees.end(), std::inserter(sinter, sinter.begin()));
			//Intersection of the vectors is the size of the input set: bipartition is in all trees
			if (sinter.size() == inputtrees.size() ){
				bipartsInAll.push_back(i);
			}
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
	std::vector<unsigned int> bipartsInNo;
	for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //for each bipartition
		vector<unsigned int> temp = ::biparttable.get_trees(i);
		std::vector<unsigned int> bipartsInNo; // vector to be returned
		for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //for each bipartition
			vector<unsigned int> temp = ::biparttable.get_trees(i); //The list of trees for bipart at index i
			set<unsigned int> sinter;
			set_intersection (temp.begin(), temp.end(), inputtrees.begin(), inputtrees.end(), std::inserter(sinter, sinter.begin()));
			//intersection of the vectors is empty: there are no trees with the bipartition
			if (sinter.size() == 0 ){
				bipartsInNo.push_back(i);
			}
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






//Returns the vector of taxa which are present in all of the input trees
std::vector<string> taxa_in_all_trees(set<unsigned int> inputtrees){
	vector<string> allTaxa = ::biparttable.lm.get_all_taxa_vect();// Vector to return
	//For each input tree
	for(std::set<unsigned int>::iterator pos = inputtrees.begin(); pos != inputtrees.end(); ++pos) {
		vector<string> temp = ::biparttable.get_taxa_in_tree(*pos);
		std::vector<string> s; //Stores intersection temporarily
		std::set_intersection(allTaxa.begin(), allTaxa.end(), temp.begin(), temp.end(),  std::inserter(s, s.end()));
		allTaxa.swap(s); //Places the value of s in allTaxa and allTaxa in s
	}
	return allTaxa;
}

//Returns a vector of the taxa which are present in none of the input trees
std::vector<string> taxa_in_no_trees(set<unsigned int> inputtrees){
	vector<string> allTaxa = ::biparttable.lm.get_all_taxa_vect(); //Vector to return
	//for each input tree	
	for(std::set<unsigned int>::iterator pos = inputtrees.begin(); pos != inputtrees.end(); ++pos) {
		vector<string> temp = ::biparttable.get_taxa_in_tree(*pos); //Taxa in current tree
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

/*
   void bitpartitions_by_frequency(set<unsigned int> inputtrees, float threshold, vector< bool * > &consensus_bs, vector< float > &consensus_branchs, vector< unsigned int> &consensus_bs_sizes){
   for (unsigned int i = 0; i < ::biparttable.biparttable_size(); i++){ //for each bipartition
   vector<unsigned int> temp = ::biparttable.get_trees(i);
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
*/

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
			taxa_in_tree.insert(::biparttable.lm.position(taxonstring));
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
			cout << " Cycle " << cycles << endl;
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
				cout << " " << ::biparttable.lm.name(dep_taxon) << "/" << taxa_info[dep_taxon]->label << " has independent trait: ";
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
				//cout << " dist_nontarget=" << dist_nontarget << endl;
			}
			else if (dep_taxa.size() > 1 && indep_taxa.size() > 1) {
				if (dep_taxon == closest_with_indep)
					dists_target.push_back(0);
				else
					dists_target.push_back(distance_to_common_ancestor(dep_taxon, closest_with_indep, tree) - 1);
				//cout << " dist_target=" << dist_target << endl;
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

set <unsigned int> get_outgroup_ids(unsigned int tree){
		
	typedef std::map< boost::dynamic_bitset<>, TreeSet >::iterator clade_it_type;
	
	//get indexes for all taxa in tree
	set<unsigned int> taxa_ids = biparttable.get_taxa_in_tree_ids(tree);
	unsigned int numtaxa = taxa_ids.size();
	unsigned int cladesused = 0;
	
	for(clade_it_type iter = biparttable.MapBenchMarks[2]; iter != biparttable.CladeMap.end(); iter++) {
		if (biparttable.contains_tree(iter, tree)){
			cladesused++;
			//boost:dynamic_bitset<> bitset = iter->first.size();
			int i = iter->first.find_first();
			while(i != string::npos){
				taxa_ids.erase(i);
				i = iter->first.find_next(i);
			}
		}
	}
	cout << "Tree has " << cladesused << " clades and " << numtaxa << " taxa" <<  endl;
	return taxa_ids;
}


