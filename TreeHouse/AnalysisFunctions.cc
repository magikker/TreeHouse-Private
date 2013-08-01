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

//Returns the number of unique bipartitions (essentially counts each bipartition once, 
//ignoring duplicates as it goes)
unsigned int num_unique_biparts(vector < vector < unsigned int > > biparts){
	unsigned int retvalue;
	set<unsigned int > tempset;

	//Add the bipartitions to a set which does not store duplicates
	for(unsigned int i = 0; i < biparts.size(); i++){//for each trees bipartitions
		for(unsigned int j = 0; j < biparts[i].size(); j++){//for ech bipartition
			tempset.insert(biparts[i][j]);
		}
	}
	//Size of the set is the number of unique bipartitions
	retvalue = tempset.size();
	return retvalue;
}


//Computes various distance measures based on bipartitions, used to compute the distance matrix of given trees through two wrapper functions
vector < vector <unsigned int> > bipart_distances(vector < vector <unsigned int> > biparts, unsigned int switch_value){
	//Return Value
	vector< vector <unsigned int> > distances;

	distances.resize(biparts.size(), vector< unsigned int >(biparts.size(), 0));

	//To set the switch since strings are intuitive to us but not switch statements


	for(unsigned int i = 0; i < biparts.size() - 1; i++){//for each tree's bipartitions
		for (unsigned int j = 1; j < biparts.size(); j++){//for all others
			//a, b, & c values come from Suzanne Matthews Dissertation
			//and the method of computing distances from bipartitions found there
			int a; //Bipartitions in both trees
			vector<unsigned int> temp;
			int b; //# Bipartions in first tree but not second
			int c; //# Bipartitions in second tree but not first
			int m; //# Unique Bipartitions
			int dist;
			//Intersection contains all shared bipartitions
			std::set_intersection(biparts[i].begin(),biparts[i].end(),
					biparts[j].begin(),biparts[j].end(),
					std::inserter(temp, temp.end()));
			//Stores the differences for each
			a = temp.size();
			b = biparts[i].size() - temp.size();
			c = biparts[j].size() - temp.size();
			m = num_unique_biparts(biparts);
			//Computes the distance and stores it based on multiple distance types
			switch (switch_value){

				case 0: //RF distance
				//	cout << "RF distance" << endl;
					dist = (b + c) / 2;
					break;
				case 1: //Euclidean distances
				//	cout << "EU dist" << endl;
					dist = sqrt(b + c);
					break;
				case 2: //Jaccard-Tanimoto distance
					//cout << "Jaccard-Tanimoto dist" << endl;
					dist = a / (a + b + c);
					break;
				case 3: //Dice distance
					dist = (2 * a) / ((2 * a) + b + c);
					break;
				case 4: //Russel-Rao distance
					dist = (a / m);
					break;
				case 5: //Sokal-Sneath
					dist = (a) / (a + (2*b) + (2 * c));
					break;
				case 6: //Ochai
					dist = a / (sqrt((a + b) * (a + c)));
					break;
				case 7: //Forbes
					dist = (a * m) / ((a + b) * (a + c));
					break;
				case 8: //Mountford 
					dist = a / (.5 * ((a*b) + (a*c)) + (b*c));
					break;
				case 9: //Sorgenfrei
					dist = (a* a) / ((a+b) * (a + c));
					break;
				case 10://Tarwid
					dist = ((m * a) - (a + b) * (a + c)) / ((m * a) + (a + b) * (a + c));
					break;
				case 11://Johnson
					dist = (a / (a + b)) + (a / ( a + c));
					break;
				case 12://Driver-Kroeber
					dist = (a/2)*((1/(a+b)) + (1/(a+c)));
					break;
				case 13://Fager-McGowan
					dist = (a / sqrt((a+b)*(a+c)))- (max(a+b,a+c)/2);
					break;
				case 14://Gilbert & Wells 
					dist = log(a) - log(m) - log((a+b)/m) - log((a + c)/m);
					break;
				case 15://Lance-Williams
					dist = (b+c)/(2*a + b + c);
					break;
				default: //No proper distance measure given
					cout << "Unknown Distance measure given.";
					break;
			}
			distances[i][j] = dist;
			distances[j][i] = dist;
		}
	}
	/*Prints the distances (for various testing purposes)
	for(unsigned int i = 0; i < distances.size(); i++){//for each tree
		cout << "Tree : " << std::setw(2) << i << ": ";
		for (unsigned int k = 0; k < i; k++){//tabs white space
			cout << "  ";
		}	
		for (unsigned int j = i; j < distances.size(); j++){//for each other tree
			cout << distances[i][j] << " ";
		}
		cout << endl;
	}*/
	return distances;
}

//The if statement which converts distance strings to integers necessary for later switch
unsigned int distance_switch(string measure){
	//Holds the string converted to a switch value
	unsigned int switch_value;

	if (measure == "rf" || measure == "RF" || measure == "Rf"){
		switch_value = 0;
	}
	else if (measure == "eu" || measure == "EU" || measure == "Eu" || measure == "euclidean"
			|| measure == "Euclidean"){
		switch_value = 1;
	}
	else if (measure == "j-t" || measure == "jaccard-tanimoto"){
		switch_value = 2;
	}
	else if (measure == "dice" || measure == "Dice"){
		switch_value = 3;
	}
	else if (measure == "r-r" || measure == "russel-rao"){
		switch_value = 4;
	}
	else if (measure == "sokal-sneath"){
		switch_value = 5;
	}
	else if (measure == "ochai"){
		switch_value = 6;
	}
	else if (measure == "forbes"){
		switch_value = 7;
	}
	else if (measure == "mountford"){
		switch_value = 8;
	}
	else if (measure == "sorgenfrei"){
		switch_value = 9;
	}
	else if (measure == "tarwid"){
		switch_value = 10;
	}
	else if (measure == "johnson"){
		switch_value = 11;
	}
	else if (measure == "driver-kroeber"){
		switch_value = 12;
	}
	else if (measure == "fager-mcgowan"){
		switch_value = 13;
	}
	else if (measure == "gilbertwells"){
		switch_value = 14;
	}
	else if (measure == "lance-williams"){
		switch_value = 15;
	}
	else if (measure == "quartet" || measure == "qd"){
		switch_value = 20;
	}
	else if (measure == "conflicting-quartet" || measure == "cqd"){
		switch_value = 21;
	}
	else if (measure == "editg" || measure == "greedy-edit"){
		switch_value = 22;
	}
	else if (measure == "editt" || measure == "total-edit"){
		switch_value = 23;
	}
	else if (measure == "editm" || measure == "minimum-edit"){
		switch_value = 24;
	}
	else if (measure == "editm" || measure == "minimum-edit"){
		switch_value = 25;
	}
	return switch_value;
}

//Computes the distance matrix for the given input and measure, input is a treeset, for when order isn't especially important
vector <vector < unsigned int> > compute_distances(set < unsigned int > treeset, string measure){
	//Return Value
	vector < vector <unsigned int> > distances;
	
	unsigned int switch_value = distance_switch(measure);
	if(switch_value < 20){
		vector < vector < unsigned int> > biparts;
	
		for(std::set<unsigned int>::iterator pos = treeset.begin(); pos != treeset.end(); ++pos){//for each tree
			biparts.push_back(::biparttable.inverted_index.at(*pos));
		}	
		distances = bipart_distances(biparts, switch_value);
	}

	else {
		//Quartet and Edit distances
		distanceWrapper(treeset, switch_value);
	}
	return distances;

}

//Computes the distance matrix for the given input and measure, input is a vector of trees, for when the given order
//is important (especially for clustering and clustering visualization).
vector < vector < unsigned int > > compute_distances(vector < unsigned int > treevect, string measure){
	//Return Value
	vector < vector <unsigned int> > distances;
	
	unsigned int switch_value = distance_switch(measure);

	if(switch_value < 20){
	vector < vector < unsigned int> > biparts;	
	
		for(unsigned int i = 0; i < treevect.size(); i++){//for each tree
			biparts.push_back(::biparttable.inverted_index.at(treevect[i]));
		}
		
		distances = bipart_distances(biparts, switch_value);
	}




	return distances;
}	

//Various tests that have been used for the distance functions
void TestDist(){
	set < unsigned int > test_trees;
	test_trees.insert(1);
	test_trees.insert(2);
	test_trees.insert(3);
	test_trees.insert(4);
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
