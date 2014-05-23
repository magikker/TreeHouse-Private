#include "distance.h"
using namespace std;



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
//vector < vector <unsigned int> > bipart_distances(vector < vector <unsigned int> > biparts, unsigned int switch_value){
vector < vector <float> > bipart_distances(vector < vector <unsigned int> > biparts, unsigned int switch_value){
	//Return Value
	vector< vector <float> > distances;

	distances.resize(biparts.size(), vector< float >(biparts.size(), 0));

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
			float dist;
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
					dist = (b + c) / 2.0;
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
				case 27://Similarity
					dist = a;
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
	else if (measure == "taxasim" || measure == "taxa-similarity"){
		switch_value = 26;
	}
	else if (measure == "sim" || measure == "similarity"){
		switch_value = 27;
	}
	else if (measure == "cp" || measure == "clade-potential"){
		switch_value = 28;
	}
	return switch_value;
}

//Computes the distance matrix for the given input and measure, input is a treeset, for when order isn't especially important
vector <vector < float> > compute_distances(set < unsigned int > treeset, string measure){
	//Return Value
	vector < vector <float> > distances;
	
	unsigned int switch_value = distance_switch(measure);
	if(switch_value < 20 || switch_value == 27){
		vector < vector < unsigned int> > biparts;
		for(std::set<unsigned int>::iterator pos = treeset.begin(); pos != treeset.end(); ++pos){//for each tree
			biparts.push_back(::biparttable.inverted_index.at(*pos));
		}	
		distances = bipart_distances(biparts, switch_value);
	}
	else if(switch_value == 26)
		distances = taxaSimilarityMatrix(treeset);
	else if(switch_value == 28){
		vector < vector <float> > taxaInCommon;
		taxaInCommon = taxaSimilarityMatrix(treeset);
		distances = taxaInCommon;
		vector < vector <float> > CladesInCommon;
		vector < vector < unsigned int> > biparts;
		for(std::set<unsigned int>::iterator pos = treeset.begin(); pos != treeset.end(); ++pos){//for each tree
			biparts.push_back(::biparttable.inverted_index.at(*pos));
		}	
		CladesInCommon = bipart_distances(biparts, 27);
	
		cout << taxaInCommon.size() << ":" << CladesInCommon.size() << endl;
	
		for(int i = 0; i < taxaInCommon.size(); i++){
			for (int j = 0; j < taxaInCommon[i].size(); j++){
				if (taxaInCommon[i][j] > 1){
					cout << "They have "<< taxaInCommon[i][j] <<" taxa in common and "<< CladesInCommon[i][j]- taxaInCommon[i][j] <<" clades in common "<< endl;
					distances[i][j] = ( (CladesInCommon[i][j]- taxaInCommon[i][j]) / (taxaInCommon[i][j] -1) ) ;
				}
				else{
					distances[i][j] = 0;
				}
			}
		}
	}	
	else {
		//Quartet and Edit distances
		distances = distanceWrapper(treeset, switch_value);
	}
	return distances;
}

//Computes the distance matrix for the given input and measure, input is a vector of trees, for when the given order
//is important (especially for clustering and clustering visualization).
vector < vector < float > > compute_distances(vector < unsigned int > treevect, string measure){
	//Return Value
	vector < vector <float> > distances;
	
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

















int sumPair(pair<int, int> x){
  return x.first + x.second;
}


vector<vector<float>> distanceWrapper(set<unsigned int> in, int mode){
  vector< unsigned int> trees(in.begin(), in.end());
  cout << "trees size is: " << trees.size() << endl;
  cout << "AT SWITCH, MODE IS " << mode << endl;
  vector<vector<float>> retVal;
  vector<float> containers;
  for(int i = 0; i < trees.size(); i++){
	retVal.push_back(containers);
	for(int j = 0; j < trees.size(); j++){
		retVal.at(i).push_back(0);
		}
	}
 

  for(int i = 0; i < trees.size(); i++){

 // int zerosToPush = 1;
  //for(unsigned int i = 0; i < trees.size() - 1; i++){
//	for(int z = zerosToPush; z > 0; z--){
//		retVal.at(i).push_back(0); //initial 0 of self to self comparison
//		}	

	for(unsigned int j = i + 1; j < trees.size(); j++){
		if(mode==20){ //quartet distance
			unsigned int distance = quartet_distance(trees.at(i), trees.at(j));		
			//retVal.at(i).push_back( quartet_distance(trees.at(i), trees.at(j)) );
			retVal.at(i).at(j) = distance;
			retVal.at(j).at(i) = distance;
			}

		else if(mode==21){ //conflicting quartet distance
			unsigned int distance = conflictingQuartetDistance(trees.at(i), trees.at(j));		
			retVal.at(i).at(j) = distance;
			retVal.at(j).at(i) = distance;
			}
		else if(mode==22){ //edit distance greedy
			unsigned int distance = edit_distance_greedy(trees.at(i), trees.at(j));		
			//retVal.at(i).push_back( edit_distance_greedy(trees.at(i), trees.at(j)) );
			retVal.at(i).at(j) = distance;		
			retVal.at(j).at(i) = distance;	
			}
		else if(mode==23){ //edit distance total
			unsigned int distance = edit_distance_total(trees.at(i), trees.at(j));		
			//retVal.at(i).push_back( edit_distance_total(trees.at(i), trees.at(j)) );
			retVal.at(i).at(j) = distance;
			retVal.at(j).at(i) = distance;			
}

		else if(mode==24){ //edit distance min
			unsigned int distance = edit_distance_minimum(trees.at(i), trees.at(j));		
			//retVal.at(i).push_back( edit_distance_minimum(trees.at(i), trees.at(j)) );
			retVal.at(i).at(j) = distance;	
			retVal.at(j).at(i) = distance;			}
		
		}
	}

//now, we have the upper triangle of a matrix. Make it square:

for(int i = 0; i < retVal.at(0).size(); i++){
	//cout << "retval's size at: " << i << " is: " << retVal.at(0).size() << endl;
	}
  return retVal;

}

void testCalculateC(){
  cout << "calculating left and right pairs:\n";
  string nw = "((A,B), (B,C))";
  nw = "( (A,B),C)";
 // string nw2 = 
  //nw = "(B,C)";
  unsigned int total = 0;
  pair<int,int> lr = numLeftRight(nw, total);
  cout << lr.first << ", " << lr.second << endl;
  cout << "total is: " << total << endl;
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
double edit_distance_total(unsigned int tree1, unsigned int tree2){
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
		//make sure that one bitstring isn't inverted
		if(XORdistance > biparttable.lm.size()/2){
			XORdistance -= biparttable.lm.size()/2;
			}
		total+=XORdistance;
		}
	}
  //note- we add to total rfDistance^2 number of times
  // since we want total to reflect the rf distance * average distance between any two given bipartitions,
  // we need to divide total by rf.size
 return (double)total / (double)rf.size();
}



unsigned int edit_distance_greedy(unsigned int tree1, unsigned int tree2){
  //NOTE- CURRENTLY ONLY WORKS ON HOMOGENEOUS DATA!

//this is average edit distance * rf distance
//in other words, an all-to-all bipartition edit distance

 unsigned int total = 0; 

 pair<set<unsigned int>, set<unsigned int>> rfs = rfDistanceSet(tree1, tree2);
 set<unsigned int> rf, rf2;
 rf = rfs.first; rf2 = rfs.second;

 //we need to make sure the rf distance is non-0 or else the program will crash!
 if(rf.size()!=0){
	 //now, calculate the XOR matrix
	 vector< vector <int> > edit = editDistanceMatrixVector(rf, rf2);
	
	  //now we have the matrix of edit distances. 
	  //keep a sorted list of all the indeces we have to explore

	set<unsigned int> indices;
	  for(unsigned int i = 0; i < edit.at(0).size(); i++){
		indices.insert(i);
		}
	  //go through the edit vector for each bipartition
	  //find the min, then remove that index from indices
	  //also add it to the total


	  int min = 0;
	  for(unsigned int i = 0; i < edit.size(); i++){
		min = edit.at(i).at(* (indices.begin()) );
		int index = *( indices.begin() ); //holds where we find the min from this row	
		//cout << "min is: " << min << " index is: " << index << endl;	
		//printVectorCompact(edit.at(i));
	
		for(set<unsigned int>::iterator j = indices.begin(); j!=indices.end(); j++){
			if(edit.at(i).at(*j) < min) {
				min = edit.at(i).at(*j);
				index = *j;
				}
			}
		//now we've gotten the min edit distance from bipartition i. 
		//add its min to to the total and remove index from indices
		total+= min;
		set<unsigned int>::iterator it = indices.find(index);
		if(it==indices.end()){
			cout << "edit_distance_greedy error: index not found\n";
			cout << "index is: " << index << endl;
			printSetCompactTwo(indices);
			printVectorCompact(edit.at(i));
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


double edit_distance_average(unsigned int tree1, unsigned int tree2){
  return edit_distance_total(tree1, tree2) /  (double)CRFDistance(tree1, tree2);
}


unsigned int edit_distance_minimum(unsigned int tree1, unsigned int tree2){

 pair<set<unsigned int>, set<unsigned int>> rfs = rfDistanceSet(tree1, tree2);
 set<unsigned int> rf, rf2;
 rf = rfs.first; rf2 = rfs.second;


  if(rf.size()!=0){

  	vector< set <int> > edit = editDistanceMatrixSet(rf, rf2);
	 //now we have the a vector of sets for each bipartition. Take the min from each set (i.e. first element)
	  unsigned int acc = 0;
	  for(vector<set<int>>::iterator i = edit.begin(); i!=edit.end(); i++){
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

double edit_distance_minimum_coverage(unsigned int tree1, unsigned int tree2){
//returns the minimum edit distance AND what proportion of the columns were covered

	pair<set<unsigned int>, set<unsigned int>> rfs = rfDistanceSet(tree1, tree2);
	set<unsigned int> rf, rf2;
	rf = rfs.first; rf2 = rfs.second;
	set<unsigned int> coverSet; //keeps track of which biparts in rf2 have had their edit distance added to total;
	  if(rf.size()!=0 && rf2.size()!=0){
		 //now, calculate the XOR matrix
		 vector< vector <int> > edit = editDistanceMatrixVector(rf, rf2);

		 //now we have the a vector of sets for each bipartition. Take the min from each set (i.e. first element)
		  for(unsigned int i = 0; i < edit.size(); i++){
			//printVectorCompact(edit.at(i));
			int min = edit.at(i).at(0);
			int minIndex = 0;		
			for(unsigned int j = 0; j < edit.at(i).size(); j++){
				if(edit.at(i).at(j) < min){
					minIndex = j;
					min = edit.at(i).at(j);
					}
				}
				coverSet.insert(minIndex);
			}
		return (double)coverSet.size() / (double)rf2.size();
		}

	  else{
		return 0;
		}
}

float CRFStringDistance(string tree1, string tree2){

	vector< vector <int > > bipartsone = compute_bitstrings_h(tree1);
	
	set<boost::dynamic_bitset<> > bitsone;
	for(unsigned int i = 0; i < bipartsone.size(); i++){
		boost::dynamic_bitset<> bi(::biparttable.lm.size());
		for(unsigned int j = 0; j < bipartsone[i].size(); j++){
			bi[j] = 1;
		}
		bitsone.insert(bi);
	}
	
	vector< vector <int > > bipartstwo = compute_bitstrings_h(tree2); 
	set<boost::dynamic_bitset<> > bitstwo;
	for(unsigned int i = 0; i < bipartstwo.size(); i++){
		boost::dynamic_bitset<> bi(::biparttable.lm.size());
		for(unsigned int j = 0; j < bipartstwo[i].size(); j++){
			bi[j] = 1;
		}
		bitstwo.insert(bi);
	}

	vector<boost::dynamic_bitset<> > rf, rf2;	
	set_difference(bitsone.begin(), bitsone.end(), bitstwo.begin(), bitstwo.end(), inserter(rf, rf.begin()));
	set_difference(bitstwo.begin(), bitstwo.end(), bitsone.begin(), bitsone.end(), inserter(rf2, rf2.begin()));
	
	return (rf.size() + rf2.size()) / 2.0;
	
} 

float CRFDistance(int tree1, int tree2){
  set<boost::dynamic_bitset<> > bitsone, bitstwo, rf, rf2;
  
  typedef std::map< boost::dynamic_bitset<>, TreeSet >::iterator clade_it_type;

  for(clade_it_type iter = biparttable.MapBenchMarks[2]; iter != biparttable.MapBenchMarks[biparttable.lm.size()]; iter++) {
	if (biparttable.contains_tree(iter, tree1)){
		bitsone.insert(iter->first);
	}
	if (biparttable.contains_tree(iter, tree2)){
		bitstwo.insert(iter->first);
	}
  }
  
  set_difference(bitsone.begin(), bitsone.end(), bitstwo.begin(), bitstwo.end(), inserter(rf, rf.begin()));
  set_difference(bitstwo.begin(), bitstwo.end(), bitsone.begin(), bitsone.end(), inserter(rf2, rf2.begin()));
  return (rf.size() + rf2.size()) / 2.0;
}


// I need to count the taxa in each tree to make sure they match. 
float BipartRFStringDistance(string tree1, string tree2){
	//I'm flipping bitstrings as needed so that the first bit in each is a one. This makes it such that all dups are accounted for and I can do a proper bipartition count across the two trees. 

	vector< vector <int > > bipartsone = compute_bitstrings_h(tree1);	
	set<boost::dynamic_bitset<> > bitsone;
	for(unsigned int i = 0; i < bipartsone.size(); i++){
		boost::dynamic_bitset<> bi(::biparttable.lm.size());
		for(unsigned int j = 0; j < bipartsone[i].size(); j++){
			bi[j] = 1;
		}
		
		if ( bi[0] == 0 ){
			bitsone.insert(bi.flip());
		}
		else{
			bitsone.insert(bi);
		}	
	}
	
	vector< vector <int > > bipartstwo = compute_bitstrings_h(tree2); 
	set<boost::dynamic_bitset<> > bitstwo;
	for(unsigned int i = 0; i < bipartstwo.size(); i++){
		boost::dynamic_bitset<> bi(::biparttable.lm.size());
		for(unsigned int j = 0; j < bipartstwo[i].size(); j++){
			bi[j] = 1;
		}
		if ( bi[0] == 0 ){
			bitstwo.insert(bi.flip());
		}
		else{
			bitstwo.insert(bi);
		}
	}

	vector<boost::dynamic_bitset<> > rf, rf2;	
	set_difference(bitsone.begin(), bitsone.end(), bitstwo.begin(), bitstwo.end(), inserter(rf, rf.begin()));
	set_difference(bitstwo.begin(), bitstwo.end(), bitsone.begin(), bitsone.end(), inserter(rf2, rf2.begin()));
	
	return (rf.size() + rf2.size()) / 2.0;
} 

//needs to be created. 
//float BipartRFDistance(int tree1, int tree2){
//  unsigned long total = 0;
//  vector<unsigned int> rf, rf2;
//  //rf is the one sided difference from t1 to t2. rf2 is from t2 to t1
//  vector<unsigned int> t1 = biparttable.inverted_index.at(tree1);
//  vector<unsigned int> t2 = biparttable.inverted_index.at(tree2);
//  set_difference(t1.begin(), t1.end(), t2.begin(), t2.end(), inserter(rf, rf.begin()));
//  set_difference(t2.begin(), t2.end(), t1.begin(), t1.end(), inserter(rf2, rf2.begin()));
//  return (rf.size() + rf2.size()) / 4.0;
//}

float BipartRFDistance(int tree1, int tree2){
  set<boost::dynamic_bitset<> > bitsone, bitstwo, rf, rf2;
  
  typedef std::map< boost::dynamic_bitset<>, TreeSet >::iterator clade_it_type;

  for(clade_it_type iter = biparttable.MapBenchMarks[2]; iter != biparttable.MapBenchMarks[biparttable.lm.size()]; iter++) {
	if (biparttable.contains_tree(iter, tree1)){
		boost::dynamic_bitset<> bi(::biparttable.lm.size());
		
		bi = iter->first;
		
		if(bi[0] == 1){
			bitsone.insert(bi);
		}
		else{
			bitsone.insert(bi.flip());
		}
	}
	if (biparttable.contains_tree(iter, tree2)){
		boost::dynamic_bitset<> bi(::biparttable.lm.size());
	
		bi = iter->first;
		if(bi[0] == 1){
			bitstwo.insert(bi);
		}
		else{
			bitstwo.insert(bi.flip());
		}
	}
  }
  set_difference(bitsone.begin(), bitsone.end(), bitstwo.begin(), bitstwo.end(), inserter(rf, rf.begin()));
  set_difference(bitstwo.begin(), bitstwo.end(), bitsone.begin(), bitsone.end(), inserter(rf2, rf2.begin()));
  return (rf.size() + rf2.size()) / 2.0;
}


void computeSimMatrixStats(vector<vector<float>> matrix){

	unsigned int nosim = 0;
	unsigned int onesim = 0;
	unsigned int twosim = 0;
	unsigned int threesim = 0;
	unsigned int fourormoresim = 0;
	unsigned int nonzerocells = 0;
	unsigned int cellsfourormore = 0;
	
	set<unsigned int> nosimtrees;
	set<unsigned int> onetothreesimtrees;
	set<unsigned int> fourormoresimtrees;
	
	
	for(size_t i = 0; i < matrix.size(); i++){
			if (*std::max_element(matrix[i].begin(),matrix[i].end()) == 0){
				nosim +=1;
				nosimtrees.insert(i);
			}
			if (*std::max_element(matrix[i].begin(),matrix[i].end()) == 1){
				onesim +=1;
				onetothreesimtrees.insert(i);
			}
			if (*std::max_element(matrix[i].begin(),matrix[i].end()) == 2){
				twosim +=1;
				onetothreesimtrees.insert(i);
			}
			if (*std::max_element(matrix[i].begin(),matrix[i].end()) == 3){
				threesim +=1;
				onetothreesimtrees.insert(i);
			}
			if (*std::max_element(matrix[i].begin(),matrix[i].end()) > 3){
				fourormoresim +=1;
				fourormoresimtrees.insert(i);
			}
			std::cout << "The largest element is "  << *std::max_element(matrix[i].begin(),matrix[i].end()) << '\n';
		
		for (size_t j = 0; j < matrix[i].size(); j++){
			if (matrix[i][j] > 0){
				nonzerocells +=1;
			}
			if (matrix[i][j] > 3){
				cellsfourormore +=1;
			}
		}
	}
	cout << "Of the trees " << nosim << " of them have no taxa in common with any other tree" << endl;
	cout << "Of the trees " << onesim << " of them have one taxon in common with another other tree" << endl;
	cout << "Of the trees " << twosim << " of them have two taxa in common with another tree" << endl;
	cout << "Of the trees " << threesim << " of them have three taxa in common with another tree" << endl;
	cout << "Of the trees " << fourormoresim << " of them have four or more taxa in common with another tree" << endl;
	cout << "There are " << nonzerocells << " non zero cells" << endl;
	cout << "There are " << cellsfourormore << " cells in the matrix with a value of four or more" << endl;
	write_tre_file(nosimtrees, "temp/nosimtrees");
	write_tre_file(onetothreesimtrees, "temp/onetothreesimtrees");
	write_tre_file(fourormoresimtrees, "temp/fourormoresimtrees");
}



//GRB
unsigned int taxaSimilarity(unsigned int tree1, unsigned int tree2){
	boost::dynamic_bitset<> temp(biparttable.lm.size());
	temp = ::biparttable.taxa_in_trees[tree1] & ::biparttable.taxa_in_trees[tree2];
	return temp.count();
}

 vector<vector<float>> taxaSimilarityMatrix(set<unsigned int> in){
   vector<vector<float>> retVal;
  
   vector< unsigned int> trees(in.begin(), in.end());

  
  vector<float> containers;
  for(unsigned int i = 0; i < trees.size(); i++){
	retVal.push_back(containers);
	for(unsigned int j = 0; j < trees.size(); j++){
		retVal.at(i).push_back(0);
	}
  }
 

  for(unsigned int i = 0; i < trees.size(); i++){

 // int zerosToPush = 1;
  //for(unsigned int i = 0; i < trees.size() - 1; i++){
//	for(int z = zerosToPush; z > 0; z--){
//		retVal.at(i).push_back(0); //initial 0 of self to self comparison
//		}	

	for(unsigned int j = i + 1; j < trees.size(); j++){
		//if(mode==20){ //quartet distance
			float distance = taxaSimilarity(trees.at(i), trees.at(j));		
			//retVal.at(i).push_back( quartet_distance(trees.at(i), trees.at(j)) );
			retVal.at(i).at(j) = distance;
			retVal.at(j).at(i) = distance;
			}
	}
	computeSimMatrixStats(retVal);
	return retVal;
}





vector<set<int>> editDistanceMatrixSet(set<unsigned int> rf, set<unsigned int> rf2){
  //calculate the all-to-all edit distance between bipartitions in rf and rf2
  vector< set<int> > edit;
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
		if(XORdistance > biparttable.lm.size()/2){
			XORdistance -= biparttable.lm.size()/2;
			}
		h.insert(XORdistance);
		}
	edit.push_back(h); //add new row to the edit matrix
	}
  return edit;
}

vector<vector<int>> editDistanceMatrixVector(set<unsigned int> rf, set<unsigned int> rf2){	 

  vector< vector <int> > edit;
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
			if(XORdistance > biparttable.lm.size()/2){
			XORdistance -= biparttable.lm.size()/2;
			}
			h.push_back(XORdistance);
			}
		edit.push_back(h); //add new row to the edit matrix
		}
	return edit;
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
 testIsBifurcating();
//testCalculateC();
//testEdit();
//testNewick();
//testDepthVariance();


}


void testEdit(){
 for(int i = 250; i < 300; i++){
	for(int j = 1000; j < 1050; j++){
		unsigned int rf, min;
		cout << "\nTesting edit for " << i << " " << j << ":\n";
                cout << "	rf distance = " << (rf = CRFDistance(i, j));
		cout << ", HDG = " << edit_distance_greedy(i,j);
	//edit_distance_greedy(i,j);
	//edit_distance_total(i,j);
	//edit_distance_total(i,j);
		cout << ", HDT = " << edit_distance_total(i,j);
		cout << ", HDA = " << edit_distance_average(i,j);
		cout << ", HDM = " << (min = edit_distance_minimum(i,j));
		cout << "\n HDM ave = " << (double)min/(double)rf;
		cout << "HDMC: " << edit_distance_minimum_coverage(i,j);
		}
	//cout << "done with " << i << endl;
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
