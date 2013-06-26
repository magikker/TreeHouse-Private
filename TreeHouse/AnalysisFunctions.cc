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

/*
//Returns the difference based on shared bipartitions from each tree to each other (RF distance)
vector < vector < unsigned int > > compute_rf_distances(vector < vector < unsigned int> > biparts){
vector < vector < unsigned int> > distances;
distances.resize(biparts.size(), vector<unsigned int>(biparts.size(),0));

for(unsigned int i = 0; i < biparts.size() - 1; i++){//for each tree's bipartitions
for (unsigned int j = 1; j < biparts.size(); j++){//for all others
vector<unsigned int> temp;
int diff1;
int diff2;
int dist;
//intersection contain all shared bipartitions
std::set_intersection(biparts[i].begin(),biparts[i].end(), 
biparts[j].begin(), biparts[j].end(),
std::inserter(temp, temp.end()));
//Stores the differences for each
diff1 = biparts[i].size() - temp.size();
diff2 = biparts[j].size() - temp.size();
//computes the actual distance the stores it
dist = (diff1 + diff2) / 2;
distances[i][j] = dist;
distances[j][i] = dist;
}
}
return distances;
}

vector < vector < unsigned int > > compute_euclidean_distances(vector < vector < unsigned int> > biparts){
vector < vector <unsigned int > > distances;
distances.resize(biparts.size(), vector<unsigned int>(biparts.size(),0));

for(unsigned int i = 0; i < biparts.size() - 1; i++){//for each tree's bipartitions
for (unsigned int j = 1; j < biparts.size() - 1; i++){//for all others
vector<unsigned int> temp;
int diff; 
int diff2;
int dist;
//Intersection contains all shared bipartitions
std::set_intersection(biparts[i].begin(),biparts[i].end(),
biparts[j].begin(), biparts[j].end(),
std::inserter(temp,temp.end()));
//Stores the differences for each
diff1 = biparts[i].size() - temp.size();
diff2 = biparts[j].size() - temp.size();
//Computes the distnace and stores it
dist = sqrt(diff1 + diff2);
distances[i][j] = dist;
distances[j][i] = dist;
}
}
}
*/

//Computes various distance measures based on the bipartitions
vector < vector < unsigned int > > compute_bipart_distances(set <unsigned int> treeset, string measure){
	//Return Value
	vector < vector < unsigned int > > distances;
	//Holds bipartitions
	vector < vector < unsigned int> > biparts;
	//Resizes to hold the proper number of elements
	//Arbitrary starting value
	unsigned int switch_value = 20;
	 int m; //# unique biparitions

	m = unique_biparts(treeset);

	//Computes the bipartitions
	for(std::set<unsigned int>::iterator pos = treeset.begin(); pos!= treeset.end(); ++pos){//for each tree
		biparts.push_back(biparts_in_tree(*pos));
	}


	distances.resize(biparts.size(), vector< unsigned int >(biparts.size(), 0));

	//To set the switch since strings are intuitive to us but not switch statements
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

	for(unsigned int i = 0; i < biparts.size() - 1; i++){//for each tree's bipartitions
		for (unsigned int j = 1; j < biparts.size(); j++){//for all others
			//a, b, & c values come from Suzanne Matthews Dissertation
			//and the method of computing distances from bipartitions found there
			int a; //Bipartitions in both trees
			vector<unsigned int> temp;
			int b; //# Bipartions in first tree but not second
			int c; //# Bipartitions in second tree but not first
			int dist;
			//Intersection contains all shared bipartitions
			std::set_intersection(biparts[i].begin(),biparts[i].end(),
					biparts[j].begin(),biparts[j].end(),
					std::inserter(temp, temp.end()));
			//Stores the differences for each
			a = temp.size();
			b = biparts[i].size() - temp.size();
			c = biparts[j].size() - temp.size();
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
				default: //No proper distance measure given
					cout << "Unknown Distance measure given.";
					break;
			}
			distances[i][j] = dist;
			distances[j][i] = dist;
		}
	}

	//Prints the distances (for various testing purposes)
	//for(unsigned int i = 0; i < distances.size(); i++){//for each tree
	//	cout << "Tree : " << std::setw(2) << i << ": ";
	//	for (unsigned int k = 0; k < i; k++){//tabs white space
	//		cout << "  ";
	//	}	
	//	for (unsigned int j = i; j < distances.size(); j++){//for each other tree
	//		cout << distances[i][j] << " ";
	//	}
	//	cout << endl;
//	}
	return distances;
}

//Various tests that have been used for the distance functions
void TestDist(){
	set < unsigned int > test_trees;
	test_trees.insert(1);
	test_trees.insert(2);
	test_trees.insert(3);
	test_trees.insert(4);

//	cout << "rf distance" << endl;
//	compute_bipart_distances(test_biparts, "rf");
//	cout << "euclidean distance" << endl;
//	compute_bipart_distances(test_biparts, "eu");
}


//Returns the average distance for each tree to the cluster containing it
vector<float> diff_to_own_clusters(vector< set <unsigned int> > inputclusts, vector < vector <unsigned int> > distances){
	unsigned int offset = 0;
	vector< float > diff;
	for(unsigned int i = 0; i < inputclusts.size(); i++){//For each cluster
		float tempsum = 0;
		float tempaverage;
		if(i > 0){
			//offset used for tracking change from cluster to cluster
			offset += inputclusts[i-1].size();
		}
		for(unsigned int j = 0 + offset; j < offset + 1; j++){//For first tree in the cluster
			//Temporarily stores sum for computing average
			for(unsigned int k = 1 + offset; k < inputclusts[i].size() + offset; k++){//For each othe tree in cluster
				tempsum += distances[j][k];
			}
			//average difference from tree to all other trees in cluster
		}

		tempaverage = tempsum / inputclusts[i].size();

		diff.push_back(tempaverage);
	}
	return diff;
}



//Returns the average distance for each tree to the cluster with the smallest average distance from the tree
vector<float> neighboring_cluster_diffs(vector  < set < unsigned int> > inputclusts, vector < vector < unsigned int> > distances){
	unsigned int offset = 0;
	vector < float > diffs;
	diffs.resize(distances[0].size(), 100);
	for(unsigned int i = 0; i < inputclusts.size(); i++){//For each cluster
		if(i > 0){
			//offset used for tracking change from cluster to cluster
			offset += inputclusts[i-1].size();
		}
		for(unsigned int j = offset; j < inputclusts[i].size() + offset; j++){//For each tree in cluster
			for(unsigned int k = 0; k < inputclusts.size(); k++){//for each cluster
				//not interested in the cluster containing the tree

				float tempsum = 0;
				//koffset used for tracking place in second cluster
				unsigned int koffset;

				if(k==0){
					koffset = 0;
				}
				else{
					koffset += inputclusts[k-1].size();
				}
				if(k == i){
					continue;
				}
				for(unsigned int l = koffset; l < inputclusts[k].size() + koffset; l++){//For each tree in cluster
					tempsum += distances[j][l];
				}

				//computes and stores the minimum average distance
				//to another cluster
				float tempaverage;
				tempaverage = tempsum /inputclusts[i].size();
				if(tempaverage < diffs[j]){
					diffs.insert(diffs.begin() + (j), tempaverage);
				}
			}
		}
	}
	return diffs;
}



//Returns the silhouette width of each tree or cluster in the input as a vector of floats
vector<float> silhouette (vector< set <unsigned int>  > inputclusts, string dist_type){
	//Return value (by tree)
	vector<float> swidth;
	//Return value (by cluster)
	vector<float> cwidth;
	//Stores the set of all trees in the clusters
	set<unsigned int> treeset;
	//Stores the average distance for each tree to its own cluster
	vector< float > diff_within_cluster;	
	//Stores the average distance for each tree to its neighboring cluster
	vector< float > neighboring_cluster;
	//Stores all of the distances
	vector < vector < unsigned int > > distances;



	//Obtains the treeset
	for(unsigned int i = 0; i < inputclusts.size(); i++){//for each cluster
		for(std::set<unsigned int>::iterator pos = inputclusts[i].begin(); pos!=inputclusts[i].end(); ++pos){//for each tree
		treeset.insert(*pos);	
		}
	}
	//Finds and stores the distances
	distances = compute_bipart_distances(treeset, dist_type);
	//For checking that silhouette properly computes the silhouette width

	//	for(unsigned int i = 0; i < distances.size(); i++){//for each tree
	//		cout << "Tree : " << i << ": "; 
	//		for (unsigned int j = i + 1; j < distances.size(); j++){//for each other tree
	//			cout << distances[i][j] << " ";
	//		}
	//		cout << endl;
	//	}

	//Computes diff_within_cluster
	diff_within_cluster = diff_to_own_clusters(inputclusts, distances);

	//	cout <<"diff_within_cluster " << endl;
	//	for (unsigned int i = 0; i < diff_within_cluster.size(); i++){
	//		cout << diff_within_cluster[i] << " ";
	//	}

	//Computes the neighboring_clusters
	neighboring_cluster = neighboring_cluster_diffs(inputclusts, distances);
	//	cout <<"neighboring_cluster " << endl;
	//	for (unsigned int i = 0; i < neighboring_cluster.size(); i++){
	//		cout << neighboring_cluster[i] << " ";
	//		}
	//		cout << endl;

	//Computes the actual silhouette width
	unsigned int soffset = 0;
	for(unsigned int i = 0; i < inputclusts.size(); i++){//for each cluster
		if(i != 0){
			soffset += inputclusts[i-1].size();
		}
		for(unsigned int j = 0 + soffset; j < inputclusts[i].size() + soffset; j++){//for each tree in cluster
			float ai = diff_within_cluster[i];
			float bi = neighboring_cluster[j];
			//			cout << "ai is " << ai << endl;
			//			cout << "bi is " << bi << endl;
			float si;
			//Prevents dividing by zero
			if(ai == bi){
				si = 0;
			}
			else{
				si = (bi - ai) /max(ai, bi);
			}
			swidth.insert(swidth.begin() + j, si);
		}
	}
	unsigned int ioffset = 0;

	//Stores the cluster average silhouette widths
	for(unsigned int i = 0; i < inputclusts.size(); i++){//for each cluster
		if (i !=0){
			ioffset += inputclusts[i - 1].size();
		}
		float caverage = 0;
		for (unsigned int j = 0 + ioffset; j < inputclusts[i].size() + ioffset; j++){//for each tree in cluster
			caverage += swidth[j];
		}
		cwidth.insert(cwidth.begin() + i, caverage / inputclusts[i].size());
	}

	//prints the clusters and their silhouette widths (by cluster)
	for(unsigned int i = 0; i < cwidth.size(); i++){//for each cluster width
		cout <<"Cluster: " << i << " has a silhouette of " << cwidth[i] << endl;
	}

	//prints the tree and their silhouette widths (by tree)
	//	unsigned int offset = 0;
	//	for(unsigned int i = 0; i < inputclusts.size(); i++){//for each cluster
	//		if(i != 0){
	//			offset += inputclusts[i - 1].size();
	//		}
	//		unsigned int j = 0;
	//		for(std::set<unsigned int>::iterator pos = inputclusts[i].begin(); pos!=inputclusts[i].end(); ++pos){//for each tree
	//			cout << "Tree: " << *pos << " has a silhouette of " << swidth[j + offset] << endl;
	//			j++;
	//		}
	//	}

	//	return swidth;

	return cwidth;
}


//Returns a pair representing the closest neighbor of a cluster, and the distances between the two clusters
pair<unsigned int, float> closest_neighbor(unsigned int clustid, set <unsigned int> checkset, vector <vector <unsigned int> > distances){
	//neighbor to clustid to be returned
	pair<unsigned int, float> retneighbor;

	float bestdist = 100;
	unsigned int best_neighbor;
	for(std::set<unsigned int>::iterator pos = checkset.begin(); pos != checkset.end(); pos++){//for each cluster in checkset
		float tempdist;

		tempdist = distances[*pos][clustid];
		if (tempdist < bestdist){
			bestdist = tempdist;
			best_neighbor = *pos;
		}
	}
	retneighbor = make_pair(best_neighbor, bestdist);

	return retneighbor;
}

//Returns a map of clusters with the two input clusters merged into a single cluster 
void  merge(unsigned int clust1, unsigned int clust2, 
		map <unsigned int, set <unsigned int > > &clusters){
	//Take the elements from clust2's set and insert them into clust1
	clusters[clust1].insert(clusters[clust2].begin(),clusters[clust2].end());
	//Remove clust2
	clusters.erase(clust2);
}

//Returns a distance matrix where the distances are recomputed for merged clusters
void recompute_distances(vector < vector < unsigned int > > &distances, map < unsigned int, set < unsigned int > > clusters, unsigned int clustid, unsigned int clustid2){

	typedef std::map < unsigned int, set <unsigned int> >::iterator it_map;
	for (it_map iterator = clusters.begin(); iterator != clusters.end(); iterator++){//for each cluster
		unsigned int current = iterator->first;
		//Average distance calculated
		float tempdist = (distances[current][clustid] + distances[current][clustid2]) / 2;

		//Stores the tempdist
		distances[current][clustid] = tempdist;
		distances[clustid][current] = tempdist;
	}
}




//Returns clusters of the input trees based on an agglomerative hierarchical method 
//(current implementation computes own distances
//could eventually take in a distance matrix or some other method)
vector < set < unsigned int > > agglo_clust (set <unsigned int> inputtrees, string dist_type){
	//The return set
	vector < set < unsigned int > > ret_clusters;
	//The map which will contain the clusters
	map <unsigned int, set < unsigned int > > clusters;
	//Clusters remaining to be compared
	set <unsigned int> remaining;
	//Contains the distances (distance matrix)	
	vector < vector < unsigned int > > distances;
	//Is used to handle finding pairs to merge more quickly (supposedly)
	//the first value of the pair is its cluster id, the second is its
	//distance/similarity measure to its parent in neighbor_stack
	stack<pair<unsigned int, float>> neighbor_stack;
	//Number of clusters we desire to end with
	unsigned int fin_num = 4;

	//Stores each tree
	unsigned int id = 0;
	for(std::set<unsigned int>::iterator pos = inputtrees.begin(); pos!=inputtrees.end(); ++pos){//for each tree
		set < unsigned int > temp;
		temp.clear();
		temp.insert(*pos);
		clusters.insert(std::make_pair(id,temp));//Creates a cluster for each tree
		remaining.insert(id);//Inserts the cluster ids into remaining
		id++;
	}


	//Computes the distances
	distances = compute_bipart_distances(inputtrees, dist_type);


	//Initialize neighbor_stack with an arbitrary element of remaining (the first one here) ((100 is arbitrary, 
	//it's meant to represent a very large distance, as there is no parent to the first element.))
	neighbor_stack.push(std::make_pair(*remaining.begin(), 100));
	remaining.erase(remaining.begin());


	//Form clusters while number of clusters is larger than the number of desired clusters
	while(clusters.size() > fin_num){
		//Find the closest neighbor of the top cluster on the stack in remaining
		pair<unsigned int, float> temp_neighbor;
		temp_neighbor = closest_neighbor(neighbor_stack.top().first,remaining,distances);


		//if temp neighbor is more similar to the current cluster than the current cluster's parent
		//then place temp neighbor on the stack and move on.
		if(temp_neighbor.second < neighbor_stack.top().second){
			neighbor_stack.push(temp_neighbor);
			remaining.erase(temp_neighbor.first);
		}
		//Else than the current cluster and it's parent on the stack are as closely related as possible
		//merge them immediately
		else{
			pair<unsigned int, float> tempcluster;
			pair<unsigned int, float> tempparent;

			tempcluster = neighbor_stack.top();
			neighbor_stack.pop();
			tempparent = neighbor_stack.top();
			neighbor_stack.pop();

			//Merges the two clusters into a single cluster
			merge(tempcluster.first, tempparent.first, clusters);
			//Adds the new cluster to the set of remaining elements
			remaining.insert(tempcluster.first);

			//Then recompute the distances matrix for the newly added cluster
			recompute_distances(distances, clusters, tempcluster.first, tempparent.first);
		}

		//If this has emptied the neighbor stack, add another element to start again.
		if (neighbor_stack.empty()){
			neighbor_stack.push(std::make_pair(*remaining.begin(),100));
			remaining.erase(remaining.begin());
		}
	}

	//Fills the ret_clusters set
	typedef std::map < unsigned int, set <unsigned int> >::iterator it_map;
	for(it_map iterator = clusters.begin(); iterator != clusters.end(); iterator++){//for each cluster
		set<unsigned int> tempset;
		tempset.clear();
		for(std::set<unsigned int>::iterator pos = iterator->second.begin(); pos != iterator->second.end(); ++pos){
			tempset.insert(*pos);
		}
		ret_clusters.push_back(tempset); //add to the return set
	}


	//Prints out the clusters map
	typedef std::map< unsigned int, set < unsigned int >>::iterator it_map;
	for(it_map iterator = clusters.begin(); iterator != clusters.end(); iterator++){
		cout << "Cluster : " << iterator->first << " Trees : " ;
		for(std::set<unsigned int>::iterator  pos = iterator->second.begin(); pos != iterator->second.end(); ++pos){
			cout << *pos << " " ;
		}
		cout << endl;
	}



	return ret_clusters;
}

//Returns the distances for each centroid to each tree (probably not the most efficient method)
void recompute_centroid_distances(vector < vector <unsigned int> > &centroid_distances, 
		vector< vector< unsigned int> > distances,
	       	unsigned int numtrees, 
		vector < set <unsigned int> > centroidcluster, vector < set <unsigned int> > oldcluster){
	
	for(unsigned int i = 0; i < centroidcluster.size(); i++){//for each centroid
		set <unsigned int> tempset = centroidcluster[i];
		set <unsigned int> oldset = oldcluster[i];
		//Used for attempting to compute only the change in centroids
		set <unsigned int> added;
		set <unsigned int> removed;
		set_difference(tempset.begin(), tempset.end(), oldset.begin(), oldset.end(),std::inserter(added, added.end()));
		set_difference(oldset.begin(), oldset.end(), tempset.begin(), tempset.end(),std::inserter(removed, removed.end()));

		for(unsigned int j = 0; j < numtrees; j++){//for each tree	
			float tempsum = 0;
			float tempdifference = 0;
			float dist = 0;
			for(std::set<unsigned int>::iterator pos = centroidcluster[i].begin(); pos != centroidcluster[i].end(); ++pos){//for each tree
				tempsum += distances[j][*pos];
			}
			//ifor(std::set<unsigned int>::iterator pos = removed.begin(); pos != removed.end(); ++pos){//for each tree removed
			//	tempdifference += distances[j][*pos];
		//	}
		//	if (added.size() > 0) {
		//	tempsum = tempsum / added.size();
			//cout << "added" << added.size();
		//	}
		//	if (removed.size() > 0) {
		//	tempdifference = tempdifference /removed.size();
		//	//cout << " removed " << removed.size() << endl;
		//	}
			
			float total = tempsum / centroidcluster[i].size();
			dist = total;
			centroid_distances[i][j] = dist;
			//cout << "dist is : " << dist << endl;

		}
	}
	//cout << "fin recompute" << endl;
}


//Returns clusters of the input trees clustered based on a k-means method
vector <set <unsigned int > > kmeans_clust(set <unsigned int> inputtrees, unsigned int k, string dist_type){
	//The return set
	vector < set < unsigned int > > ret_clusters;
	//The map which will contain the clusters from the previous loop
	vector <set < unsigned int > > oclusters;
	//The map which will contain the clusters from the new loop
	vector < set < unsigned int > > nclusters;
	//History of the centroids
//	vector< map <unsigned int, set <unsigned int> > > centroid_history;
	//Clusters remaining to be compared
	set <unsigned int> centroids;
	//Contains the distances (distance matrix)	
	vector < vector < unsigned int > > tree_distances;
	//Contains the distance from each tree to each centroid
	vector < vector < unsigned int > > centroid_distances;

	//Computes the tree-distances matrix
	tree_distances = compute_bipart_distances(inputtrees, dist_type);

	//Creates a history of centroids
//	for(unsigned int i = 0; i < k; i++){//k times
//		map <unsigned int, set <unsigned int>> tempclust;
//		centroid_history.push_back(tempclust);
//	}


	//Arbitrarily assign original centroids (shuffles a set randomly and selects the first k)

	vector < unsigned int > shuffled;
	for(std::set<unsigned int>::iterator pos = inputtrees.begin(); pos != inputtrees.end(); ++pos){//for each tree
		shuffled.push_back(*pos);
	}
	random_shuffle(shuffled.begin(), shuffled.end());
	for(unsigned int i = 0; i < k; i++){//while we do not have k clusters
		set<unsigned int> temp;
		temp.clear();
		temp.insert(shuffled[i]);
		cout << "shuffled[i] = " << shuffled[i] << endl;
		nclusters.push_back(temp);
		centroids.insert(i);
	}


	//Initializes the centroid_distance matrix
	centroid_distances.resize(centroids.size(), vector <unsigned int>(inputtrees.size(), 0));
	for(unsigned int i = 0; i < centroid_distances.size(); i++){//for each centroid
		for(unsigned int j = 0; j < inputtrees.size(); j++){//for each tree
			centroid_distances[i].insert(centroid_distances[i].begin() + j, tree_distances[i][j]);
		}
	}

	//While there is a change in a loop
	while(nclusters != oclusters){
		//Currently unused centroid history piece
	//	if(std::find(centroid_history.begin() + 1, centroid_history.end(), centroid_history[0]) != centroid_history.end()){
	//		cout << "breaking " << endl;
	//			break;
	//			}
	//centroid_history.resize(k);
	//	for(unsigned int i = k - 1; i > 0; i--){//for each history
	//			centroid_history[i] = centroid_history[i - 1];
	//		}
	//	centroid_history[1] = centroid_history[0];
		
		oclusters = nclusters;
		nclusters.clear();
	//	centroid_history[0].clear();
	

		//fill nclusters with the centroids
		for(std::set<unsigned int>::iterator pos = centroids.begin(); pos!= centroids.end(); ++pos){//for each centroid
			set <unsigned int> temp;
			nclusters.push_back(temp);
		}
		//Find the nearest neighbor in centroids of each tree and add that tree to the centroid in nclusters
		unsigned int count = 0;
		for(std::set<unsigned int>::iterator pos = inputtrees.begin(); pos != inputtrees.end(); ++pos){//for each tree
				unsigned int closest_centroid;
				closest_centroid = closest_neighbor(count,centroids,centroid_distances).first;
				//Add the tree to it's closest centroid
				nclusters[closest_centroid].insert(count);
				count++;
		}
		//Recompute the centroid distance matrix
		recompute_centroid_distances(centroid_distances, tree_distances, inputtrees.size(), nclusters, oclusters);
	}
	//Adds the cluster information to the return set
//	typedef std::map < unsigned int, set <unsigned int> >::iterator it_map;
//	for(it_map iterator = nclusters.begin(); iterator != nclusters.end(); iterator++){//for each cluster
//		set<unsigned int> tempset;
//		tempset.clear();
//		for(std::set<unsigned int>::iterator pos = iterator->second.begin(); pos != iterator->second.end(); ++pos){// for each tree in the cluster
//			set<unsigned int>::iterator it = inputtrees.begin();
//			std::advance(it, *pos);
//			tempset.insert(*it);
//		}
//		ret_clusters.push_back(tempset); //add to the return set
//	}

	ret_clusters = nclusters;

	//Prints out the clusters map (maybe less than useful on very large datasets)
//	typedef std::map< unsigned int, set < unsigned int >>::iterator it_map;
	for(unsigned int i = 0; i < ret_clusters.size(); i++){//for each cluster
		cout << "Cluster : " << i << " Trees : " ;
		for(std::set<unsigned int>::iterator  pos = ret_clusters[i].begin(); pos != ret_clusters[i].end(); ++pos){//for each tree in the cluster
			set<unsigned int>::iterator it = inputtrees.begin();
			std::advance(it, *pos);
			cout << *it << " " ;
		}
		cout << endl;
	}


	return ret_clusters;
}

//various tests that have been used for clusters
void TestClust(){

	set <unsigned int> testset;
	for(unsigned int i = 0; i < 11; i++){
		testset.insert(i);
	}
	silhouette(kmeans_clust(testset, 3, "rf"), "rf");
	silhouette(agglo_clust(testset, "rf"), "rf");

	set <unsigned int> testset2;
	for(unsigned int i = 3; i < 10; i++){
		testset2.insert(i);
	}
	kmeans_clust(testset2, 3, "rf");
	cout << endl;
	agglo_clust(testset2, "rf");

	//Tests silhouette
	//Uncommenting print statements (in silhouette) will show the data being used
	//to calculate the silhouette
	//
	//Based on manual calculations cluster silhouettes should be
	//Cluster 1: -.15
	//Cluster 2: .6333(repeating)
	//Cluster 3: -.2
	//set <unsigned int> testset3;
	//testset3.insert(1);
	//testset3.insert(4);
	//set <unsigned int> testset4;
	//testset4.insert(3);
	//testset4.insert(8);
	//set <unsigned int> testset5;
	//testset5.insert(10);
	//testset5.insert(6);
	//vector <set < unsigned int > > test_tree_vect;
	//test_tree_vect.push_back(testset3);
	//test_tree_vect.push_back(testset4);
	//test_tree_vect.push_back(testset5);
	//silhouette(test_tree_vect);
}


//Returns the vector of taxa in are present in all of the input trees
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
