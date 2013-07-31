#include "clustering.h"

using namespace std;
//----------------------Cluster Analysis------------------------------------
//Returns the average distance for each tree to the cluster containing it (silhouette helper)
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



//Returns the average distance for each tree to the cluster with the smallest average distance from the tree (silhouette helper)
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
	distances = compute_distances(treeset, dist_type);

	//Computes diff_within_cluster
	diff_within_cluster = diff_to_own_clusters(inputclusts, distances);

	//Computes the neighboring_clusters
	neighboring_cluster = neighboring_cluster_diffs(inputclusts, distances);

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

//Returns a set of pairs representing the trees which share a cluster in the input cluster set
set <pair <unsigned int, unsigned int> > get_cluster_pairs (vector < set < unsigned int> > clusters){
	set <pair <unsigned int, unsigned int> > retset;
	
	//Inserts the pairings of trees which share a cluster into the return set
	for(unsigned int i = 0; i < clusters.size(); i++){//for each cluster
		for(std::set<unsigned int>::iterator pos = clusters[i].begin(); pos != clusters[i].end(); ++pos){//for each tree in cluster
			for(std::set<unsigned int>::iterator pos2 = clusters[i].begin(); pos2 != clusters[i].end(); ++pos2){//for each other tree in cluster
				if(*pos != *pos2){
				retset.insert(make_pair(*pos, *pos2));
				}
			}
		}
	}
	return retset;
}

//Returns a set of pairs representing the trees which do not share a cluster in the input cluster set
//(will contain both a pair and its mirror)
set <pair <unsigned int, unsigned int> > get_unclustered_pairs (vector < set < unsigned int> > clusters){
	set < pair < unsigned int, unsigned int > > retset;

	//Inserts the pairings of tree which do no share a cluster into the return set
	for (unsigned int i = 0; i < clusters.size(); i++){//for each cluster
		for(unsigned int j = 0; j < clusters.size(); j++){//for each other cluster
			if(i != j){
			for(std::set<unsigned int>::iterator pos = clusters[i].begin(); 
					pos != clusters[i].end(); ++pos){//for each tree in the cluster
				for(std::set<unsigned int>::iterator pos2 = clusters[j].begin(); 
						pos2 != clusters[j].end(); ++pos2){//for each tree in the other cluster
					retset.insert(make_pair(*pos, *pos2));
					retset.insert(make_pair(*pos2, *pos));
				}
			}
			}
		}
	}
	return retset;
}

//Returns a number between 0 and 1 representing how similar 2 clusterings are, where 1 is identical and 0 is completely different
float rand_index (vector < set < unsigned int > > clusters1, vector < set < unsigned int > > clusters2){
	float retvalue;
	//# of Pairs of trees that are in the same cluster in clusters1, and the same cluster in clusters2
	unsigned int a;
	//# of Pairs of trees that are in different clusters in clusters1, and different clusters in cluststers2
	unsigned int b;
	//# of Pairs of trees that are clustered together in clusters1, and clustered apart in clusters2
	unsigned int c;
	//# of Pairs of trees that are clustered apart in clusters1, and clustered together in clusters2
	unsigned int d;
	//Paired trees in clusters1
	set <pair <unsigned int, unsigned int> > same1;
	//Paired trees in clusters2
	set <pair <unsigned int, unsigned int> > same2;
	//Separated trees in clusters1
	set <pair <unsigned int, unsigned int> > separate1;
	//Separated trees in clusters2
	set <pair <unsigned int, unsigned int> > separate2;

	//Fill the sets
	same1 = get_cluster_pairs(clusters1);
	same2 = get_cluster_pairs(clusters2);
	separate1 = get_unclustered_pairs(clusters1);
	separate2 = get_unclustered_pairs(clusters2);
	
	//Obtain the value to be used for the calculation (a, b, c, d)
	set<pair <unsigned int, unsigned int> > tempset;
	set_intersection(same1.begin(), same1.end(), same2.begin(), same2.end(),std::inserter(tempset, tempset.end()));
	a = tempset.size();
	tempset.clear();

	set_intersection(separate1.begin(), separate1.end(), separate2.begin(), separate2.end(),std::inserter(tempset, tempset.end()));
	b = tempset.size();
	tempset.clear();

	set_intersection(same1.begin(), same1.end(), separate2.begin(), separate2.end(),std::inserter(tempset, tempset.end()));
	c = tempset.size();
	tempset.clear();

	set_intersection(separate1.begin(), separate1.end(), same2.begin(), same2.end(),std::inserter(tempset, tempset.end()));
	d = tempset.size();
	tempset.clear();

	//Compute the rand index for the clusterings
	float temp1 = a + b;
	float temp2 = a + b + c + d;

	retvalue = temp1 / temp2;

	cout << "Rand Index is: " << retvalue << endl;
	return retvalue;
}

//Computes n choose k for the given n and k
long long nchoosek(unsigned long long n, unsigned long long k) {
	if ( k > n ) {
		return 0;
	}
	unsigned long long r = 1;
	for(unsigned int d = 1; d <= k; d++){
		r *= n--;
		r /= d;
	}
	return r;
}


float adjusted_rand_index(vector < set < unsigned int > > clusters1, vector < set < unsigned int > > clusters2){
	float retvalue;
	//Holds the number of elements in common between the clusters in each clustering
	vector < vector < unsigned int> > contingency_table;
	contingency_table.resize(clusters1.size(), vector <unsigned int> (clusters2.size(), 0));
	//Holds the sum values for computation
	vector <float> asums;
	vector <float> bsums;
	unsigned int num_trees = 0;

	for(unsigned int i = 0; i < clusters1.size(); i++){//for each cluster in clusters1
		num_trees += clusters1[i].size();
	}

	//Fills the contingency table
	for(unsigned int i = 0; i < clusters1.size(); i++){//for each cluster in clusters1
		float tempsum = 0;
		for(unsigned int j = 0; j < clusters2.size(); j++){//for each cluster in clusters2
			set < unsigned int > tempset;
			tempset.clear();
			set_intersection(clusters1[i].begin(), clusters1[i].end(), clusters2[j].begin(), clusters2[j].end(), std::inserter(tempset, tempset.end()));
			contingency_table[i][j] = tempset.size();
			contingency_table[j][i] = tempset.size();
			tempsum += tempset.size();
		}
		asums.push_back(tempsum);
	}

	//Fills the bsums vector
	for(unsigned int i = 0; i < clusters2.size(); i++){//for each cluster in clusters2
		float tempsum = 0;
		for(unsigned int j = 0; j < clusters1.size(); j++){//for each cluster in clusters1
			tempsum += contingency_table[j][i];
		}	
		bsums.push_back(tempsum);
	}

	//Computes the necessary sum from the contingency_table
	long long nsum = 0;
	for(unsigned int i = 0; i < contingency_table.size(); i++){//for each from clusters1
		for(unsigned int j = 0; j < contingency_table[i].size(); j++){//for each from clusters2
			nsum += nchoosek(contingency_table[i][j], 2);
		}
	}
/*
	for(unsigned int i = 0; i < contingency_table.size(); i++){
		cout << "Row: " << i << ": ";
		for (unsigned int k = 0; k < i; k++){
			cout << "  ";
		}
		for (unsigned int j = i; j < contingency_table.size(); j++){
			cout << contingency_table[i][j] << " ";
		}
		cout << endl;
	}
*/
	//Computes the sums totals
	long long asum = 0;
	long long bsum = 0;
	for(unsigned int i = 0; i < asums.size(); i++){//for each asum
		asum += nchoosek(asums[i], 2);
	}
	for(unsigned int i = 0; i < bsums.size(); i++){//for each bsum
		bsum += nchoosek(bsums[i], 2);
	}
	
	//Computes the nchoose 2 of the number of trees
	long long n = nchoosek(num_trees, 2);


	retvalue = ((nsum - (asum * bsum / n)) / ((asum + bsum) / 2.0 - (asum * bsum / n))) * 1.0;

//	cout << "Nsum: " << nsum << " asum : " << asum << " bsum : " << bsum << " n : " << n << endl;

	cout << "Adjusted Rand Index: " << retvalue << endl;

	return retvalue;
}




//--------------------------Forming Clusters----------------------------------

//Returns a pair representing the closest neighbor of a cluster, and the distances between the two clusters (agglo_clust helper)
pair<unsigned int, float> closest_neighbor(unsigned int clustid, set <unsigned int> checkset, vector <vector <unsigned int> > distances){
	//neighbor to clustid to be returned
	pair<unsigned int, float> retneighbor;

	float bestdist = 100;
	unsigned int best_neighbor = 0;
	for(std::set<unsigned int>::iterator pos = checkset.begin(); pos != checkset.end(); pos++){//for each cluster in checkset
		float tempdist;
		//If the distance from clustid to anything in the checkset is
		//better than the best distance
		tempdist = distances[*pos][clustid];
		if (tempdist < bestdist){
			bestdist = tempdist;
			//Best neighbor equal current item from checklist
			best_neighbor = *pos;
		}
	}
	retneighbor = make_pair(best_neighbor, bestdist);

	return retneighbor;
}

//Returns a map of clusters with the two input clusters merged into a single cluster (agglo_clust helper)
void  merge(unsigned int clust1, unsigned int clust2, 
		map <unsigned int, set <unsigned int > > &clusters){
	//Take the elements from clust2's set and insert them into clust1
	clusters[clust1].insert(clusters[clust2].begin(),clusters[clust2].end());
	//Remove clust2
	clusters.erase(clust2);
}

//Returns a distance matrix where the distances are recomputed for merged clusters (agglo_clust helper)
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
vector < set < unsigned int > > agglo_clust (set <unsigned int> inputtrees, unsigned int numclusts, string dist_type){
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
	unsigned int fin_num = numclusts;

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
	distances = compute_distances(inputtrees, dist_type);


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

/*
	//Prints out the clusters map
	typedef std::map< unsigned int, set < unsigned int >>::iterator it_map;
	for(it_map iterator = clusters.begin(); iterator != clusters.end(); iterator++){
		cout << "Cluster : " << iterator->first << " Trees : " ;
		for(std::set<unsigned int>::iterator  pos = iterator->second.begin(); pos != iterator->second.end(); ++pos){
			cout << *pos << " " ;
		}
		cout << endl;
	}
*/


	return ret_clusters;
}

//Returns the distances for each centroid to each tree (probably not the most efficient method) (kmeans_clust helper)
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
//	tree_distances = compute_distancesv(inputtrees, dist_type);

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
	tree_distances = compute_distances(shuffled, dist_type);
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
	for(unsigned int i = 0; i < nclusters.size(); i++){//for each cluster
		set<unsigned int> tempset;
		tempset.clear();
		for(std::set<unsigned int>::iterator pos = nclusters[i].begin(); pos != nclusters[i].end(); ++pos){//for each tree in cluster
		       set<unsigned int>::iterator it = inputtrees.begin();
		       std::advance(it, *pos);
		       tempset.insert(*it);
		}
		ret_clusters.push_back(tempset); //add to the return set
	}

	//Prints out the clusters (maybe less than useful on very large datasets)
	/*
	for(unsigned int i = 0; i < ret_clusters.size(); i++){//for each cluster
		cout << "Cluster : " << i << " Trees : " ;
		for(std::set<unsigned int>::iterator  pos = ret_clusters[i].begin(); pos != ret_clusters[i].end(); ++pos){//for each tree in the cluster
			cout << *pos << " " ;
		}
		cout << endl;
	}
	*/

	return ret_clusters;
}

//Returns a mapping of the input trees to all trees which are no more than eps distance away / dissimilar from them
//(dbscan_clust helper)
map <unsigned int, vector <unsigned int> > compute_eps_neighborhoods(set <unsigned int> treeset, unsigned int eps, string dist_type){
	
	map <unsigned int, vector < unsigned int> > retmap;
	vector < vector < unsigned int> > distances;
	//Computes the distance matrix to be used
	distances = compute_distances(treeset, dist_type);

	unsigned int i = 0; //Iterates over the trees while also storing the value for the matrix
	for(std::set<unsigned int>::iterator pos = treeset.begin(); pos != treeset.end(); ++pos){//for each tree
		
		vector <unsigned int> temp;
		temp.clear();
		retmap.insert(make_pair(*pos, temp));
		
		unsigned int j = 0; //Iterates over the trees while also storing the value for the matrix
		for(std::set<unsigned int>::iterator pos1 = treeset.begin(); pos1 != treeset.end(); ++pos1){//for each other tree
			if (i != j){
				//if the distance is small enough, add the mapping from the current tree to the inner loop tree
				if( distances[i][j] <= eps){
					retmap[*pos].push_back(*pos1);
				}
			
			}

			j++;
		}

		i++;
	}
		
	return retmap;
}

//Returns a cluster starting from the input tree and including all members of the input
//tree's eps neighborhood who are not already in a cluster as well as the members of subsequent
//eps neighborhoods when those trees belong to dense neighborhoods (dbscan_clust helper)
set <unsigned int> grow_cluster(unsigned int inputtree, map <unsigned int, vector < unsigned int> > eps_neighborhoods, set<unsigned int> &visited,
		vector < unsigned int > temphood, unsigned int eps, unsigned int minpts,set<unsigned int> &clustered){
	
	set <unsigned int> retset;

	//insert the current tree into the growing cluster and into the clustered set
	retset.insert(inputtree);
	clustered.insert(inputtree);

	for (unsigned int j = 0; j < temphood.size(); j++){//for each tree 
		//set to hold eps_neighborhood for each tree
		vector <unsigned int> tempvect;
		tempvect.clear();
		//If tree is not visited mark it as such
		if(visited.find(temphood[j]) == visited.end()){
			visited.insert(temphood[j]);
			tempvect = eps_neighborhoods[temphood[j]];
			//If the tempset is of appropriate size
			if(tempvect.size() >= minpts){
				//Merge the two while avoiding duplicates
				for(unsigned int i = 0; i < tempvect.size(); i++){//for each tree in tempset
					if(std::find(temphood.begin(), temphood.end(), tempvect[i]) != temphood.end()){
						temphood.push_back(tempvect[i]);
					}
				}
			}
			//If the tree is not already clustered, add it to the current cluster
			if(clustered.find(temphood[j]) == clustered.end()){
				retset.insert(temphood[j]);
				clustered.insert(temphood[j]);
			}
		}
	}
	return retset;
}


//A clustering algorithm that takes in a set of trees, a given distance eps which is the minimum distance
//two trees can be from each other and still be considered, and minpts, the minimum number of pts required for a cluster
vector < set < unsigned int > >dbscan_clust(set <unsigned int> treeset, unsigned int eps, unsigned int minpts, string dist_type){
	vector < set < unsigned int> > retset;
	map <unsigned int,  vector < unsigned int > > eps_neighborhoods;
	set <unsigned int> visited;
	set <unsigned int> noise;
	set <unsigned int> clustered;


	//Computes the eps neighborhoods all at once and stores them for use 
	eps_neighborhoods = compute_eps_neighborhoods(treeset, eps, dist_type); //Currently unwritten


	for(std::set<unsigned int>::iterator pos = treeset.begin(); pos != treeset.end(); ++pos){//for each tree
		//If the tree has not been visited, add it to the visited set and acquire it's eps_neighborhood
		if(visited.find(*pos) == visited.end()){
			visited.insert(*pos);
			vector <unsigned int> temphood = eps_neighborhoods[*pos];
			//If the neighborhood is not large enough, consider the point noise
			if(temphood.size() < minpts){
				noise.insert(*pos);
			}
			//Otherwise grow a cluster starting from this tree and place it in the return set
			else{
				set<unsigned int> tempclust;
				tempclust.clear();
				tempclust = grow_cluster(*pos, eps_neighborhoods, visited, temphood, eps, minpts, clustered);
				retset.push_back(tempclust);
			}
		}
	}	
	//Adds the noise cluster to the return set
	retset.push_back(noise);

	//Prints out the clusters (maybe less than useful on very large datasets)
	for(unsigned int i = 0; i < retset.size(); i++){//for each cluster
		cout << "Cluster : " << i << " Trees : " ;
		for(std::set<unsigned int>::iterator  pos = retset[i].begin(); pos != retset[i].end(); ++pos){//for each tree in the cluster
			cout << *pos << " " ;
		}
		cout << endl;
	}
	
	return retset;
}

//Recomputes the distance for burnin_clust (burnin_clust helper)
void recompute_distances(vector < vector < unsigned int > > distances, 
				vector < pair < float, float > > &rem_distances, vector < pair < float, float> > &burned_distances,
				unsigned int burnedSize, unsigned int remSize, unsigned int treenum, unsigned int setsize){

	for(unsigned int i = treenum; i < setsize; i++){//for each tree after the treenum
	
		//Recompute the remaining distances for each tree
		rem_distances[i].second -= distances[i][treenum];
		rem_distances[i].first = rem_distances[i].second / remSize;
		//Recompute the burned distances for each tree
		burned_distances[i].second += distances[i][treenum];
		burned_distances[i].first = burned_distances[i].second / burnedSize;
	}
}


//A function for determining the trees to throw out after a run of a program such as Mr. Bayes 
//Clustering approach being attempted as alternative method
//(in the works, currently theoretical and largly untested) 
vector < set < unsigned int > > burnin_clust (set <unsigned int> treeset, string dist_type){
	vector < set < unsigned int> > retset;
	//Stores and computes the distance matrix
	vector < vector < unsigned int > > distances;
	distances = compute_distances(treeset, dist_type);
	set <unsigned int > burnedset;
	set <unsigned int > keepset = treeset;

	//Stores the average distance from each tree to the remaining trees
	//and the total distance before averaging
	vector < pair < float, float > > rem_distances;

	//Stores the average distance from each tree to the already clustered trees
	//and the toal distance before averaging
	vector < pair < float, float > > burned_distances;

	//Fills the rem_distances vector (seems to make function very slow - should be revised)
	//Possibly run on subset of the input dataset
	for(unsigned int i = 0; i < treeset.size(); i++){//For each tree
		float tempsum = 0;
		for(unsigned int j = 0; j < treeset.size(); j++){//For each other tree
			if(i != j){
				tempsum += distances[i][j];
			}
		}

		rem_distances.push_back(make_pair(tempsum / treeset.size(), tempsum));
	}
	cout << "finished rem_distances" << endl;

	//Fills the past_distances vector to begin
	for (unsigned int i = 0; i < treeset.size(); i++){//for each tree
		burned_distances.push_back(make_pair(0,0));
	}


	unsigned int treenum = 0;

	for(std::set<unsigned int>::iterator pos = treeset.begin(); pos != treeset.end(); ++pos){//for each tree
		//If closer to the remaining trees than to the burned set
		//we're done
		cout << "Distances : " << rem_distances[treenum].first << " " << burned_distances[treenum].first << endl;
		if(rem_distances[treenum].first < burned_distances[treenum].first){	
			break;
		}
		//Else, burn the current tree and recompute_distances
		else{
			burnedset.insert(*pos);
			keepset.erase(*pos);
			recompute_distances(distances, rem_distances, burned_distances, burnedset.size(), keepset.size(), treenum, treeset.size());
		}
		treenum++;
	}


	retset.push_back(burnedset);
	retset.push_back(keepset);

	//Prints out the clusters (maybe less than useful on very large datasets)
	
	for(unsigned int i = 0; i < retset.size(); i++){//for each cluster
		cout << "Cluster : " << i << " Trees : " ;
		for(std::set<unsigned int>::iterator  pos = retset[i].begin(); pos != retset[i].end(); ++pos){//for each tree in the cluster
			cout << *pos << " " ;
		}
		cout << endl;
	}
	

	return retset;

}



//various tests that have been used for clusters
void TestClust(){

	set <unsigned int> testset;
	//vector < set <unsigned int> > testvect;
	//set <unsigned int> piece;
	//set <unsigned int> piece2;

	for(unsigned int i = 0; i < 10000; i++){
		testset.insert(i);
	}

//	burnin_clust(testset, "rf");

/*	
	piece.insert(5);
	piece.insert(0);
	piece.insert(2);
	piece.insert(8);
	piece2.insert(3);
	piece2.insert(4);
	piece2.insert(6);
	piece2.insert(7);
	piece2.insert(9);
	piece2.insert(10);
	testvect.push_back(piece);
	testvect.push_back(piece2);
*/


	//kmeans_clust(testset, 2, "rf");
//	cout<< "adjusted_rand_index" << endl;
//	adjusted_rand_index(kmeans_clust(testset, 2, "rf"),kmeans_clust(testset,2,"rf"));
//	adjusted_rand_index(kmeans_clust(testset, 3, "rf"),kmeans_clust(testset,2,"rf"));
//	adjusted_rand_index(kmeans_clust(testset, 2, "rf"),testvect);
//	cout << "rand_index" << endl;
//	rand_index(kmeans_clust(testset, 2, "rf"),kmeans_clust(testset,2,"rf"));
//	rand_index(kmeans_clust(testset, 3, "rf"),kmeans_clust(testset,2,"rf"));
	//cout << "First test" << endl;
	//dbscan_clust(testset, 1, 2, "rf");
	//cout << "All trees singular test" << endl;
	//dbscan_clust(testset, 0, 0, "rf");
	//cout << "All trees noise test" << endl;
	//dbscan_clust(testset, 0, 2, "rf");


	//silhouette(kmeans_clust(testset, 3, "rf"), "rf");
	//silhouette(agglo_clust(testset, "rf"), "rf");

	//set <unsigned int> testset2;
	//for(unsigned int i = 3; i < 10; i++){
	//	testset2.insert(i);
	//}
	//kmeans_clust(testset2, 3, "rf");
	//cout << endl;
	//agglo_clust(testset2, "rf");

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


//-------------------Visualizing Trees and Clusters----------------------------

//Writes a distance matrix to a file to be read by another program
void write_dmatrix(vector <vector < unsigned int > > dmatrix, string filename, bool overwrite){
	ofstream cfile;
	//Checks to see if necessary to overwrite any existing file. Only doesn't overwrite in the case 
	//of write_dmatrix_dred, which has already overwritten any existing file with a single line we
	//want to keep
	if(overwrite){
		cfile.open("./temp/" + filename, ios::out | ios::trunc);
	}
	else{
		cfile.open("./temp/" + filename, ios::out | ios::app);
	}
	for (unsigned int i = 0; i < dmatrix.size(); i++){//for each tree
		for(unsigned int j = 0; j < dmatrix.size(); j++){//for each other tree
			if(i == j){
			       cfile << 0 << " ";
			}	
			else{
				cfile << dmatrix[i][j] << " ";
			}
		}
		cfile << endl;
	}
	cfile.close();
	cout << "created " << filename << endl;
}

//Wrapper function for write_dmatrix that adds an extra line used by dredviz
void write_dmatrix_dred(vector <vector <unsigned int > > dmatrix, string filename){
	ofstream cfile;
	cfile.open("./temp/" + filename, ios::out | ios::trunc);
	cfile << dmatrix.size() <<endl; //Adds necessary value to top of file for dredviz
	cfile.close();
	write_dmatrix(dmatrix, filename, false);
}

//Displays a heatmap for the given set of trees (and stores it to a filename) based on a distance type
void display_heatmap(set <unsigned int> treeset, string filename, string dist_type){

	//Stores the distance matrix and prints it to a file
	vector <vector < unsigned int> > distances = compute_distances(treeset, dist_type);
	string dmatrixfile;
	dmatrixfile = filename + "_matrix";
	write_dmatrix(distances, dmatrixfile, true);

	//Values necessary in the heatmap
	unsigned int maxdist = 0;
	float range = treeset.size() - .5;

	//Obtain the maximum distance in the distance matrix
	for(unsigned int i = 0; i < distances.size(); i++){//for each tree
		for(unsigned int j = 0; j < distances[i].size(); j++){//for each other tree
			if(distances[i][j] > maxdist){
				maxdist = distances[i][j];
			}
		}
	}

	//Opens a stream to write to a file
	ofstream cfile;
	cfile.open("./temp/dis_heatmap.p", ios::out | ios::trunc);
//Commands for gnuplot script to run
	//Sets file output type and name
	cfile << "# set terminal png transparent nocrop enhanced font arial 8 size 420, 320" << endl;	
	cfile << "set output '" + filename + ".png'" << endl;
	//Establishes basic heatmp format
	cfile << "set bar 1.000000" << endl;
	cfile << "set style rectangle back fc lt -3 fillstyle solid 1.00 border -1" << endl;
	cfile << "unset key" << endl;
	cfile << "set view map" << endl;
	cfile << "set xtics border in scale 0,0 mirror norotate offset character 0, 0, 0" << endl;
	cfile << "set ytics border in scale 0,0 mirror norotate offset character 0, 0, 0" << endl;
	cfile << "set ztics border in scale 0,0 mirror norotate offset character 0, 0, 0" << endl;
	cfile << "set nocbtics" << endl;
	//Sets a title to be used
	cfile << "set title \"" + filename + "\"" << endl;
	//Sets the ranges and labels for display
	cfile << "set rrange [ * : * ] noreverse nowriteback # (currently [20.0000:0.0000])" << endl;
	cfile << "set trange [ * : * ] noreverse nowriteback # (currently [-5.0000:0.0000])" << endl;
	cfile << "set urange [ * : * ] noreverse nowriteback # (currently [-5.0000:0.0000])" << endl;
	cfile << "set vrange [ * : * ] noreverse nowriteback # (currently [-5.0000:0.0000])" << endl;
	//xrange and later yrange are based on the number of trees
	cfile << "set xrange [ -0.500000 : " << range << " ] noreverse nowriteback" << endl;
	cfile << "set ylabel offset character 0, 0 , 0 font \"\" textcolor lt -1 rotate by 90" << endl;
	cfile << "set y2label offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90" << endl;
	cfile << "set yrange [ -0.500000 : " << range << " ] noreverse nowriteback" << endl;
	cfile << "set cblabel \"Distance\"" << endl;
	cfile << "set cblabel offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90" << endl;
	//The maximum distance will be one end of the spectrum in the heatmap
	cfile << "set cbrange [ 0.00000 : " << maxdist << " ] noreverse nowriteback" << endl;
	cfile << "set locale \"C\"" << endl;
	//Determines the color palette, this is a fairly good standard for a heatmap I think
	cfile << "set palette rgb 33, 13, 10" << endl;
	cfile << "splot './temp/" + dmatrixfile +  "' matrix with image" << endl;
	//Pauses to allow the data to be displayed
	cfile << "pause -1" << endl;
	cfile.close();
	//String runs gnuplot on the script file and displays
	string cmdstring = "gnuplot \"./temp/dis_heatmap.p\"";
	cout << cmdstring << endl;
	system(cmdstring.c_str());
}
	
//Computes the multidimensional scaling of the input trees based on the given scaling type (nerv or lmds)
//takes in a vector for its treeset because order matters
void compute_mds(vector <unsigned int> treevect, string type, string filename, string dist_type){
	//computes the distances, function should probably accept a distance type input
	vector < vector < unsigned int > > distances = compute_distances(treevect, dist_type);
	string dmatrixfile;
	//Creates a filename for the matrix file that must be written
	dmatrixfile = filename + "_matrix";
	write_dmatrix_dred(distances, dmatrixfile);
	//Runs dredviz to compute the multidimensional scaling
	string cmdstring = "../dredviz/dredviz-1.0.2/" + type + " --inputdist " + "./temp/" + dmatrixfile + " --outputfile " + "./temp/" + filename;
	cout << cmdstring << endl;
	system(cmdstring.c_str());
}

//Displays a multidimensional scaling for a given input vector of trees
void display_mds(vector <unsigned int> treevect, string type, string dist_type){
	//Computes the mds
	compute_mds(treevect, type, "display_mds", dist_type);
	//Writes a gnuplot script file for the data
	ofstream cfile;
	cfile.open("./temp/graph.p", ios::out | ios::trunc);
	cfile << "set autoscale" << endl;
	cfile << "set xtic auto" << endl;
	cfile << "set ytic auto" << endl;
	cfile << "set title \" MDS of Input Trees \"" << endl;
	cfile << "plot \"temp/display_mds\" every ::1 using 1:2" << endl;
	cfile << "pause -1" << endl;
	//Runs gnuplot and the script file to display the mds
	string cmdstring = "gnuplot \"./temp/graph.p\"";
	cout << cmdstring << endl;
	system(cmdstring.c_str());
}

//Helper function to create a vector of trees in the order of their clustering
vector<unsigned int> trees_by_cluster(vector <set <unsigned int> > clusters){
	//Return vector
	vector <unsigned int> clustered_trees;

	for(unsigned int i = 0; i < clusters.size(); i++){//for each cluster
		for(std::set<unsigned int>::iterator pos = clusters[i].begin(); pos != clusters[i].end(); ++pos){//for each tree in cluster
			clustered_trees.push_back(*pos);
		}
	}
	return clustered_trees;
}


//Displays a clustering of trees using MDS and gnuplot
void display_clusters(string type, string dist_type, vector < set < unsigned int> > clusters){
	vector<unsigned int> cluster_trees;
	//Adds the trees to a vector based on clustering
	cluster_trees = trees_by_cluster(clusters);
	//Computes MDS to a file
	compute_mds(cluster_trees, type, "clustersDisplay", "rf");
	
	ofstream cfile;
	//Writes a gnuplot script file for the data
	cfile.open("./temp/clustergraph.p", ios::out | ios::trunc);
	cfile << "set autoscale" << endl;
	cfile << "set xtic auto" << endl;
	cfile << "set ytic auto" << endl;
	cfile << "set title \" MDS of Input Clusters \"" << endl;
	unsigned int offset = 1;
	for(unsigned int i = 0; i < clusters.size(); i++){//For each cluster
		if(i == 0){
			//Initial plot command
			cfile << "plot \"temp/clustersDisplay\" every ::0::";
			//Ending line
			cfile << clusters[i].size();
			//Which columns to use
			cfile << " using 1:2 ";
			//increment offset
			offset += clusters[i].size();
			//Title the data piece
			cfile << "title 'Cluster : ";
			cfile << i;
			cfile << "'";

		}
		else if(i == clusters.size() -1){
			cfile << " , ";
			cfile << " \"temp/clustersDisplay\"";
			cfile << " every ::";
			//Starting line
			cfile << offset;
			cfile << "::";
			cfile << offset + clusters[i].size();
			cfile << " using 1:2";
			cfile << " title 'Cluster : ";
			cfile << i;
			cfile << "'";
		}
		else{
			cfile << " , ";
			cfile << " \"temp/clustersDisplay\"";
			cfile << " every ::";
			cfile << offset;
			cfile << "::";
			cfile << offset + clusters[i].size();
			cfile << " using 1:2 ";
			offset += clusters[i].size();
			cfile << "title 'Cluster : ";
			cfile << i;
			cfile << "'";

		}
	}
	//Keeps gnuplot from closing
	cfile << endl;
	cfile << "pause -1" << endl;
	cfile.close();
	//Runs gnuplot and displays the data
	string cmdstring = "gnuplot \"./temp/clustergraph.p\"";
	system(cmdstring.c_str());
}

//Various tests used in setting up cluster display and MDS
void mdsTests(){
	vector< vector < unsigned int> > distances;
	vector <unsigned int> treevect;
	for (unsigned int i = 0; i < 11; i++){
		treevect.push_back(i);
	}

	set <unsigned int> treeset;
	for (unsigned int i = 0; i < 11; i++){
		treeset.insert(i);
	}

//	display_heatmap(treeset,"test", "rf");
//	vector <unsigned int> vect1;
//	vector<unsigned int> vect2;

//	vect1.push_back(0);
//	vect1.push_back(1);
//	vect1.push_back(2);
//	for(unsigned int i = 0; i < vect1.size(); i++){
//		cout << vect1[i] <<", ";
//	}
//	cout << endl;

//	write_dmatrix(compute_bipart_distances(treevect, "rf"), "heatmaptest");

//	vect2.push_back(2);
//	vect2.push_back(0);
//	vect2.push_back(1);

//	for(unsigned int i = 0; i < vect2.size(); i++){
//		cout << vect2[i] <<", ";
//	}
//	cout << endl;


//	write_dmatrix(compute_bipart_distances(vect2, "rf"), "vect2");
//	write_dmatrix(compute_bipart_distances(treevect, "rf"), "demo");
//	vector < set < unsigned int> > clusters = kmeans_clust(treeset, 2, "rf");
//	cout << "clusters defined" << endl;
//	vector < set <unsigned int> > clusters = kmeans_clust(treeset, 2, "eu");
//	vector <set <unsigned int> > cluster = kmeans_clust(treeset, 2, "rf");
//	distances = compute_bipart_distances(treeset, "rf");
//	write_dmatrix(distances, "testwrite");
//	compute_mds(treeset, "nerv", "test2");
//	display_mds(treeset, "lmds");
//	silhouette(clusters,"rf");
//	cout << "silhouette run" << endl;
//	display_clusters("lmds", clusters);
//	cout << "Display clusters fin" << endl;	
}
