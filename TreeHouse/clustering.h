#ifndef _CLUSTERING_H_
#define _CLUSTERING_H_

#include <vector>
#include <stack>
#include <iomanip> //Formatting output

//Global Vars like the hashtable.
#include "THGlobals.h"
#include "UtilityFunctions.h"

#include "global.h"
#include "distance.h"

//The label map class. 
#include "label-map.hh"

#include "newick.h"
#include "SearchFunctions.h"
#include "AnalysisFunctions.h"
using namespace std;
//computes the silhouette widths of your input clusters for cluster analysis
vector<float> silhouette(vector < set < unsigned int >> inputclusts, string dist_type);

//Computes the rand index of the two given clusterings for a set of trees
float rand_index(vector < set < unsigned int > > cluster1, vector < set < unsigned int > > cluster2);

//Forms clusters of the input trees based on a distance type using the agglomerative method
vector< set < unsigned int > > agglo_clust(set <unsigned int > inputtrees,unsigned int numclusts, string dist_type);

//Forms clusters of the input trees based on a distance type using the kmeans method
vector< set < unsigned int > > kmeans_clust(set <unsigned int> inputtrees, unsigned int k, string dist_type);

//Forms a cluster of the input trees based on a distance type using the dbscan method
vector < set < unsigned int> > dbscan_clust(set<unsigned int> treeset, unsigned int eps, unsigned int minpts, string dist_type);

//Tests for visualization
void mdsTests();

//Uses mds to display clusters in gnuplot
void display_clusters(string type,string dist_type, vector <set < unsigned int > >clusters);

//Displays a heatmap of the distances of the given treeset based on the given dist_type
void display_heatmap(set <unsigned int> treeset, string filename, string dist_type);

//Test for clustering functions
void TestClust();

#endif
