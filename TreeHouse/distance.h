#include <vector>
#include <iostream>
#include <set>
#include <algorithm>
#include <string>
#include <numeric>
#include <cmath>

#include <stdlib.h>
#include <string.h>

#include "label-map.hh"
#include "BipartitionTable.h"
#include "THGlobals.h"
#include "global.h"
#include "UtilityFunctions.h"

#include "pqlsymbol.h"
#include "newick.h"
//#include "hungarian.h"


vector < vector < float>> compute_distances(set <unsigned int> treeset, string dist_type);

vector < vector <float>> compute_distances(vector <unsigned int> treevect, string dist_type);

 vector<vector<float>> taxaSimilarityMatrix(set<unsigned int> in);

float CRFDistance(int tree1, int tree2);
float BipartRFDistance(int tree1, int tree2);
float BipartRFStringDistance(string tree1, string tree2);
float CRFStringDistance(string tree1, string tree2);


pair<int, int> numLeftRight(string nw, unsigned int &total );
int sumPair(pair<int, int> x);
bool isBifurcating(string nw);
bool isBifurcating(string nw, int depth);

vector < vector < float > > distanceWrapper (set < unsigned int >, int);

unsigned int edit_distance_greedy(unsigned int tree1, unsigned int tree2);
double edit_distance_total(unsigned int tree1, unsigned int tree2);
unsigned int edit_distance_minimum(unsigned int tree1, unsigned int tree2);
double edit_distance_minimum_coverage(unsigned int tree1, unsigned int tree2);
double edit_distance_average(unsigned int tree1, unsigned int tree2);

pair<set<unsigned int>, set<unsigned int>> rfDistanceSet(int tree1, int tree2);
vector< set<int> > editDistanceMatrixSet(set<unsigned int> rf, set<unsigned int> rf2);
vector<vector<int>> editDistanceMatrixVector(set<unsigned int> rf, set<unsigned int> rf2);
float rfDistance(int tree1, int tree2);
void printRFset(int tree1, int tree2);

void testDistance();
vector<unsigned int> checkForDuplicateBitstrings();
void testDepthVariance();
void testAverageDepth();
void testNewick();
void testEdit();
void testCalculateC();
void testIsBifurcating();
