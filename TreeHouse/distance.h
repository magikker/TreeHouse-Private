#include <vector>
#include <iostream>
#include <set>
#include <algorithm>
#include <string>
#include <numeric>

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

double average_ancestral_distance(unsigned int taxon1, unsigned int taxon2);
double average_ancestral_distance(unsigned int taxon1, unsigned int taxon2, set<unsigned int> treeSet);
int distance_between_taxa(unsigned int taxon1, unsigned int taxon2, unsigned int tree);
int distance_to_common_ancestor(unsigned int taxon1, unsigned int taxon2, unsigned int tree);
unsigned int distance_to_root(unsigned int taxon1, unsigned int tree);
double average_depth(unsigned int tree);
double expected_average_depth(unsigned int n);
double depth_variance(unsigned int tree);

double average_distance_between_taxa(unsigned int taxon1, unsigned int taxon2);
unsigned int hamming_distance_greedy(unsigned int tree1, unsigned int tree2);
double hamming_distance_total(unsigned int tree1, unsigned int tree2);
unsigned int hamming_distance_minimum(unsigned int tree1, unsigned int tree2);
double hamming_distance_average(unsigned int tree1, unsigned int tree2);

pair<set<unsigned int>, set<unsigned int>> rfDistanceSet(int tree1, int tree2);

unsigned int rfDistance(int tree1, int tree2);
void printRFset(int tree1, int tree2);

void testDistance();
vector<unsigned int> checkForDuplicateBitstrings();
void testDepthVariance();
void testAverageDepth();
void testNewick();
void testHamming();
