#include <vector>
#include <iostream>
#include <set>
#include <algorithm>
#include <string>

#include <stdlib.h>
#include <string.h>

#include "label-map.hh"
#include "BipartitionTable.h"
#include "THGlobals.h"
#include "global.h"
#include "UtilityFunctions.h"

#include "pqlsymbol.h"
#include "newick.h"

double average_ancestral_distance(unsigned int taxon1, unsigned int taxon2);
double average_ancestral_distance(unsigned int taxon1, unsigned int taxon2, set<unsigned int> treeSet);
int distance_between_taxa(unsigned int taxon1, unsigned int taxon2, unsigned int tree);
int distance_to_common_ancestor(unsigned int taxon1, unsigned int taxon2, unsigned int tree);
unsigned int distance_to_root(unsigned int taxon1, unsigned int tree);
double average_depth(unsigned int tree);
double average_distance_between_taxa(unsigned int taxon1, unsigned int taxon2);

void testDistance();
vector<unsigned int> checkForDuplicateBitstrings();
