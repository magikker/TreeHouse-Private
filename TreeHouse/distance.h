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

#include "pqlsymbol.h"
#include "newick.h"

int distance_between_taxa(unsigned int taxon1, unsigned int taxon2, unsigned int tree);
int distance_to_common_ancestor(unsigned int taxon1, unsigned int taxon2, unsigned int tree);
unsigned int distance_to_root(unsigned int taxon1, unsigned int tree);
double averageDepth(unsigned int tree);
double average_distance_between_taxa(unsigned int taxon1, unsigned int taxon2);

void testDistance();
