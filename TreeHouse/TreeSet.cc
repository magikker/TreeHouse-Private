//file TreeSet.cc
#include "TreeSet.h"
using namespace std;

TreeSet::TreeSet(set<unsigned int> treeids){
	tree_ids = treeids;
	is_inverse = false;
	no_branch_lengths = true;
}

TreeSet::TreeSet(set<unsigned int> treeids, bool inverseids){
	tree_ids = treeids;
	is_inverse = inverseids;
	no_branch_lengths = true;
}

TreeSet::TreeSet(set<unsigned int> treeids, bool inverseids, vector<float> branchlengths){
	tree_ids = treeids;
	is_inverse = inverseids;
	branch_lengths = branchlengths;
	no_branch_lengths = false;
}



//set<unsigned int> TreeSet::get_trees(){
//	return tree_ids;
//}

