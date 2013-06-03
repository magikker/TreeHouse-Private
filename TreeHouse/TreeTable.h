#ifndef _TREETABLE_H
#define _TREETABLE_H

#include <map>
#include <vector>
#include <stdexcept>
#include <string>

using namespace std;

class TreeTable {

public:
	
	//For each tree
	vector< vector<float> > branches; 			    //stores the branches (in terms of inverted index)		Don't need?
	vector< vector<bool* > > bipartitions; 		    // 															???
	vector< vector<unsigned int> > bs_sizes; 		//The sizes of all the bipartitions 


	//unsigned NUM_TREES;// = 0; // number of trees
	//unsigned NUM_TAXA;//  = 0; // number of taxa
	//vector < vector < unsigned int > > searchtable;
	//unsigned int * hash_lengths;
	std::map <int, vector<int> > dups;


	void print_dups(){
		cout << "dups: " << endl;
		for (unsigned int i = 0; i < dups.size(); i++){
			cout << i << ": ";
			for (unsigned int j = 0; j < dups[i].size(); j++){
				cout << dups[i][j] << " ";
			}
			cout << endl;
		}
	}
	
};

#endif
