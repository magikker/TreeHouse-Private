#ifndef _BIPARTITION_H
#define _BIPARTITION_H

//#include <map>
#include <vector>
//#include <stdexcept>
//#include <string>
#include <iostream>
//#include "global.h"


#include <boost/dynamic_bitset.hpp>

using namespace std;

class Bipartition {
public:
	Bipartition(int len);

	void set(int place, bool val);
	void add_tree(int tree, float branch_len);
	bool in_bitstring(int len);
	int number_of_ones();
	int number_of_zeros();
	bool is_zero(int position);
	bool same_bitstring_value(int position1, int position2);
	bool same_bitstring_value(vector<int> positions);
	bool is_one(int position);
	void print_bitstring(bool endline);
	void print_trees(bool endline);
	void print_branches(bool endline);
	void print_line();

private:
	boost::dynamic_bitset<> bitstring;
	vector<unsigned int> trees;
	vector<float> branchlengths;
};

#endif

