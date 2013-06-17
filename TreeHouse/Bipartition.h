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
	//constructor
	Bipartition(int len);
	Bipartition(boost::dynamic_bitset<> bs );
	Bipartition(boost::dynamic_bitset<> bs, vector<unsigned int> t, vector<float> branchl );
	Bipartition(boost::dynamic_bitset<> bs, vector<unsigned int> t, float branchl );
	
	//iterators
	std::vector<unsigned int>::iterator trees_begin();
	std::vector<unsigned int>::iterator trees_end();
	std::vector<unsigned int>::reverse_iterator trees_rbegin();
	std::vector<unsigned int>::reverse_iterator trees_rend();
	std::vector<unsigned int>::const_iterator trees_cbegin() const;
	std::vector<unsigned int>::const_iterator trees_cend() const;
	std::vector<unsigned int>::const_reverse_iterator trees_crbegin() const;
	std::vector<unsigned int>::const_reverse_iterator trees_crend() const;

	std::vector<float>::iterator branchlengths_begin();
	std::vector<float>::iterator branchlengths_end();
	std::vector<float>::reverse_iterator branchlengths_rbegin();
	std::vector<float>::reverse_iterator branchlengths_rend();
	std::vector<float>::const_iterator branchlengths_cbegin() const;
	std::vector<float>::const_iterator branchlengths_cend() const;
	std::vector<float>::const_reverse_iterator branchlengths_crbegin() const;
	std::vector<float>::const_reverse_iterator branchlengths_crend() const;

	//access
	bool get_bit(int bitindex);
	void set(int place, bool val);
	void add_tree(int tree, float branch_len);
	vector<unsigned int> get_trees();
	unsigned int get_tree(int index);
	boost::dynamic_bitset<> get_bitstring();
	float get_ave_branchlength();
	vector<float> get_branchlengths();
	vector<unsigned int> get_ones_indices();
	bool in_bitstring(int len);
	int number_of_ones();
	int number_of_zeros();
	bool is_zero(int position);
	bool same_bitstring_value(int position1, int position2);
	bool same_bitstring_value(vector<int> positions);
	bool is_one(int position);
	//size
	unsigned int bitstring_size();
	unsigned int trees_size();
	unsigned int branchlengths_size();
	//print
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

