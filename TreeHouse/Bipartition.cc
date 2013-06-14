//file Bipartition.cc
#include "Bipartition.h"

Bipartition::Bipartition(int len){
	bitstring = boost::dynamic_bitset<>(len);
}

Bipartition::Bipartition(boost::dynamic_bitset<> bs ){
	bitstring = bs;	
}

Bipartition::Bipartition(boost::dynamic_bitset<> bs, vector<unsigned int> t, vector<float> branchl ){
	bitstring = bs;
	trees = t;
	branchlengths = branchl;
}

Bipartition::Bipartition(boost::dynamic_bitset<> bs, vector<unsigned int> t, float branchl ){
	bitstring = bs;
	trees = t;
	for (int i = 0; i < t.size(); i++){
		branchlengths.push_back(branchl);
	}
}

std::vector<unsigned int>::iterator Bipartition::trees_begin() { 
	return trees.begin(); 
}
std::vector<unsigned int>::iterator Bipartition::trees_end() { 
	return trees.end(); 
}
std::vector<unsigned int>::reverse_iterator Bipartition::trees_rbegin() { 
	return trees.rbegin(); 
}
std::vector<unsigned int>::reverse_iterator Bipartition::trees_rend() { 
	return trees.rend(); 
}
std::vector<unsigned int>::const_iterator Bipartition::trees_cbegin() const{
	return trees.cbegin();
}
std::vector<unsigned int>::const_iterator Bipartition::trees_cend() const{
	return trees.cend();
}
std::vector<unsigned int>::const_reverse_iterator Bipartition::trees_crbegin() const{
	return trees.crbegin();
}
std::vector<unsigned int>::const_reverse_iterator Bipartition::trees_crend() const{
	return trees.crend();
}

bool Bipartition::get_bit(int bitindex){
	return bitstring[bitindex];
}

float Bipartition::get_ave_branchlength(){
	float ave;
	for (unsigned int i = 0; i < branchlengths.size(); i++){
		ave += branchlengths[i];
	}
	ave /= branchlengths.size();
	return ave;
}

vector<float> Bipartition::get_branchlengths(){
	return branchlengths;
}

unsigned int Bipartition::num_trees(){
	return trees.size();
}
	
vector<unsigned int> Bipartition::get_trees(){
	return trees;
}
unsigned int Bipartition::get_tree(int index){
	return trees[index];
}

boost::dynamic_bitset<> Bipartition::get_bitstring(){
	return bitstring;
}

void Bipartition::set(int place, bool val){
	if(bitstring.size() > (unsigned)place){
		bitstring[place] = val;
	}
	else{	
		cout << "Tried to set a bitstring place that was out of range";
	}
}

void Bipartition::add_tree(int tree, float branch_len){
	trees.push_back(tree);
	branchlengths.push_back(branch_len);

}


bool Bipartition::in_bitstring(int len){
	if(bitstring.size() > (unsigned)len){
		return true;
	}
		return false;
}

int Bipartition::number_of_ones(){
	int count = 0;
	for(unsigned int i = 0; i < bitstring.size(); i ++ ){
		if (bitstring[i] == true){
			count += 1;
		}
	}
	return count;
}

int Bipartition::number_of_zeros(){
	int count = 0;
	for(unsigned int i = 0; i < bitstring.size(); i ++ ){
		if (bitstring[i] == false){
			count += 1;
		}
	}
	return count;
}

bool Bipartition::is_zero(int position){
	if( (! in_bitstring(position)) || bitstring[position] == false){
		return true;
	}
	else{ 
		return false;
	}
}

bool Bipartition::same_bitstring_value(int position1, int position2){
	bool value;
	if(in_bitstring(position1))
		value = bitstring[position1];
	else
		value = false;
	
	if(in_bitstring(position2))
		return value == bitstring[position2];
	else
		return value == false;
	cout << "shouldn't get here in function::same_bitstring_value" << endl;
	return false;
}

bool Bipartition::same_bitstring_value(vector<int> positions){
	if (is_zero(positions[0])){
		for (unsigned int i = 0; i < positions.size(); i++){
			if(is_one(positions[i])){
				return false;
			}
		}
		return true;
	}
		
	for (unsigned int i = 0; i < positions.size(); i++){
		if(is_zero(positions[i])){
			return false;
		}
	}
	return true;
}

bool Bipartition::is_one(int position){
	if(in_bitstring(position) && bitstring[position] == true){
		return true;
	}
	else{ 
		return false;
	}
}

unsigned int Bipartition::bitstring_size(){
	return bitstring.size();
}
unsigned int Bipartition::trees_size(){
	return trees.size();
}


void Bipartition::print_bitstring(bool endline){
	for (unsigned int i = 0; i < bitstring.size(); i++){
		cout << bitstring[i];
	}
	if (endline){
		cout << endl;
	}
}

void Bipartition::print_trees(bool endline){
	for (unsigned int i = 0; i < trees.size(); i++){
		cout << trees[i] << " ";
	}
	if (endline){
		cout << endl;
	}
}

void Bipartition::print_branches(bool endline){
	for (unsigned int i = 0; i < branchlengths.size(); i++){
		cout << branchlengths[i] << " ";
	}
	if (endline){
		cout << endl;
	}
}

void Bipartition::print_line(){
	print_bitstring(false);
	cout << " --> [";
	print_trees(false);
	cout << "] : [";
	print_branches(false);
	cout << "]" << endl;
}


