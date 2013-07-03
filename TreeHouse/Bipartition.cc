//file Bipartition.cc
#include "Bipartition.h"

//Constructor
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

//Iterators
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

std::vector<float>::iterator Bipartition::branchlengths_begin() { 
	return branchlengths.begin(); 
}
std::vector<float>::iterator Bipartition::branchlengths_end() { 
	return branchlengths.end(); 
}
std::vector<float>::reverse_iterator Bipartition::branchlengths_rbegin() { 
	return branchlengths.rbegin(); 
}
std::vector<float>::reverse_iterator Bipartition::branchlengths_rend() { 
	return branchlengths.rend(); 
}
std::vector<float>::const_iterator Bipartition::branchlengths_cbegin() const{
	return branchlengths.cbegin();
}
std::vector<float>::const_iterator Bipartition::branchlengths_cend() const{
	return branchlengths.cend();
}
std::vector<float>::const_reverse_iterator Bipartition::branchlengths_crbegin() const{
	return branchlengths.crbegin();
}
std::vector<float>::const_reverse_iterator Bipartition::branchlengths_crend() const{
	return branchlengths.crend();
}

bool Bipartition::get_bit(int bitindex){
	return bitstring[bitindex];
}

float Bipartition::get_ave_branchlength(){
	float ave = 0;
	for (unsigned int i = 0; i < branchlengths.size(); i++){
		ave += branchlengths[i];
	}
	if(branchlengths.size() > 0){
		ave /= branchlengths.size();
		return ave;
	}
	return 0.0;
}

vector<float> Bipartition::get_branchlengths(){
	return branchlengths;
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

vector<unsigned int> Bipartition::get_ones_indices(){
	vector<unsigned int> indices;
	size_t index = bitstring.find_first();
	while(index != boost::dynamic_bitset<>::npos){
		indices.push_back(index);
        index = bitstring.find_next(index);
	}
	return indices;
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
	return bitstring.count();
}

int Bipartition::number_of_zeros(){
	return bitstring.size() - bitstring.count();
}


bool Bipartition::is_one(int position){
	if(in_bitstring(position) && bitstring[position] == true){
		return true;
	}
	else{ 
		return false;
	}
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
	if(is_one(position1)){
		if(is_one(position2)){
			return true;
		}
		else{
			return false;
		}
	}
	else{
		if(is_one(position2)){
			return false;
		}
	}
	return true;
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

std::set<unsigned int> Bipartition::getOnes(){
  std::set<unsigned int> retSet;
  for(unsigned int i = 0; i < bitstring.size(); i++){
	if(bitstring[i]){
		retSet.insert(i);
		}
	}
  return retSet;
}

std::set<unsigned int> Bipartition::getZeros(){
  std::set<unsigned int> retSet;
  for(unsigned int i = 0; i < bitstring.size(); i++){
	if(!bitstring[i]){
		retSet.insert(i);
		}
	}
  return retSet;
}

//size
unsigned int Bipartition::bitstring_size(){
	return bitstring.size();
}
unsigned int Bipartition::trees_size(){
	return trees.size();
}
unsigned int Bipartition::branchlengths_size(){
	return branchlengths.size();
}

//print
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


