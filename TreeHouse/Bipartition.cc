//file Bipartition.cc
#include "Bipartition.h"

Bipartition::Bipartition(int len){
	bitstring = boost::dynamic_bitset<>(len);
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


