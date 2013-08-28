// Copyright (C) 2003, 2004 by BiRC -- Bioinformatics Research Center
//                             University of Aarhus, Denmark
//                             Contact: Thomas Mailund <mailund@birc.dk>

#include "label-map.hh"
//#include "global.h"
#include <iostream>
#include <algorithm>
#include <boost/algorithm/string/predicate.hpp>

using namespace std;

size_t
LabelMap::push(string label) throw(AlreadyPushedEx)
{
    if (_map.find(label) != _map.end()) throw AlreadyPushedEx(label);
    _names.push_back(label);  
    return _map[label] = _count++;
}

size_t
LabelMap::add(string label) {
  if (_map.find(label) == _map.end()){
    _names.push_back(label);  
    return _map[label] = _count++;
  }
  else 
    return _map[label];
}

//void
//LabelMap::del(string label) {
//  if (_map.find(label) == _map.end()){
//  }
//  else{
//	  
//		_names.erase(std::remove(_names.begin(), _names.end(), label), _names.end());
//		//_names.erase(label);
//		_map.erase(label);
//		_count--;
//  }
//}

size_t
LabelMap::operator[](string label) throw(UnkownLabelEx){
    map<string, size_t>::const_iterator i = _map.find(label);
    if (i == _map.end()) {
      //if (!HETERO || !HCHECK){
		cout << "WARNING: Used the currently unsafe LabelMap::operator[](string label)" <<endl;
		//HETERO = true;
		//NUM_TAXA++;
		size_t loc = push(label);	
		return loc;
      //}
      std::cout << "cannot find label : " << label << std::endl; 
      throw UnkownLabelEx(label);
    }
    return i->second;
}


int
LabelMap::position(string label){
    map<string, size_t>::const_iterator i = _map.find(label);
    if (i == _map.end()) {
		return -1;
    }
    return i->second;
}


std::string
LabelMap::name(unsigned int idx) const throw(std::out_of_range)
{
    return _names.at(idx);
}

void 
LabelMap::rename(string currentl, string newl, unsigned int loc)
{
  _names[loc] = newl;
  _map[newl] = _map[currentl];
}

void
LabelMap::sortTaxa(){
  _names.clear();
  size_t pos = 0;
  string label;
  map<string, size_t>::iterator i = _map.begin();
  while (i  != _map.end()){
    label = i->first;
    i->second = pos;
    _names.push_back(label);
    pos++;
    i++;
  }
}

void 
LabelMap::lookUpLabels(vector<string> names, vector<int> &numbers){
	for (unsigned int i = 0; i < names.size(); i++){
       numbers.push_back(position(names[i]));
    }
}

vector<int> 
LabelMap::lookUpLabels(vector<string> names){
	vector<int> numbers;
	for (unsigned int i = 0; i < names.size(); i++){
       numbers.push_back(position(names[i]));
    }
	return numbers;
}

vector<int>
LabelMap::lookUpLabels(string name){
	vector<int> numbers;
    numbers.push_back(position(name));
	return numbers;
}

bool 
LabelMap::areLegalTaxa(vector<string> taxanames){
	//Returns true if each taxa is already in the label map. AKA the user didn't typo a taxa name
	for (unsigned int i = 0; i < taxanames.size(); i++){
		if(position(taxanames[i]) == -1){
			return false;
		}
	}
	return true;
}

vector<string> 
LabelMap::catchDeclaredTaxa(vector<string> taxanames){
	vector<string> misses;
	for (unsigned int i = 0; i < taxanames.size(); i++){
		if(position(taxanames[i]) == -1){
			misses.push_back(taxanames[i]);
		}
	}
	return misses;
}

vector<string> LabelMap::get_all_taxa_vect(){
	vector<string> taxavect;
	for(unsigned int i = 0; i < _names.size(); i ++ ){
			taxavect.push_back(_names[i]);
	}
	return taxavect;
}

set<string> LabelMap::get_all_taxa_set(){
	set<string> taxaset;
	for(unsigned int i = 0; i < _names.size(); i ++ ){
			taxaset.insert(_names[i]);
	}
	return taxaset;
}

//upper/lower checking for taxa labels.... Mark's ancestral distance stuff uses it. . 

int LabelMap::index_in_labelmap(string label) {
  int ind = -1;
  for (unsigned int i = 0; i < _names.size(); i++) {
    if (boost::iequals(_names[i], label)) {
      ind = i;
      break;
    }
  }
  return ind;
}


void 
LabelMap::printMap(){
  map<string, size_t>::iterator i = _map.begin();
  while (i  != _map.end()){
    cout << i->first << ": " << i->second << endl;
    i++;
  }
  cout << "labels are:" << endl;
  for (unsigned int i =0; i < _names.size(); i++)
    cout << "_names[" << i << "]: " << _names[i] << endl;
}

void 
LabelMap::clear(){
  _names.clear();
  _map.clear();
}
