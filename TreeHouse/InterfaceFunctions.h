#ifndef _INTERFACE_FUNCTIONS_HH_
#define _INTERFACE_FUNCTIONS_HH_

#include <iostream>
#include <vector>


//Global Vars like the hashtable.
#include "THGlobals.h"

//The label map class. 
#include "label-map.hh"

#include "UtilityFunctions.h"

using namespace std;

int HashCS(vector<int> input_from_int);
int HashCS(int x);
int phlash(vector<int> input_from_int);
int phlash(set<unsigned int> input_from_int);
int phlash(int x);
int quick_quartet(vector<int> input_from_int);
int quick_quartet(set<unsigned int> input_from_int);
int quick_quartet(int x);

#endif
