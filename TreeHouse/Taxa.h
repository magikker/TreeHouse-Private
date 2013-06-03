#ifndef _TAXA_H
#define _TAXA_H

#include <map>
#include <vector>
#include <stdexcept>
#include <string>
#include <iostream>

using namespace std;

class Taxa {

public:
//Taxonomic ranks

unsigned int bitstringPosition;

string kingdom;
string phylum;
string class;
string order;
string family;
string genus;
string species;

};

#endif
