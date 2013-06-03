#ifndef _TAXON_H_
#define _TAXON_H_

#include <string>

using namespace std;

string to_lower(string str); //THGlobals.h includes Taxon.h before to_lower is defined.

class Taxon {

 public:
  string label;
  string d; //domain
  string k; //kingdom
  string p; //phylum
  string c; //class
  string o; //order
  string f; //family
  string g; //genus

  int pos; //position in bitstrings

  vector<int> traits;

  bool is_member(string group) {
    group = to_lower(group);
    return (to_lower(g) == group ||
	    to_lower(f) == group ||
	    to_lower(o) == group ||
	    to_lower(c) == group ||
	    to_lower(p) == group ||
	    to_lower(k) == group ||
	    to_lower(d) == group);
  }

  
	

  void print() {
    cout << label << ":";
    cout << "\n\tDomain: " << d 
	 << "\n\tKingdom: " << k 
	 << "\n\tPhylum: " << p 
	 << "\n\tClass: " << c
	 << "\n\tOrder: " << o
	 << "\n\tFamily: " << f
	 << "\n\tGenus: " << g 
	 << endl;
  }

};

#endif
