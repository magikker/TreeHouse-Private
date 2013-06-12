#ifndef _TAXON_H_
#define _TAXON_H_

#include <string>
#include <boost/algorithm/string/predicate.hpp>

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
    //group = to_lower(group);
    return (boost::iequals(g, group) ||
	    boost::iequals(f,group) ||
	    boost::iequals(o,group) ||
	    boost::iequals(c,group) ||
	    boost::iequals(p,group) ||
	    boost::iequals(k,group) ||
	    boost::iequals(d,group) );
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
