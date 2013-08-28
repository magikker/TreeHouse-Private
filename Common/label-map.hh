// Copyright (C) 2003, 2004 by BiRC -- Bioinformatics Research Center
//                             University of Aarhus, Denmark
//                             Contact: Thomas Mailund <mailund@birc.dk>

#ifndef LABEL_MAP_HH
#define LABEL_MAP_HH

#include <map>
#include <vector>
#include <set>
#include <stdexcept>
#include <string>

class LabelMap {
  size_t _count;
  std::map<std::string, size_t> _map;
  std::vector<std::string>  _names;

public:
  LabelMap() : _count(0) {};

  struct AlreadyPushedEx {
    std::string label;
    AlreadyPushedEx(std::string l) : label(l) {}
  };
    
  struct UnkownLabelEx {
    std::string label;
    UnkownLabelEx(std::string l) : label(l) {}
  };
  //void del(std::string label);
  size_t push(std::string label) throw(AlreadyPushedEx);
  size_t push2(std::string label);
  size_t add(std::string label);
  void rename(std::string oldl, std::string newl, unsigned int loc);
  size_t size() const { return _count; }
  void sortTaxa();
  void clear();
  
  void lookUpLabels(std::vector<std::string> names, std::vector<int> &numbers);
  std::vector<int> lookUpLabels(std::vector<std::string> names);
  std::vector<int> lookUpLabels(std::string name);
  bool areLegalTaxa(std::vector<std::string> taxanames);
  std::vector<std::string> catchDeclaredTaxa(std::vector<std::string> taxanames);
  
  std::vector<std::string> get_all_taxa_vect();
  std::set<std::string> get_all_taxa_set();
  int index_in_labelmap(std::string label);

  void printMap();
  int position(std::string label);
  size_t operator[](std::string label) throw(UnkownLabelEx);
  std::string name(unsigned int idx) const throw(std::out_of_range);
};

#endif // LABEL_MAP_HH
