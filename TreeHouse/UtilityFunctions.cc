#include "UtilityFunctions.h"

using namespace std;

//I need some work on this one. 
set <unsigned int> duplicates(int treein){
	//Returns the unique trees in a tree set. 
	set<unsigned int> duplicatelist;
	
	if (::tree_dups[treein].size() > 0 ){
		if(::tree_dups[treein][0] < treein){
			vector<int> tempvect = ::tree_dups[::tree_dups[treein][0]];
			copy(tempvect.begin(), tempvect.end(), inserter(duplicatelist, duplicatelist.end()));
			//duplicatelist = ::dups[dups[treein][0]];
			duplicatelist.insert(::tree_dups[treein][0]);
		}
		
		else{
			vector<int> tempvect = ::tree_dups[treein];
			copy(tempvect.begin(), tempvect.end(), inserter(duplicatelist, duplicatelist.end()));
			//duplicatelist = ::dups[treein];
		}
	}
	return duplicatelist;
}

//Returns a random sampling of the trees in input set of the size requested
set <unsigned int> sample_trees (set <unsigned int> treeset, unsigned int numtrees){
	set <unsigned int> retset;
	//Used for randomizing the set
	vector<unsigned int> shuffled;
	for(std::set<unsigned int>::iterator pos = treeset.begin(); pos != treeset.end(); ++pos){
		shuffled.push_back(*pos);
	}
	//Shuffles the shuffled vector and then adds the first numtrees elements to the return set
	random_shuffle(shuffled.begin(), shuffled.end());
	for(unsigned int i = 0; i < numtrees; i++){//While we do still need trees
		retset.insert(shuffled[i]);
	}

	return retset;
}


vector<string> to_newick(vector<int> input_from_int) 
{
   vector<string> temp;
   for (unsigned int x = 0; x < input_from_int.size(); x++){
     //temp.push_back(compute_tree(::lm, ::treetable.bipartitions[x], ::treetable.branches[x], x, 0, ::treetable.bs_sizes[x]));
     temp.push_back(th_compute_tree(::biparttable, x, 0));
   }
   return temp;
}

vector<string> to_newick(set<unsigned int> inval) 
{
   vector<string> temp;
	for(set<unsigned int>::const_iterator pos = inval.begin(); pos != inval.end(); ++pos){
	  //temp.push_back(compute_tree(::lm, ::treetable.bipartitions[*pos], ::treetable.branches[*pos], *pos, 0, ::treetable.bs_sizes[*pos]));
	  temp.push_back(th_compute_tree(::biparttable, *pos, 0));
	}
   return temp;
}

string to_newick(int x) {
  //return compute_tree(::lm, ::treetable.bipartitions[x], ::treetable.branches[x], x, 0, ::treetable.bs_sizes[x]);
  return th_compute_tree(::biparttable, x, 0);
}

void show_newick(string nwstr, string title, string cssfile, string mode) {
  if (nwstr[nwstr.size()-1] != ';')
    nwstr += ";";
  if (cssfile != "")
    cssfile = "-c "+cssfile;
  if (mode == "text") {
    cout << title << ":" << endl;
    string dispcmd = "echo \""+nwstr+"\" | ../NwUtils/src/nw_display -";
    system(dispcmd.c_str());
  }
  else if (mode == "ortho") {
    string dispcmd = "echo \"<td>"+title+"\" >> temp/svg.html && "+"echo \""+nwstr+"\" | ../NwUtils/src/nw_display -s -W 10 -w 500 -l 'font-size:12' -v 20 "+cssfile+" - >> temp/svg.html";
    system(dispcmd.c_str());
  } 
  else if (mode == "radial") {
    int nlabels = num_labels_in_newick(nwstr);
    int pixwidth = 500;
    if (nlabels > 55)
	pixwidth = nlabels * 9;
    string dispcmd = "echo \"<td>"+title+"\" >> temp/svg.html && "+"echo \""+nwstr+"\" | ../NwUtils/src/nw_display -s -r -W 10 -w "+to_string(pixwidth)+" -l 'font-size:12' "+cssfile+" - >> temp/svg.html";
    system(dispcmd.c_str());
  }
}

void show_newick(vector<string> nwvect, string mode) {
  for (unsigned int i = 0; i < nwvect.size(); i++) 
    show_newick(nwvect[i], "Newick Tree "+to_string(i), "", mode);
}

void export_newick(string nwstr, string title, string path, string cssfile, string mode) {
  if (nwstr[nwstr.size()-1] != ';')
    nwstr += ";";
  if (cssfile != "")
    cssfile = "-c "+cssfile;
  if (mode == "text") {
    cout << "Creating " << path << "/" << title << ".txt" << endl;
    string dispcmd = "echo \""+nwstr+"\" | ../NwUtils/src/nw_display - > "+path+"/"+title+".txt";
    system(dispcmd.c_str());
  }
  else if (mode == "ortho") {
    cout << "Creating " << path << "/" << title << ".svg" << endl;
    string dispcmd = "echo \""+nwstr+"\" | ../NwUtils/src/nw_display -s -W 10 -w 500 -l 'font-size:12' -v 20 "+cssfile+" - > "+path+"/"+title+".svg";
    system(dispcmd.c_str());
  } 
  else if (mode == "radial") {
    cout << "Creating " << path << "/" << title << ".svg" << endl;
    int nlabels = num_labels_in_newick(nwstr);
    int pixwidth = 500;
    if (nlabels > 55)
	pixwidth = nlabels * 9;
    string dispcmd = "echo \""+nwstr+"\" | ../NwUtils/src/nw_display -s -r -W 10 -w "+to_string(pixwidth)+" -l 'font-size:12' "+cssfile+" - > "+path+"/"+title+".svg";
    system(dispcmd.c_str());
  }
}

void export_newick(vector<string> nwvect, string path, string mode) {
  for (unsigned int i = 0; i < nwvect.size(); i++) 
    export_newick(nwvect[i], "nwtree"+to_string(i), path, "", mode);
}

void show(set<unsigned int> inval, string mode) {
  for(set<unsigned int>::const_iterator pos = inval.begin(); pos != inval.end(); ++pos){
    show(*pos, mode);
  }
}

void show(vector< int> input_from_int, string mode) {
  for (unsigned int x = 0; x < input_from_int.size(); x++){
    show(input_from_int[x], mode);
  }
}

void show(int x, string mode) {
  string nwstr = to_newick(x);
  show_newick(nwstr, ("Tree "+to_string(x)), "", mode);
}

void exportt(set<unsigned int> inval, string path, string mode) {
  for(set<unsigned int>::const_iterator pos = inval.begin(); pos != inval.end(); ++pos){
    exportt(*pos, path, mode);
  }
}

void exportt(vector< int> input_from_int, string path, string mode) {
  for (unsigned int x = 0; x < input_from_int.size(); x++){
    exportt(input_from_int[x], path, mode);
  }
}

void exportt(int x, string path, string mode) {
  string nwstr = to_newick(x);
  system(("mkdir "+path+" -p").c_str());
  export_newick(nwstr, "tree"+to_string(x), path, "", mode);
}

void classification(vector<string> taxa) {
  for (unsigned int i = 0; i < taxa.size(); i++)
    classification(taxa[i]);
}

void classification(string taxon) {
  int ind = ::biparttable.lm.index_in_labelmap(taxon);
  if (ind == -1) {
    cerr << "Error: taxon name '" << taxon << "' does not exist in labelmap." << endl;
  } 
  else {
    ::taxa_info[ind]->print();
  }
}

void write_classifications(string filename) {
  ifstream isfile(filename);
  if (isfile.good()) {
    cout << "'" << filename << "' exists; overwrite? (y/n): ";
    string yn;
    getline (cin, yn);
    if (yn != "y" && yn != "Y") {
      cout << "Not overwriting '" << filename << "'. Operation aborted." << endl;
      isfile.close();
      return;
    }
  }
  isfile.close();
  ofstream cfile;
  cfile.open(filename, ios::out | ios::trunc);
  cfile << "<taxa>\n";
  for(unsigned int i = 0; i < ::taxa_info.size(); i++) {
    cfile << "  <taxon label=\"" << ::taxa_info[i]->label << "\">\n"
	  << "    <d>" << ::taxa_info[i]->d << "</d>\n"
	  << "    <k>" << ::taxa_info[i]->k << "</k>\n"
	  << "    <p>" << ::taxa_info[i]->p << "</p>\n"
	  << "    <c>" << ::taxa_info[i]->c << "</c>\n"
	  << "    <o>" << ::taxa_info[i]->o << "</o>\n"
	  << "    <f>" << ::taxa_info[i]->f << "</f>\n"
	  << "    <g>" << ::taxa_info[i]->g << "</g>\n"
	  << "    <tr>";
    for (unsigned int j=0; j<(::taxa_info[i]->traits.size()); j++) {
    cfile << ::taxa_info[i]->traits[j];
    }
    cfile << "</tr>\n"
	  << "  </taxon>\n";
  }
  cfile << "</taxa>";
  cfile.close();
}

void edit_classification(vector<string> taxa, string rank, string info) {
  std::transform(rank.begin(), rank.end(), rank.begin(), ::tolower);
  string ranks[] = {"d","k","p","c","o","f","g","domain","kingdom","phylum","class","order","family","genus"};
  set<string> rankset(ranks, ranks+14);
  if (rankset.find(rank) == rankset.end()) {
    cerr << "Error: '" << rank << "' does not denote a valid rank." << endl;
  } 
  else {
    for (unsigned int i = 0; i < taxa.size(); i++) {
      edit_classification(taxa[i], rank, info);
    }
  }
}

void edit_classification(string taxon, string rank, string info) {
  std::transform(rank.begin(), rank.end(), rank.begin(), ::tolower);
  //rank = to_lower(rank);
  int ind = ::biparttable.lm.index_in_labelmap(taxon);
  if (ind == -1) {
    cerr << "Error: taxon name '" << taxon << "' does not exist in labelmap." << endl;
    return;
  }
  if (rank == "d" || rank == "domain") ::taxa_info[ind]->d = info;
  else if (rank == "k" || rank == "kingdom") ::taxa_info[ind]->k = info;
  else if (rank == "p" || rank == "phylum") ::taxa_info[ind]->p = info;
  else if (rank == "c" || rank == "class") ::taxa_info[ind]->c = info;
  else if (rank == "o" || rank == "order") ::taxa_info[ind]->o = info;
  else if (rank == "f" || rank == "family") ::taxa_info[ind]->f = info;
  else if (rank == "g" || rank == "genus") ::taxa_info[ind]->g = info;
  else {
    cerr << "Error: '" << rank << "' does not denote a valid rank." << endl;
    return;
  }
}

bool compare_padded_taxon_names(Taxon* taxon1, Taxon* taxon2) {
  string name1 = taxon1->label;
  string name2 = taxon2->label;
  while (name1.size() < 100)
    name1 = "0" + name1;
  while (name2.size() < 100)
    name2 = "0" + name2;
  return name1 < name2;
}
  
void load_trait_matrix(string file) {
  //if numerical labels are used, get the ordering of the taxa used in the matrix file
  vector<Taxon *> taxa_info_ordered = ::taxa_info;
  std::sort(taxa_info_ordered.begin(), taxa_info_ordered.end(), compare_padded_taxon_names);
  ifstream stream;
  bool foundmatrix = 0;
  string line;
  int taxon = 0;
  stream.open(file);
  if (stream.is_open()) {
    while (stream.good() && line != ";" && taxon < taxa_info_ordered.size()) {
      getline (stream,line);
      cout << line << endl;
      if (line == "\tMATRIX") {
	foundmatrix = 1;
	continue;
      }
      if (foundmatrix) {
	for (std::string::reverse_iterator rit=line.rbegin(); rit!=line.rend(); ++rit) {
	  if (*rit == ' ') {
	    break;
	  }
	  string tempstr;
	  tempstr.push_back(*rit);
	  cout << tempstr << endl;
	  taxa_info_ordered[taxon]->traits.insert(taxa_info_ordered[taxon]->traits.begin(), stoi(tempstr));
	}
	taxon++;
      }
    }
  }
}
	    
	  
      
  

vector<string> taxa_in_group(string group) {
  vector<string> taxa;
  for (unsigned int i = 0; i < ::taxa_info.size(); i++) {
    if (::taxa_info[i]->is_member(group))
      taxa.push_back(taxa_info[i]->label);
  }
  return taxa;
}

vector<string> groups_at_level(string rank) {
  std::transform(rank.begin(), rank.end(), rank.begin(), ::tolower);
  //rank = to_lower(rank);
  vector<string> groupvect;
  set<string> groupset;
  for (unsigned int i = 0; i < ::taxa_info.size(); i++) {
    if (rank == "d" || rank == "domain") groupset.insert(::taxa_info[i]->d);
    else if (rank == "k" || rank == "kingdom") groupset.insert(::taxa_info[i]->k);
    else if (rank == "p" || rank == "phylum") groupset.insert(::taxa_info[i]->p);
    else if (rank == "c" || rank == "class") groupset.insert(::taxa_info[i]->c);
    else if (rank == "o" || rank == "order") groupset.insert(::taxa_info[i]->o);
    else if (rank == "f" || rank == "family") groupset.insert(::taxa_info[i]->f);
    else if (rank == "g" || rank == "genus") groupset.insert(::taxa_info[i]->g);
    else {
      cerr << "Error: '" << rank << "' does not denote a valid rank." << endl;
      return groupvect;
    }
  }
  for (set<string>::const_iterator pos = groupset.begin(); pos != groupset.end(); pos++) {
    groupvect.push_back(*pos);
  }
  return groupvect;
}
    
vector<string> parse_newick(string newick) {
  vector<string> parsed_newick;
  string label;
  for (unsigned int i = 0; i < (newick.size()-1); i++) {
    if (newick[i] == '(') {
      if (label != "") {
	parsed_newick.push_back(label);
	label = "";
      }
      parsed_newick.push_back("(");
    } 
    else if (newick[i] == ')') {
      if (label != "") {
	parsed_newick.push_back(label);
	label = "";
      }
      parsed_newick.push_back(")");
    } 
    else if (newick[i] == ',') {
      parsed_newick.push_back(label);
      label = "";
    }
    else {
      label.push_back(newick[i]);
    }
  }
  return parsed_newick;
}

int num_labels_in_newick(string newick) {
  vector<string> parsed_newick = parse_newick(newick);
  int num_labels = 0;
  for (unsigned int i = 0; i < (parsed_newick.size()-1); i++) {
    if (parsed_newick[i] != "(" && parsed_newick[i] != ")")
      num_labels++;
  }
  return num_labels;
}

bool is_elem_of_vector(string item, vector<string> vect) {
  for (unsigned int i = 0; i < vect.size(); i++) {
    if (vect[i] == item)
      return true;
  }
  return false;
}

string get_clade_of_group(string group, string newick) {
  string taxastring;
  vector<string> taxa = taxa_in_group(group);
  for (unsigned int i = 0; i < taxa.size(); i++) {
    taxastring += (" " + taxa[i]);
  }
  system(("echo \""+newick+"\" | ../NwUtils/src/nw_clade -"+taxastring+" > temp/clade.txt").c_str());
  ifstream cladefile ("temp/clade.txt");
  string clade;
  getline (cladefile, clade);
  cladefile.close();
  return clade;
}

void highlight_group_css(string group, vector<string> taxa, string newick, int treeid, map<string, vector <int> > &out_of_place) {
  string taxastring;
  for (unsigned int i = 0; i < taxa.size(); i++) {
    taxastring += (" " + taxa[i]);
  }
  string clade = get_clade_of_group(group, newick);
  vector<string> parsed_clade = parse_newick(clade);
  vector<string> foreign_taxa;
  for (unsigned int i = 0; i < parsed_clade.size(); i++) {
    if (parsed_clade[i] != "(" && parsed_clade[i] != ")" && !is_elem_of_vector(parsed_clade[i], taxa))
      foreign_taxa.push_back(parsed_clade[i]);
  }
  system(("echo '\"stroke-width:2; stroke:blue\" Clade"+taxastring+"' > temp/group.css").c_str());
  for (unsigned int i = 0; i < foreign_taxa.size(); i++) {
    system(("echo '\"stroke-width:2; stroke:red\" Individual "+foreign_taxa[i]+"' >> temp/group.css").c_str());
    out_of_place[foreign_taxa[i]].push_back(treeid);
  }
}

int highlight_group_taxa( string group, vector<string> taxa, int tree, string mode, map<string, vector <int> > &out_of_place, bool isolate/*=false*/)  {
  string nwstr;
  if (isolate)
    nwstr = get_clade_of_group(group, to_newick(tree));
  else
    nwstr = to_newick(tree);
  highlight_group_css(group, taxa, nwstr, tree, out_of_place);
  show_newick(nwstr, (group+" in Tree "+to_string(tree)), "temp/group.css", mode);
  return 0;
}

int highlight_group(string group, set<unsigned int> treeset, string mode, bool isolate/*=false*/) {
  vector<string> taxa = taxa_in_group(group);
  if (taxa.size() == 0) {
    cerr << "Error: no taxa in group '"+group+"'" << endl;
    return 1;
  }
  map<string, vector <int> > out_of_place;
  for (set<unsigned int>::const_iterator pos = treeset.begin(); pos != treeset.end(); pos++) {
    highlight_group_taxa(group, taxa, *pos, mode, out_of_place, isolate);
    system("rm temp/group.css");
  }
  //put the out-of-place map pairs accrued in highlight_group into a vector.
  vector< pair <string, vector <int> > > out_of_place_vect;
  map<string, vector <int> >::iterator it;
  for (it = out_of_place.begin(); it != out_of_place.end(); it++) {
    out_of_place_vect.push_back(*it);
  }
   //sort the vector by its pairs' second.size()
  sort(out_of_place_vect.begin(), out_of_place_vect.end(),
       [](const pair< string, vector <int> > &a, const pair< string, vector <int> > &b) -> bool {
	 return a.second.size() > b.second.size();
	   });
  for (unsigned int i = 0; i < out_of_place_vect.size(); i++) {
    string report = out_of_place_vect[i].first+" was placed in "+group+" by [";
    for (unsigned int j = 0; j < out_of_place_vect[i].second.size()-1; j++) {
      report += (to_string(out_of_place_vect[i].second[j]) + ", ");
    }
    report += (to_string(out_of_place_vect[i].second.back()) + "] ");
    double pct = ((double)(out_of_place_vect[i].second.size()) / (double)(treeset.size())) * 100;
    string pctstr = to_string(pct);
    report += ("(" + pctstr + "%)");
    system("echo \"</table>\" >> temp/svg.html");
    system(("echo '"+report+"' '</br>' >> temp/svg.html").c_str());
  }
  return 0;
 }
 
int highlight_group(string group, vector<int> treevect, string mode, bool isolate/*=false*/) {
  set<unsigned int> treeset(treevect.begin(), treevect.end());
  return highlight_group(group, treeset, mode, isolate);
 }

//int show_clade(vector<string> taxa, set<unsigned int> treeset, string mode, bool isolate/*=false*/) {
//string taxastring = "[";
//for (unsigned int i = 0; i < taxa.size() - 1; i++)
//  taxastring += (taxa[i]+", ");
//taxastring += (taxa.back()+"]");
//ret highlight_group_taxa(taxastring, taxa, 

int show_level(string rank, int tree, string mode) {
  vector<string> groups = groups_at_level(rank);
  if (groups.size() == 0)
    return 1;
  string nwstr = to_newick(tree);
  vector<string> clades;
  for (unsigned int i = 0; i < groups.size(); i++) {
    string clade = get_clade_of_group(groups[i], nwstr);
    clade.resize(clade.size() - 1); //erase semicolon at end of string
    clades.push_back(clade);
  }
  for (unsigned int i = 0; i < clades.size(); i++) {
    int pos = nwstr.find(clades[i]);
    if (pos == -1) {
      cerr << "Error: Tree "+to_string(tree)+" has conflicts at the "+rank+" level." << endl;
      return 1;
    }
    nwstr.replace(pos, clades[i].length(), groups[i]);
  }
  show_newick(nwstr, ("Tree "+to_string(tree)+" at "+rank+" level"), "", mode);
  return 0;
}

int show_level(string rank, set<unsigned int> treeset, string mode) {
  for(set<unsigned int>::const_iterator pos = treeset.begin(); pos != treeset.end(); pos++){
   int ret = show_level(rank, *pos, mode);
   if (ret == 1)
     return 1;
  }
  return 0;
}

int show_level(string rank, vector<int> treevect, string mode) {
  for(unsigned int i = 0; i < treevect.size(); i++) {
    int ret = show_level(rank, treevect[i], mode);
    if (ret == 1)
      return 1;
  }
  return 0;
}

int show_only_taxa(string group, vector<string> taxa, int tree, string mode) {
  set<int> all_positions;
  for (unsigned int i = 0; i < ::biparttable.lm.size(); i++) {
    all_positions.insert(i);
  }
  
  set<int> in_positions;
  for (unsigned int i = 0; i < taxa.size(); i++) {
    in_positions.insert(::biparttable.lm.position(taxa[i]));
  }
  
  vector<int> out_positions;
  std::set_difference(all_positions.begin(), all_positions.end(), in_positions.begin(), in_positions.end(), std::back_inserter(out_positions));

  vector<float> branches;
  vector<boost::dynamic_bitset<> > tree_bipartitions;
  ::get_tree_data(tree, tree_bipartitions, branches);
  vector<boost::dynamic_bitset<> > bipartitions;

  for (unsigned int i = 0; i < tree_bipartitions.size(); i++) {
    // copy bipartition from treetable into new bool array; push onto bipartitions vector
    boost::dynamic_bitset<> bitstring(tree_bipartitions[i].size());
    for (unsigned int j = 0; j < tree_bipartitions[i].size(); j++) {
      bitstring[j] = tree_bipartitions[i][j];
    }
    bipartitions.push_back(bitstring);
    // mark non-group taxa with 0 in all bipartitions
    for (unsigned int j = 0; j < out_positions.size(); j++) {
      if (out_positions[j] < bipartitions[i].size())
	bipartitions[i][out_positions[j]] = 0;
    }
  }
  //remove bipartitions that are now all zeroes or duplicates
  set<string> bipartset;
  //vector<unsigned int>::reverse_iterator bs_size = bs_sizes.rbegin();
  vector<float>::reverse_iterator branch = branches.rbegin();
  for (vector< boost::dynamic_bitset<> >::reverse_iterator bipart = bipartitions.rbegin();
       bipart < bipartitions.rend();
       bipart++, branch++) {
    string bitstring;
    bool has_one = false;
    for (unsigned int j = 0; j < (*bipart).size(); j++) {
      bitstring += to_string((*bipart)[j]);
      if ((*bipart)[j] == 1) {
	     has_one = true;
      }
    }
    if (has_one) { //remove trailing zeroes
      while (bitstring[bitstring.size() - 1] == '0') {
    	bitstring.resize(bitstring.size() - 1);
      }
    }
    if (!has_one || !(bipartset.insert(bitstring)).second) {
      bipartitions.erase(bipart.base() - 1);
      if (branches.size() > 0)
	branches.erase(branch.base() - 1);
    }
  }
  
  string nwstr = compute_tree(::biparttable.lm, bipartitions, branches, tree, 0);
  cout << nwstr << endl;
  show_newick(nwstr, (group+" in Tree "+to_string(tree)), "", mode);
  return 0;
}
	      	       
int show_only(string group, set<unsigned int> treeset, string mode) {
  vector<string> taxa = taxa_in_group(group);
  if (taxa.size() == 0) {
    cerr << "Error: no taxa in group "+group+"." << endl;
    return 1;
  }
  for(set<unsigned int>::const_iterator pos = treeset.begin(); pos != treeset.end(); pos++){
    show_only_taxa(group, taxa, *pos, mode);
  }
  return 0;
}

int show_only(string group, vector<int> treevect, string mode) {
  set<unsigned int> treeset(treevect.begin(), treevect.end());
  return show_only(group, treeset, mode);
}

int show_only(vector<string> groupvect, set<unsigned int> treeset, string mode) {
  set<string> taxaset;
  string groupname = "[";
  for (unsigned int i = 0; i < groupvect.size(); i++) {
    vector<string> taxa = taxa_in_group(groupvect[i]);
    if (taxa.size() == 0) {
      cerr << "Error: no taxa in group "+groupvect[i]+"." << endl;
      return 1;
    }
    for(unsigned int j = 0; j < taxa.size(); j++) {
      taxaset.insert(taxa[j]);
    }
    if (i < groupvect.size() - 1) 
      groupname += (groupvect[i]+", ");
    else 
      groupname += (groupvect[i]+"]");
  }
  vector<string> all_taxa(taxaset.begin(), taxaset.end());
  for(set<unsigned int>::const_iterator pos = treeset.begin(); pos != treeset.end(); pos++){
    show_only_taxa(groupname, all_taxa, *pos, mode);
  }
  return 0;
}

int show_only(vector<string> groupvect, vector<int> treevect, string mode) {
  set<unsigned int> treeset(treevect.begin(), treevect.end());
  return show_only(groupvect, treeset, mode);
}
//This is being replaced to a great degree... 
/*
void taxa_filter (vector<string> taxa, unsigned int tree) {
  unsigned int new_tree_id = ::NUM_TREES;
  set<int> all_positions;
  for (unsigned int i = 0; i < ::NUM_TAXA; i++) {
    all_positions.insert(i);
  }
  set<int> in_positions;
  for (unsigned int i = 0; i < taxa.size(); i++) {
    in_positions.insert(::lm.position(taxa[i]));
  }
  vector<int> out_positions;
  std::set_difference(all_positions.begin(), all_positions.end(), in_positions.begin(), in_positions.end(), std::back_inserter(out_positions));
  vector<unsigned int> bs_sizes;
  vector<float> branches;
  vector<bool *> tree_bipartitions;
  vector<unsigned int> bipart_indices = ::get_tree_data(tree, tree_bipartitions, bs_sizes, branches);
  vector<bool *> bipartitions;
  for (unsigned int i = 0; i < tree_bipartitions.size(); i++) {
    // copy bipartition from treetable into new bool array; push onto bipartitions vector
    bool * bitstring  = new bool [bs_sizes[i]];
    for (unsigned int j = 0; j < bs_sizes[i]; j++) {
      bitstring[j] = tree_bipartitions[i][j];
    }
    bipartitions.push_back(bitstring);
    // mark non-group taxa with 0 in all bipartitions
    for (unsigned int j = 0; j < out_positions.size(); j++) {
      if (out_positions[j] < bs_sizes[i])
	bipartitions[i][out_positions[j]] = 0;
    }
  }
  //remove bipartitions that are now all zeroes or duplicates
  set<string> bipartset;
  vector<unsigned int>::reverse_iterator bs_size = bs_sizes.rbegin();
  vector<float>::reverse_iterator branch = branches.rbegin();
  vector<unsigned int>::reverse_iterator bipart_ind = bipart_indices.rbegin();
  for (vector<bool *>::reverse_iterator bipart = bipartitions.rbegin();
       bipart < bipartitions.rend();
       bipart++, bs_size++, branch++, bipart_ind++) {
    string bitstring;
    bool has_one = false;
    for (unsigned int j = 0; j < *bs_size; j++) {
      bitstring += to_string((*bipart)[j]);
      if ((*bipart)[j] == 1) {
	has_one = true;
      }
    }
    if (has_one) { //remove trailing zeroes
      while (bitstring[bitstring.size() - 1] == '0') {
	bitstring.resize(bitstring.size() - 1);
	(*bs_size)--;
      }
    }
    if (!has_one || !(bipartset.insert(bitstring)).second) {
      bipartitions.erase(bipart.base() - 1);
      bs_sizes.erase(bs_size.base() - 1);
      if (branches.size() > 0)
	branches.erase(branch.base() - 1);
      bipart_indices.erase(bipart_ind.base() - 1);
    }
  }
  ::biparttable.tree_branches.push_back(branches);
  for (unsigned int i = 0; i < bipartitions.size(); i++) {
    bool new_bitstring = false;
    unsigned int insert_ind = bipart_indices[i];
    if (bs_sizes[i] != ::biparttable.bitstring_size(bipart_indices[i]) )
      new_bitstring = true;
    if (bs_sizes[i] == ::biparttable.bitstring_size(bipart_indices[i])) {
      for (unsigned int j = 0; j < bs_sizes[i]; j++) {
	if (bipartitions[i][j] != ::biparttable.get_bit(bipart_indices[i], j) ) {
	  new_bitstring = true;
	  break;
	}
      }
    }
    if (new_bitstring && bipart_indices[i] > 0 && bs_sizes[i] == ::biparttable.bitstring_size(bipart_indices[i]-1)) {
      //created bipartitions go directly above their original counterparts; we have to check those too
      new_bitstring = false;
      insert_ind--;
      for (unsigned int j = 0; j < bs_sizes[i]; j++) {
	if (bipartitions[i][j] != ::biparttable.get_bit(bipart_indices[i]-1, j)) {
	  new_bitstring = true;
	  insert_ind++;
	  break;
	}
      }
    }
    if (new_bitstring) {
      ::biparttable.bipartitions.insert(::biparttable.bipartitions.begin() + insert_ind, bipartitions[i]);
      ::biparttable.length_of_bitstrings.insert(::biparttable.length_of_bitstrings.begin() + insert_ind, bs_sizes[i]);
      vector<unsigned int> search_entry;
      search_entry.push_back(new_tree_id);
      ::biparttable.searchtable.insert(::biparttable.searchtable.begin() + insert_ind, search_entry);
      for (unsigned int j = 0; j < bipart_indices.size(); j++) {
	bipart_indices[j]++;
      }
      bool * treestr = new bool [::NUM_TREES];
      ::biparttable.treetable.insert(::biparttable.treetable.begin() + insert_ind, treestr);
    }
    if (!new_bitstring) {
      delete[] bipartitions[i];
      ::biparttable.searchtable[insert_ind].push_back(new_tree_id);
    }
  }
  ::NUM_TREES++;
  for (unsigned int i = 0; i < ::biparttable.bipartitions.size(); i++) {
      taxa_in_tree[i] = 0;
  }
  ::taxa_in_trees.push_back(taxa_in_tree);
  //recalculate tree_dups
  vector<int> new_dups;
  for (unsigned int i = 0; i < ::NUM_TREES - 1; i++) {
    bool is_dup = true;
    for (unsigned int j = 0; j < ::biparttable.bipartitions.size(); j++) {
      is_dup = !(::biparttable.treetable[j][i] ^ ::biparttable.treetable[j][new_tree_id]);
      if (!is_dup)
	break;
    }
    if (is_dup) {
      new_dups.push_back(i);
      ::tree_dups[i].push_back(new_tree_id);
    }
  }
  ::tree_dups[new_tree_id] = new_dups;
  ::all_trees.insert(new_tree_id);
  symbol_table["trees"] = new pqlsymbol(all_trees, (unsigned int)::NUM_TREES);
  cout << "Created Tree " << new_tree_id << "." << endl;
}

void taxa_filter (vector<string> taxa, set<unsigned int> treeset) {
  for (set<unsigned int>::const_iterator pos = treeset.begin(); pos != treeset.end(); pos++)
    taxa_filter(taxa, *pos);
}

void taxa_filter(vector<string> taxa, vector<int> treevect) {
  for (unsigned int i = 0; i < treevect.size(); i++)
    taxa_filter(taxa, treevect[i]);
}

void group_filter(vector<string> groups, unsigned int tree) {
  set<string> all_taxa_set;
  for (unsigned int i = 0; i < groups.size(); i++) {
    vector<string> taxa = taxa_in_group(groups[i]);
    if (taxa.size() == 0)
      cerr << "Warning: no taxa found in group '"+groups[i]+"'." << endl;
    for (unsigned int j = 0; j < taxa.size(); j++) {
      all_taxa_set.insert(taxa[j]);
    }
  }
  vector<string> all_taxa(all_taxa_set.begin(), all_taxa_set.end());
  taxa_filter(all_taxa, tree);
}

void group_filter (vector<string> groups, set<unsigned int> treeset) {
  for (set<unsigned int>::const_iterator pos = treeset.begin(); pos != treeset.end(); pos++)
    group_filter(groups, *pos);
}

void group_filter(vector<string> groups, vector<int> treevect) {
  for (unsigned int i = 0; i < treevect.size(); i++)
    group_filter(groups, treevect[i]);
}

void delete_tree(unsigned int tree) {
  delete[] ::taxa_in_trees[tree];
  ::taxa_in_trees.erase(::taxa_in_trees.begin() + tree);
  ::biparttable.tree_branches.erase(::biparttable.tree_branches.begin() + tree);
  ::tree_dups.erase(tree);
  for (map<int, vector<int> >::iterator pos = ::tree_dups.begin(); pos != tree_dups.end(); pos++) {
    for (unsigned int i = 0; i < (*pos).second.size(); i++) {
      if ((*pos).second[i] == tree)
	(*pos).second.erase((*pos).second.begin() + i);
    }
  }
  for (unsigned int i = 0; i < ::biparttable.bipartitions.size(); i++) {
    bool increment = true;
    for (unsigned int j = tree; j+1 < ::NUM_TREES; j++) {
      ::biparttable.treetable[i][j] = ::biparttable.treetable[i][j+1];
    }
    for (unsigned int j = 0; j < ::biparttable.searchtable[i].size(); j++) {
      if (::biparttable.searchtable[i][j] == tree) { //if this tree had bipartition i
    delete[] ::biparttable.treetable[i];
    bool * treestr = new bool [::NUM_TREES];
    for (unsigned int j = 0; j < ::NUM_TREES; j++) {
      treestr[j] = 0;
    }
    for (unsigned int j = 0; j < ::biparttable.searchtable[i].size(); j++) {
      treestr[::biparttable.searchtable[i][j]] = 1;
    }
    ::biparttable.treetable[i] = treestr;
  }
  bool * taxa_in_tree = new bool [::NUM_TAXA];
  for (unsigned int i = 0; i < ::NUM_TAXA; i++) {
    if (in_positions.find(i) != in_positions.end())
      taxa_in_tree[i] = 1;
    else
      taxa_in_tree[i] = 0;
  }
  ::taxa_in_trees.push_back(taxa_in_tree);
  //recalculate tree_dups
  vector<int> new_dups;
  for (unsigned int i = 0; i < ::NUM_TREES - 1; i++) {
    bool is_dup = true;
    for (unsigned int j = 0; j < ::biparttable.bipartitions.size(); j++) {
      is_dup = !(::biparttable.treetable[j][i] ^ ::biparttable.treetable[j][new_tree_id]);
      if (!is_dup)
	break;
    }
    if (is_dup) {
      new_dups.push_back(i);
      ::tree_dups[i].push_back(new_tree_id);
    }
  }
  ::tree_dups[new_tree_id] = new_dups;
  ::all_trees.insert(new_tree_id);
  symbol_table["trees"] = new pqlsymbol(all_trees, (unsigned int)::NUM_TREES);
  cout << "Created Tree " << new_tree_id << "." << endl;
}

void taxa_filter (vector<string> taxa, set<unsigned int> treeset) {
  for (set<unsigned int>::const_iterator pos = treeset.begin(); pos != treeset.end(); pos++)
    taxa_filter(taxa, *pos);
}

void taxa_filter(vector<string> taxa, vector<int> treevect) {
  for (unsigned int i = 0; i < treevect.size(); i++)
    taxa_filter(taxa, treevect[i]);
}

void group_filter(vector<string> groups, unsigned int tree) {
  set<string> all_taxa_set;
  for (unsigned int i = 0; i < groups.size(); i++) {
    vector<string> taxa = taxa_in_group(groups[i]);
    if (taxa.size() == 0)
      cerr << "Warning: no taxa found in group '"+groups[i]+"'." << endl;
    for (unsigned int j = 0; j < taxa.size(); j++) {
      all_taxa_set.insert(taxa[j]);
    }
  }
  vector<string> all_taxa(all_taxa_set.begin(), all_taxa_set.end());
  taxa_filter(all_taxa, tree);
}

void group_filter (vector<string> groups, set<unsigned int> treeset) {
  for (set<unsigned int>::const_iterator pos = treeset.begin(); pos != treeset.end(); pos++)
    group_filter(groups, *pos);
}

void group_filter(vector<string> groups, vector<int> treevect) {
  for (unsigned int i = 0; i < treevect.size(); i++)
    group_filter(groups, treevect[i]);
}

void delete_tree(unsigned int tree) {
  delete[] ::taxa_in_trees[tree];
  ::taxa_in_trees.erase(::taxa_in_trees.begin() + tree);
  ::biparttable.tree_branches.erase(::biparttable.tree_branches.begin() + tree);
  ::tree_dups.erase(tree);
  for (map<int, vector<int> >::iterator pos = ::tree_dups.begin(); pos != tree_dups.end(); pos++) {
    for (unsigned int i = 0; i < (*pos).second.size(); i++) {
      if ((*pos).second[i] == tree)
	(*pos).second.erase((*pos).second.begin() + i);
    }
  }
  for (unsigned int i = 0; i < ::biparttable.bipartitions.size(); i++) {
    bool increment = true;
    for (unsigned int j = tree; j+1 < ::NUM_TREES; j++) {
      ::biparttable.treetable[i][j] = ::biparttable.treetable[i][j+1];
    }
    for (unsigned int j = 0; j < ::biparttable.searchtable[i].size(); j++) {
      if (::biparttable.searchtable[i][j] == tree) { //if this tree had bipartition i
	//set the searchtable and treetable to no longer include the tree
	::biparttable.searchtable[i].erase(::biparttable.searchtable[i].begin() + j);
	j--;
	if (::biparttable.searchtable[i].size() == 0) { //if no trees have this bipartition any longer
	  //remove bipartition entirely
	  delete[] ::biparttable.bipartitions[i];
	  ::biparttable.bipartitions.erase(::biparttable.bipartitions.begin() + i);
	  ::biparttable.length_of_bitstrings.erase(::biparttable.length_of_bitstrings.begin() + i);
	  ::biparttable.searchtable.erase(::biparttable.searchtable.begin() + i);
	  delete[] ::biparttable.treetable[i];
	  ::biparttable.treetable.erase(::biparttable.treetable.begin() + i);
	  if (biparttable.branches.size() > i) 
	    ::biparttable.branches.erase(::biparttable.branches.begin() + i);
	  i--;
	  break;
	}
      }
      else if (::biparttable.searchtable[i][j] > tree) //change searchtable to account for new tree ids
	::biparttable.searchtable[i][j]--;
    }
  }
  ::NUM_TREES--;
  all_trees.clear();
  for (unsigned int i = 0; i < ::NUM_TREES; i++) {
    all_trees.insert(i);
  }
  symbol_table["trees"] = new pqlsymbol(all_trees, (unsigned int)::NUM_TREES);
  if (original_trees.erase(tree)) {
    ::NUM_TREES_INIT--;
    for (unsigned int i = 0; i < ::NUM_TREES_INIT; i++) {
      original_trees.insert(i);
    }
    symbol_table["original_trees"] = new pqlsymbol(original_trees, (unsigned int)::NUM_TREES_INIT);
  }
  cout << "Tree erased." << endl;
}

void delete_tree(set<unsigned int> treeset) {
  int offset = 0;
  for (set<unsigned int>::const_iterator pos = treeset.begin(); pos != treeset.end(); pos++) {
    cout << (*pos) << " -> ";
    delete_tree(*pos - offset);
    offset++;
  }
}

void delete_tree(vector<int> treevect) {
  set<unsigned int> treeset(treevect.begin(), treevect.end());
    delete_tree(treeset);
}
*/
void write_trz(vector<string> nwvect, string filename) {
  filename += ".tre";
  ifstream isfile(filename);
  if (isfile.good()) {
    cout << "'" << filename << "' exists; overwrite? (y/n): ";
    string yn;
    getline (cin, yn);
    if (yn != "y" && yn != "Y") {
      cout << "Not overwriting '" << filename << "'. Operation aborted." << endl;
      isfile.close();
      return;
    }
  }
  isfile.close();
  ofstream cfile;
  cfile.open(filename, ios::out | ios::trunc);
  for (unsigned int i = 0; i < nwvect.size(); i++) {
    cfile << nwvect[i] << endl;
  }
  cfile.close();
  cout << "Created " << filename << endl;
  system(("../TreeZip/treezip "+filename).c_str());
}

void write_trz(int tree, string filename) {
  vector<string> nwvect;
  nwvect.push_back(to_newick(tree));
  write_trz(nwvect, filename);
}

void write_trz(vector<int> treevect, string filename) {
  vector<string> nwvect;
  for (unsigned int i = 0; i < treevect.size(); i++)
    nwvect.push_back(to_newick(treevect[i]));
  write_trz(nwvect, filename);
}

void write_trz(set<unsigned int> treeset, string filename) {
  vector<string> nwvect;
  for (set<unsigned int>::const_iterator pos = treeset.begin(); pos != treeset.end(); pos++)
    nwvect.push_back(to_newick(*pos));
  write_trz(nwvect, filename);
}

set <unsigned int> unique(set< unsigned int> treesin){
	//cout << "unique has been called" <<endl;
	//Returns the unique trees in a tree set. 
	set<unsigned int> uniquetrees = treesin;
	
	for(set<unsigned int>::const_iterator pos = treesin.begin(); pos != treesin.end(); ++pos){
		for (unsigned int j = 0; j < ::biparttable.tree_dups[*pos].size(); j++){
			if (*pos < ::biparttable.tree_dups[*pos][j]){ // if the dup is greater in id number than the tree we are on delete it
				uniquetrees.erase(::biparttable.tree_dups[*pos][j]);
			}
			else if (*pos > ::biparttable.tree_dups[*pos][j]){ // if the dup is less than the tree we are on the we need to do some special ops... we get a new list based on the lowest number and we delete only those that are greater than the current value. 
				//cout << "in the else if" << endl;
				vector<int> tempvect = ::biparttable.tree_dups[::biparttable.tree_dups[*pos][j]];
				for (unsigned int k = 0; k < tempvect.size(); k++){
					if (*pos < tempvect[k]){ // if the dup is greater in id number than the tree we are on delete it
						uniquetrees.erase(tempvect[k]);
					}
				}
			}			
		}
	}
	return uniquetrees;
}

//GRB we don't need to pass around vectors of bool * to get a number, we can use bipartition ids
//we can also just use the inverted index to get the bipartition ids directly. 
int unique_biparts(set< unsigned int > treesin){
//	cout << "unique biparts has been called" << endl;
	//returns the unique bipartitions contained in a set of trees
	set<unsigned int> retSet; //the set we are going to return

	for(int i = 0; i < treesin.size(); i++){ //iterate through all of the input trees
		vector<unsigned int> tree_biparts = ::biparttable.inverted_index[i];
		for(int j = 0; j < tree_biparts.size(); j++){			
			retSet.insert(tree_biparts[j]);	  
		}
	}
	return retSet.size();
}
//when returning trees, I'm going to try to return sets in general. So lets phase the vector versions
/*
vector <unsigned int> unique(vector<unsigned int> treesin){
	//Returns the unique trees in a tree set. 
	vector<unsigned int> uniquetrees = treesin;
		
	for (unsigned int i = 0; i < treesin.size(); i++){ // for each tree in the input set
		for (unsigned int j = 0; j < ::biparttable.tree_dups[treesin[i]].size(); j++){ // check to see if it has dups. 
			//cout << "treesin[i] = " << treesin[i] << endl;
			//cout << "::dups[treesin[i]][j] = " << ::dups[treesin[i]][j] << endl;

			if (treesin[i] < ::biparttable.tree_dups[treesin[i]][j]){ // if the dup is greater in id number than the tree we are on delete it
				uniquetrees.erase(std::remove(uniquetrees.begin(), uniquetrees.end(), ::biparttable.tree_dups[treesin[i]][j]), uniquetrees.end());
			}
			
			else if (treesin[i] > ::biparttable.tree_dups[treesin[i]][j]){ // if the dup is less than the tree we are on the we need to do some special ops... we get a new list based on the lowest number and we delete only those that are greater than the current value. 
				//cout << "in the else if" << endl;
				vector<int> tempvect = ::biparttable.tree_dups[::biparttable.tree_dups[treesin[i]][j]];
				for (unsigned int k = 0; k < tempvect.size(); k++){
					if (treesin[i] < tempvect[k]){ // if the dup is greater in id number than the tree we are on delete it
						uniquetrees.erase(std::remove(uniquetrees.begin(), uniquetrees.end(), tempvect[k]), uniquetrees.end());
					}
				}
			}
    //cout << endl;
		}
	}
	//~ for (unsigned int i = 0; i < uniquetrees.size(); i++){
		//~ cout << uniquetrees[i];
	//~ }
	//~ cout << endl;

	return uniquetrees;
}
*/
/*
vector <int> unique(vector<int> treesin){
	//cout << "unique has been called" <<endl;
	//Returns the unique trees in a tree set. 
	vector<int> uniquetrees = treesin;
		
	for (unsigned int i = 0; i < treesin.size(); i++){ // for each tree in the input set
		for (unsigned int j = 0; j < ::tree_dups[treesin[i]].size(); j++){ // check to see if it has dups. 
			//cout << "treesin[i] = " << treesin[i] << endl;
			//cout << "::dups[treesin[i]][j] = " << ::dups[treesin[i]][j] << endl;
			if (treesin[i] < ::tree_dups[treesin[i]][j]){ // if the dup is greater in id number than the tree we are on delete it
				uniquetrees.erase(std::remove(uniquetrees.begin(), uniquetrees.end(), ::tree_dups[treesin[i]][j]), uniquetrees.end());
			}
	
			else if (treesin[i] > ::tree_dups[treesin[i]][j]){ // if the dup is less than the tree we are on the we need to do some special ops... we get a new list based on the lowest number and we delete only those that are greater than the current value. 
				//cout << "in the else if" << endl;
				vector<int> tempvect = ::tree_dups[::tree_dups[treesin[i]][j]];
				for (unsigned int k = 0; k < tempvect.size(); k++){
					if (treesin[i] < tempvect[k]){ // if the dup is greater in id number than the tree we are on delete it
						uniquetrees.erase(std::remove(uniquetrees.begin(), uniquetrees.end(), tempvect[k]), uniquetrees.end());
					}
				}
			}
    //cout << endl;
		}
	}
	//~ for (unsigned int i = 0; i < uniquetrees.size(); i++){
		//~ cout << uniquetrees[i];
	//~ }
	//~ cout << endl;
	return uniquetrees;
}
*/

vector<string> split(const char *str, char c = ' '){
    vector<string> result;

    while(1){
        const char *begin = str;

        while(*str != c && *str){
			str++;
		}
        result.push_back(string(begin, str));
        if(0 == *str++)
            break;
    }
    return result;
}


int pANTLR3_COMMON_TOKEN_to_int(pANTLR3_COMMON_TOKEN tok){
//	  pANTLR3_COMMON_TOKEN token = $i2;
      ANTLR3_MARKER start = tok->getStartIndex(tok);
      ANTLR3_MARKER end = tok->getStopIndex(tok);
      std::string id((const char *)start, end-start+1 );
      int val = atoi(id.c_str());
	  return val;
}

string pANTLR3_COMMON_TOKEN_string_lit_to_string(pANTLR3_COMMON_TOKEN tok){

      ANTLR3_MARKER start = tok->getStartIndex(tok);
      ANTLR3_MARKER end = tok->getStopIndex(tok);
      std::string id((const char *)start, (end-start+1) );
      string val = id.c_str();
	
	  val.erase( remove( val.begin(), val.end(), '\"' ), val.end() );

	  //cout << "in the transformation function the id is = " << val << endl;
	  return val;
}

string pANTLR3_COMMON_TOKEN_to_string(pANTLR3_COMMON_TOKEN tok){

      ANTLR3_MARKER start = tok->getStartIndex(tok);
      ANTLR3_MARKER end = tok->getStopIndex(tok);
      std::string id((const char *)start, (end-start+1) );
      string val = id.c_str();
	  //cout << "in the transformation function the id is = " << val << endl;
	  return val;
}

double pANTLR3_COMMON_TOKEN_to_double(pANTLR3_COMMON_TOKEN tok){

      ANTLR3_MARKER start = tok->getStartIndex(tok);
      ANTLR3_MARKER end = tok->getStopIndex(tok);
      std::string id((const char *)start, end-start+1 );
	  double val = atof(id.c_str());
	  return val;
}

vector<int> pANTLR3_COMMON_TOKEN_to_intvect(pANTLR3_COMMON_TOKEN tok){

	vector<int> val;
		
	ANTLR3_MARKER begin = tok->getStartIndex(tok);
	ANTLR3_MARKER end = tok->getStopIndex(tok);
	std::string id((const char *)begin, end-begin+1 );
    
    //cout << "in id" <<endl;
    //cout << id << endl;
    
	vector<string> startandstop = split(id.c_str(), '.'); 
	int start =  atoi(startandstop[0].c_str());
	int stop = atoi(startandstop[2].c_str());

    //cout << "in start and stop" <<endl;	
	//print_vector_of_strings(startandstop);
	
	//cout << "in start" <<endl;
	//cout << start << endl;
	//cout << "in stop" <<endl;
	//cout << stop << endl;
	
	
	for (int i = start; i <= stop; i++)
	{
		//cout << "what's in val " << i << endl;
		val.push_back(i);
	}		  
		
	return val;
}

void printBipartTier(int tier){
	for (unsigned int j = 0; j < ::biparttable.biparttable_size(); j++) { // for each bipartition in the hashtable
		if ( tier == ::biparttable.number_of_ones(j)){
			//cout << "Number of ones = "<< ::biparttable.number_of_ones(j) << endl;
			::biparttable.print_bitstring(j);
		}
		if ( tier == ::biparttable.lm.size() - ::biparttable.number_of_ones(j)){
			::biparttable.print_bitstring(j);			
		}
	}
}

void printTaxaTierTuples(int tier, int taxa){
	for (unsigned int j = 0; j < ::biparttable.biparttable_size(); j++) { // for each bipartition in the hashtable
		if ( tier == ::biparttable.number_of_ones(j) && ::biparttable.is_one(j, taxa)){
			//cout << "Number of ones = "<< ::biparttable.number_of_ones(j) << endl;
			::biparttable.print_bitstring(j);
		}
		if ( tier == ::biparttable.lm.size() - ::biparttable.number_of_ones(j) && ::biparttable.is_zero(j, taxa) ){
			::biparttable.print_bitstring(j);
		}
	}
	
}

int TaxaTierTuplesCount(int tier, int taxa){
	int count = 0;
	for (unsigned int j = 0; j < ::biparttable.biparttable_size(); j++) { // for each bipartition in the hashtable
		if ( tier == ::biparttable.number_of_ones(j) && ::biparttable.is_one(j, taxa)){
			//cout << "Number of ones = "<< ::biparttable.number_of_ones(j) << endl;
			//::biparttable.print_bitstring(j);
			count += 1;
		}
		if ( tier == ::biparttable.lm.size() - ::biparttable.number_of_ones(j) && ::biparttable.is_zero(j, taxa) ){
			count += 1;
			//::biparttable.print_bitstring(j);
		}
	}
	return count;
}

float TaxaTierEntropy(vector<int> counts){
	int total = 0;
	float E = 0;
	for (unsigned int i = 0; i < counts.size(); i++){
		total += counts[i];
	}
	
	for (unsigned int i = 0; i < counts.size(); i++){
		E += ((counts[i] / float(total)) * log(counts[i] / float(total)));
	}
	E *= -1; 
	return E;
}

void rateTierRogueness(int tier){
	for (unsigned int i = 0; i < ::biparttable.lm.size(); i++){
		cout << "taxa " << i << " has  " << TaxaTierTuplesCount(tier, i) << " permutations of tier " << tier << endl;
	}
}

void printBipartition(vector<string> leftside, vector<string> rightside)
{
	cout << "r(";
	for(int i = 0; i < (int)leftside.size(); i++)
	{
		cout << leftside[i] << " ";
	}

	cout << "|";
	for(int i = 0; i < (int)rightside.size(); i++)
	{
		cout << rightside[i] << " ";
	}
	cout << ");" << endl;
}

void printSetOfTrees(set<unsigned int> trees)
{
	set<unsigned int>::iterator it = trees.begin();
	if (it == trees.end())
	{
		cout << "Empty set of trees" << endl;
		return;
	}
	

	for (; it!=trees.end(); it++)
	{
		cout << *it << " ";
	}
}

void printBipartition(vector<unsigned int> leftside, vector<unsigned int> rightside)
{
	cout << "r(";
	for(unsigned int i = 0; i < leftside.size(); i++)
	{
		cout << leftside[i] << " ";
	}

	cout << "|";
	for(unsigned int i = 0; i < rightside.size(); i++)
	{
		cout << rightside[i] << " ";
	}
	cout << ")r" << endl;
}

//testing
void printBipartition(int bipartID){
	cout << "[";
	for (unsigned int i = 0; i < ::biparttable.lm.size(); i++){
		if (::biparttable.is_one(bipartID, i)){
			cout << ::biparttable.lm.name(i) << " ";
		}
	}
	cout << "]" << endl;
}


void print_bitstring(bool * bitstring, unsigned int length)
{
  for (unsigned int i = 0; i < length; i++)
  {
    cout << bitstring[i];
  }
  cout << endl;
}


void print_vector_of_bs(vector< bool * > bitstrings, unsigned int length_of_bitstrings){
  cout << "We show bitstrings:" << endl;
  //keep in mind that hashtable and hash_lengths are global variables
  for (unsigned int i = 0; i < bitstrings.size(); i++){
    for (unsigned int j = 0; j < length_of_bitstrings; j++){
      cout << bitstrings[i][j];
    }
    cout << endl;
  }
}

/*template <class T>
void printVector(vector<T> in){
  for(int i = 0; i < in.size(); i++){
	cout << in.at(i) << endl;
	}  
}
*/
/*
void print_list_bs(vector< bool * > list_bs){
  cout << "We show bitstring reps of the bipartitions:" << endl << endl;
  //keep in mind that hashtable and hash_lengths are global variables
  for (unsigned int i = 0; i < NUMBIPART; i++){
    for (unsigned int j = 0; j < ::biparttable.lm.size(); j++){
      cout << list_bs[i][j];
    }
    cout << endl;
  }
}
*/

//~ void print_hashtable(){
  //~ cout << "We show the hash table below, with the corresponding bitstring reps of the bipartitions:" << endl << endl;
  //~ //keep in mind that hashtable and hash_lengths are global variables
  //~ unsigned int mycurrsize = 0;
  //~ for (unsigned int i = 0; i < ::NUMBIPART; i++){
    //~ for (unsigned int j = 0; j < ::length_of_bitstrings[i]; j++){
      //~ cout << ::list_bs[i][j];
    //~ }
    //~ cout << " --> [ ";
    //~ mycurrsize = ::hash_lengths[i];
    //~ for (unsigned int j = 0; j < mycurrsize; j++){
      //~ cout << ::hashtable[i][j] << " ";
    //~ }
    //~ cout << "]" << endl;
  //~ }
//~ }

void print_set(set<unsigned int> t)
{
  set<unsigned int>::iterator it;
  for (it=t.begin(); it!=t.end(); it++)
  	cout << *it << " ";
  cout << endl;    

}


void print_vector_of_strings(vector< string > bitstrings){
  //keep in mind that hashtable and hash_lengths are global variables
  for (unsigned int i = 0; i < bitstrings.size(); i++)
  {
    cout << bitstrings[i] << endl;
  }
}

/*
//TODO- make using template
* Please uses STL's find. 
bool isInVector(vector<int> toSearch, int x){

	for(int i = 0; i < toSearch.size(); i++)
		{
		if(toSearch.at(i)==x){
			return true;
			}
		}
	return false;
}*/


unsigned int factorial(int n){
	return (n<=1) ? 1 : (n * factorial(n-1));
	}
	
void printSetCompactTwo(set<unsigned int> in){
  for(set<unsigned int>::iterator it = in.begin(); it!=in.end(); it++){
	cout << *it << ", ";}
	cout << endl;
}


