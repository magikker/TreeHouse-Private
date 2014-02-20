//file: TreeParsing.cc
#include "TreeParsing.h"
#include <boost/lexical_cast.hpp>

void populate_integrals(unsigned int * hold_integrals, string branches_int, unsigned int encode_size){
  unsigned int integral_size = branches_int.size();
  unsigned int intpos = 0;
  unsigned int treeloc, m;
  int treeval, n;
  unsigned char v;
  treeloc = 0;
  treeval = 0;
  m = 0;
  for (unsigned int i  =0; i < ::biparttable.NumTrees; i++)
    hold_integrals[i] = 0;
  while (intpos < integral_size){
    while (m < encode_size){ //while we're less than the encoded size
      v = branches_int[intpos];
      treeloc+= v;
      treeloc-=33;
      m++;
      if (m < encode_size)
	treeloc *= 100;
      intpos++;
    }
    treeval = branches_int[intpos]; //next, get the tree value
    treeval-=48;
    intpos++;
    if (treeval < 0){
      n = 0;
      treeval = 0;
      while ( (n <= 9) && (n >= 0)){
	treeval*=10;
	n = branches_int[intpos];
	n-=48;
	treeval+= n;
	intpos++;
	n = branches_int[intpos];
	n-=48;
      }
      intpos++;
    }
    hold_integrals[treeloc] = treeval;
    m = 0;
    treeloc = 0;
    treeval = 0;
  }
}

void decompress_branch(unsigned int * hold_integrals, vector<unsigned int> my_set_of_ids, Bipartition &B, string branches_frac){
  unsigned int myplace = 0;
  unsigned int numer = 0;
  unsigned int denom = 100;
  float weight = 0;
  unsigned char v;
  unsigned int val;
  unsigned int count = my_set_of_ids.size();
  if (count < ::biparttable.NumTrees){
    for (unsigned int i = 0; i < count; ++i){
      unsigned int temp = my_set_of_ids[i];
      while (numer != 3){
	    v = branches_frac[myplace];
	    val = v;
	    val -= 33;
	    weight += (float)val/denom;
	    denom*=100;
	    myplace++;
	    numer++;
      }
      weight+=hold_integrals[temp];
      //printf("%f\n", weight);
      hold_integrals[temp] = 0;
      //bitstrings_branches.push_back(weight);
      B.add_tree(temp, weight);
      ::biparttable.tree_branches[temp].push_back(weight);
      numer = 0;
      denom = 100;
      weight = 0;
    }
  }
  else{
    for (unsigned int i = 0; i < ::biparttable.NumTrees; ++i){
      while (numer != 3){
	    v = branches_frac[myplace];
	    val = v;
	    val -= 33;
	    weight += (float)val/denom;
	    denom*=100;
	    myplace++;
	    numer++;
      }
      weight+=hold_integrals[i];
      //printf("%f\n", weight);
      hold_integrals[i] = 0;
      //bitstrings_branches.push_back(weight);
      B.add_tree(i, weight);
      ::biparttable.tree_branches[i].push_back(weight);
      numer = 0;
      denom = 100;
      weight = 0;
    }
  }
}


unsigned int get_bitstring_length(string bitstring){ //determines the size of the bitstring 
  short bitstring_size = bitstring.size();
  int x = 0;
  int ypos = 0;
  int y = 0;
  unsigned int count = 0;
  int val;
  while ( x < bitstring_size ) {
    char type = bitstring[x];
    if ((type == 'A') || (type == 'B')) //if it is either A or B, just increment by one
      count++;
    else if ((type == 'K') || (type == 'L')){ //
      ypos = x + 1;
      val = -1;
      y = 0;
      while (val < 0) { 
        val = bitstring[ypos];
        val-=65;
        if (val < 0) { 
          ++ypos;
          ++y;
        }
        else{
          ypos--;
        }
        if (ypos == bitstring_size)
          break;
      }
      assert(ypos > -1);
      
      string amt = bitstring.substr((x+1),y);
      x = ypos;
      int iamt = atoi(amt.c_str());
      count += iamt;
    }
    x++;
  } //end go through bitstring
  return count;
} 

//parses the ntrees line of the TRZ file
unsigned int get_ntrees(string str){
  unsigned int trees = 0;
  int pos = str.find_first_of(" ");
  string line_type = str.substr(0, pos);
  if (line_type != "NTREES"){
    cerr << "Error! Number of trees field not found. Exiting..\n";
    exit(2);
  }
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  pos = str.find_first_of(" ");
  if (pos != -1){
    string ntrees = str.substr(0, pos);
    str = str.substr(pos+1);
    if ( (str == "H") && (!::biparttable.hetero) ){
      cerr << "Warning! Detected Heterogeneous collection of trees!" << endl;
      ::biparttable.hetero = true;
    }
    trees = atoi(ntrees.c_str());
  }
  else
    trees = atoi(str.c_str()); //this is the total number of trees
  return trees;
} 

//parses the unique line of the TRZ file
unsigned int get_unique(string str, unsigned int ntrees){
  int pos = str.find_first_of(" ");
  string line_type = str.substr(0, pos);
  unsigned int num_unique = ntrees;
  assert(num_unique != 0);
  if (line_type != "UNIQUE_T"){
    cerr << "Error! Unique trees line not found! Exiting...\n";
    exit(2);
  }
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  if (str != "all")
    num_unique = atoi(str.c_str());
  return num_unique;
}

//used for parsing nbipartitions and duplicates lines
void parse_and_get(string str, string check, unsigned int & var){
  int pos;
  pos = str.find_first_of(" ");
  string type =str.substr(0, pos);
  if (type != check){
    cerr << "Error! " << check << " field not found!" << endl;
    cout << str << endl;
    exit(1);
  }
  str = str.substr(pos+1);
  pos = str.find_first_of("\n");
  str = str.substr(0, pos);
  var = atoi(str.c_str());
}


unsigned int decode(string encoded, unsigned int * found){
  //transform string into another string that is "intermediary"
  //cout << "string to decode is: " << encoded << endl;
  unsigned int ids_length = encoded.size();
  char *tmp = new char[std::numeric_limits <int>::digits10 + 2];
  string myids;

  for (unsigned int a = 0; a < encoded.size()-1; ++a) {
    int trueval = encoded[a];
    int trueval2 = encoded[a+1];
    trueval -= 65;
    trueval2 -= 65;
    if (trueval >= 0) { //we have a number!
      if (trueval < 10) { 
	sprintf(tmp, "%d", trueval);
	myids+=tmp;         
      }
      else{
	myids+="1:";
      }
    }
    else{
      myids+=encoded[a];
    }
    if (trueval2 >=0){ 
      myids+=" ";
    }
  }
  int trueval = encoded[ids_length-1];
  trueval-=65;
  if (trueval >= 0){
    sprintf(tmp, "%d", trueval);
    myids+=tmp;
  }
  else{
    myids+=encoded[ids_length-1];
  }
  delete[] tmp;

  //decode the string my ids:
  stringstream splitarray(myids);
  int pos;

  unsigned place = 0;
  string myid;
  unsigned int current = 0;
  while(splitarray >> myid) { 
    pos = myid.find_first_of(":");
    if (pos == -1){ 
      unsigned int temp = atoi(myid.c_str()); 
      temp += current;            
      found[place] = temp;
      place++;
      current = temp;
    }
    else{
      string id_str = myid.substr(0,pos);
      unsigned int id = atoi(id_str.c_str());
      myid = myid.substr(pos+1);
      unsigned int times = atoi(myid.c_str());
      for (unsigned int i = 0; i < times; ++i) {
	found[place] = current+id;
	current+=id;
	place++;
      }
    }
  }
  return place;
}


void decode_bitstring(string bitstring, boost::dynamic_bitset<> &bs, unsigned int maxLength){
  short bitstring_size = bitstring.size();
  unsigned int j = 0;
  int x = 0;
  bool typeb = false;
  int ypos = 0;
  int y = 0;
  while ( x < bitstring_size ) {
    int val = bitstring[x];
    val-= 65;
    if (val < 2) { 
      bs[j] = (bool)val;
      ++j;
    }
    else{
      char type = bitstring[x];
      if (type == 'K')
	typeb = true;
      if (type == 'L')
	typeb = false;
      
      ypos = x + 1;
      val = -1;
      y = 0;
      
      while (val < 0) { 
	val = bitstring[ypos];
	val-=65;
	if (val < 0) { 
	  ++ypos;
	  ++y;
	}
	else{
	  ypos--;
	}
	if (ypos == bitstring_size)
	  break;
      }
      assert(ypos > -1);
      
      string amt = bitstring.substr((x+1),y);
      x = ypos;
      int iamt = atoi(amt.c_str());
      for (short z = 0; z < iamt; ++z) { 
	bs[j] = typeb;
	++j;
      }
    }
    ++x;
  }
  assert(j == maxLength);
} 




string compute_tree_bl( LabelMap lm, vector< boost::dynamic_bitset<> > my_bitstrings, vector< float > my_branches) {
	vector<vector<SCNode*> > vvec_distinctClusters2;

	//update distinct clusters
	for (unsigned int i = 0; i < my_bitstrings.size(); ++i) {
		vector<SCNode*> vec_nodes2;
		for (unsigned int j = 0; j < my_bitstrings[i].size(); j++) {
			if (my_bitstrings[i][j]) {
				SCNode* aNode = new SCNode();
				aNode->name = lm.name(j);
				vec_nodes2.push_back(aNode);
			}
		}
		vvec_distinctClusters2.push_back(vec_nodes2);
	}
  

  for (unsigned int npos = 0; npos < vvec_distinctClusters2.size(); ++npos) { 
    for (unsigned int x = 0; x < vvec_distinctClusters2[npos].size(); ++x) { 
      cout << vvec_distinctClusters2[npos][x]->name << " ";
    }
    cout << endl;
  }
  
  //exit(0);
  
  //if (vvec_distinctClusters2.size() != (NUM_TAXA-2)){
  //  fprintf(stderr, "ERROR!: Tree %d has only %d bipartitions!!\n Exiting...\n", id, vvec_distinctClusters2.size());
  //  exit(0);
  //  }
  
  
	SCTree *scTree = new SCTree();
	bool addedRoot = false;
	unsigned int intNodeNum = 0;
	vector<SCNode*> vec_garbageCan; // for collecting pointers
    
	for (unsigned int pos = 0; pos < vvec_distinctClusters2.size(); ++pos) {
		if (!addedRoot) {
			// The first cluster has all the taxa.
			// This constructs a star tree with all the taxa.
			// 1. Dangle all the taxa as root's children by adjusting parent link.
			// 2. Push all the node* in the root's children
			// 3. Push all the nodes in the tree's nodelist.
			// 4. Push all the nodes' parent in the tree's parentlist.
			cout << "are we assuming the first vvec_distinctClusters2 has all the taxa in it?" << endl;
			cout << "vvec_distinctClusters2[pos].size() = "<< vvec_distinctClusters2[pos].size() << endl;
			for (unsigned int i=0; i < vvec_distinctClusters2[pos].size(); ++i) {
				vvec_distinctClusters2[pos][i]->parent = scTree->root;
				scTree->root->children.push_back(vvec_distinctClusters2[pos][i]);
				scTree->nodelist.push_back(vvec_distinctClusters2[pos][i]);
				assert(scTree->nodelist[0]->name == "root");
				scTree->parentlist2.insert(map<string,int>::value_type(vvec_distinctClusters2[pos][i]->name, 0));
			}
			addedRoot = true;
		}
		else {
			// For the next biggest cluster,
			// 1. Find node list to move (= vvec_distinctClusters2[itr->second]) and
			//    Get the parent node of the to-go nodes.
			// 2. Make an internal node.
			// 3. Insert the node in the nodelist of the tree and update the parentlist accordingly
			// 4. Adjust the to-go nodes' parent link.
			// 5. Adjust the parent node's link to children (delete the moved nodes from children).
      
			// 1. --------------------------------------------------------------------------
			SCNode* theParent = NULL;
			theParent = scTree->nodelist[scTree->parentlist2[vvec_distinctClusters2[pos][0]->name]];
			assert(theParent != NULL);
			assert(theParent->name != "");
      
			// 2. --------------------------------------------------------------------------
			//		string newIntNodeName = "int" + itostr(intNodeNum, 10);
			
			string newIntNodeName = "int" + std::to_string(intNodeNum);
			SCNode* newIntNode = new SCNode();
			vec_garbageCan.push_back(newIntNode);
			newIntNode->name = newIntNodeName;
			newIntNode->parent = theParent;
			
			newIntNode->bl = my_branches[pos];
			
			// 3. --------------------------------------------------------------------------
			assert(newIntNodeName.size() != 0);
			scTree->nodelist.push_back(newIntNode);
			assert(scTree->nodelist[scTree->nodelist.size()-1]->name == newIntNode->name);
      
			scTree->parentlist2.insert(map<string, unsigned>::value_type(newIntNodeName, scTree->nodelist.size()-1));
      
			for (unsigned int i=0; i<vvec_distinctClusters2[pos].size(); ++i) {
				// 4. --------------------------------------------------------------------------
				vvec_distinctClusters2[pos][i]->parent = newIntNode;
	
				// We have to update parentlist in the tree.
				assert(vvec_distinctClusters2[pos][i]->parent->name == scTree->nodelist[scTree->nodelist.size()-1]->name);
	
				scTree->parentlist2[vvec_distinctClusters2[pos][i]->name] = scTree->nodelist.size()-1;
				newIntNode->children.push_back(vvec_distinctClusters2[pos][i]);
	
				// 5. --------------------------------------------------------------------------
				// Delete the moved nodes from parent's children.
				vector<SCNode*>::iterator itr2;
	
				for (itr2 = theParent->children.begin(); itr2 != theParent->children.end(); ++itr2) {
					if (vvec_distinctClusters2[pos][i]->name == (*itr2)->name) {
						theParent->children.erase(itr2);
						break;
					}
				}
			}
			theParent->children.push_back(newIntNode);
			intNodeNum++;
		}
	}
  
	//GET STRING
  
	string mytree;

	mytree = scTree->GetTreeString(true);

	
	//CLEAN UP
	
	//scTree->DrawOnTerminal(true);
   
	for (unsigned int i = 0; i < vvec_distinctClusters2.size(); ++i) {
		for (unsigned int j=0; j<vvec_distinctClusters2[i].size(); ++j) {
			SCNode * tmp = vvec_distinctClusters2[i][j];
			delete tmp;
			vvec_distinctClusters2[i][j] = NULL;
		}
	}

	for (unsigned i=0; i<vec_garbageCan.size(); ++i) {
		if (vec_garbageCan[i]) {
			SCNode* temp = vec_garbageCan[i];
			delete temp;
			vec_garbageCan[i] = NULL;
		}
	}
  
	if (scTree->root) {
		SCNode * tmp = scTree->root;
		delete tmp;
		scTree->root = NULL;
	}
	
	scTree->nodelist.clear();
	scTree->parentlist2.clear();
	delete scTree;

	return mytree;
		
}












string compute_tree( LabelMap lm, vector< boost::dynamic_bitset<> > my_bitstrings, vector< float > my_branches, unsigned id, bool branch) {
	vector<vector<SCNode*> > vvec_distinctClusters2;

	boost::dynamic_bitset<> alltaxa(::biparttable.lm.size());
	
	for(size_t i = 0; i < my_bitstrings.size(); i++){
		alltaxa |= my_bitstrings[i]; 
	}
	my_bitstrings.insert(my_bitstrings.begin(),alltaxa);
	
	//update distinct clusters
	for (unsigned int i = 0; i < my_bitstrings.size(); ++i) {
		vector<SCNode*> vec_nodes2;
		
		//Trival clades don't seem to matter
		if (my_bitstrings[i].count() == 1){
			continue;
		}
		
		//I need the bit strings to come in sorted by size biggest first. The first bit string needs to set the clades for the whole tree. 
		
		for (unsigned int j = 0; j < my_bitstrings[i].size(); j++) {
			if (my_bitstrings[i][j]) {
				SCNode* aNode = new SCNode();
				aNode->name = lm.name(j);
				vec_nodes2.push_back(aNode);
			}
		}
		vvec_distinctClusters2.push_back(vec_nodes2);
	}
  

  for (unsigned int npos = 0; npos < vvec_distinctClusters2.size(); ++npos) { 
    for (unsigned int x = 0; x < vvec_distinctClusters2[npos].size(); ++x) { 
      cout << vvec_distinctClusters2[npos][x]->name << " ";
    }
    cout << endl;
  }
  
  //exit(0);
  //if (vvec_distinctClusters2.size() != (NUM_TAXA-2)){
  //  fprintf(stderr, "ERROR!: Tree %d has only %d bipartitions!!\n Exiting...\n", id, vvec_distinctClusters2.size());
  //  exit(0);
  //  }
  
  
	SCTree *scTree = new SCTree();
	bool addedRoot = false;
	unsigned int intNodeNum = 0;
	vector<SCNode*> vec_garbageCan; // for collecting pointers
    
	for (unsigned int pos = 0; pos < vvec_distinctClusters2.size(); ++pos) {
		if (!addedRoot) {
			// The first cluster has all the taxa.
			// This constructs a star tree with all the taxa.
			// 1. Dangle all the taxa as root's children by adjusting parent link.
			// 2. Push all the node* in the root's children
			// 3. Push all the nodes in the tree's nodelist.
			// 4. Push all the nodes' parent in the tree's parentlist.
			cout << "are we assuming the first vvec_distinctClusters2 has all the taxa in it?" << endl;
			cout << "vvec_distinctClusters2[pos].size() = "<< vvec_distinctClusters2[pos].size() << endl;

			for (unsigned int i=0; i < vvec_distinctClusters2[pos].size(); ++i) {
				vvec_distinctClusters2[pos][i]->parent = scTree->root;
				scTree->root->children.push_back(vvec_distinctClusters2[pos][i]);
				scTree->nodelist.push_back(vvec_distinctClusters2[pos][i]);
				assert(scTree->nodelist[0]->name == "root");
				scTree->parentlist2.insert(map<string,int>::value_type(vvec_distinctClusters2[pos][i]->name, 0));
			}
			addedRoot = true;
		}
		else {
			// For the next biggest cluster,
			// 1. Find node list to move (= vvec_distinctClusters2[itr->second]) and
			//    Get the parent node of the to-go nodes.
			// 2. Make an internal node.
			// 3. Insert the node in the nodelist of the tree and update the parentlist accordingly
			// 4. Adjust the to-go nodes' parent link.
			// 5. Adjust the parent node's link to children (delete the moved nodes from children).
      
			// 1. --------------------------------------------------------------------------
			SCNode* theParent = NULL;
			theParent = scTree->nodelist[scTree->parentlist2[vvec_distinctClusters2[pos][0]->name]];
			assert(theParent != NULL);
			assert(theParent->name != "");
      
			// 2. --------------------------------------------------------------------------
			//		string newIntNodeName = "int" + itostr(intNodeNum, 10);
			
			string newIntNodeName = "int" + std::to_string(intNodeNum);
			SCNode* newIntNode = new SCNode();
			vec_garbageCan.push_back(newIntNode);
			newIntNode->name = newIntNodeName;
			newIntNode->parent = theParent;
			if (::biparttable.weighted){
				newIntNode->bl = my_branches[pos];
			}
			// 3. --------------------------------------------------------------------------
			assert(newIntNodeName.size() != 0);
			scTree->nodelist.push_back(newIntNode);
			assert(scTree->nodelist[scTree->nodelist.size()-1]->name == newIntNode->name);
      
			scTree->parentlist2.insert(map<string, unsigned>::value_type(newIntNodeName, scTree->nodelist.size()-1));
      
			for (unsigned int i=0; i<vvec_distinctClusters2[pos].size(); ++i) {
				// 4. --------------------------------------------------------------------------
				vvec_distinctClusters2[pos][i]->parent = newIntNode;
	
				// We have to update parentlist in the tree.
				assert(vvec_distinctClusters2[pos][i]->parent->name == scTree->nodelist[scTree->nodelist.size()-1]->name);
	
				scTree->parentlist2[vvec_distinctClusters2[pos][i]->name] = scTree->nodelist.size()-1;
				newIntNode->children.push_back(vvec_distinctClusters2[pos][i]);
	
				// 5. --------------------------------------------------------------------------
				// Delete the moved nodes from parent's children.
				vector<SCNode*>::iterator itr2;
	
				for (itr2 = theParent->children.begin(); itr2 != theParent->children.end(); ++itr2) {
					if (vvec_distinctClusters2[pos][i]->name == (*itr2)->name) {
						theParent->children.erase(itr2);
						break;
					}
				}
			}
			theParent->children.push_back(newIntNode);
			intNodeNum++;
		}
	}
  
	//GET STRING
  
	string mytree;

	if (branch && ::biparttable.weighted){
		mytree = scTree->GetTreeString(true);
	}
	else{
		mytree = scTree->GetTreeString(false);
	}
	
	//CLEAN UP
	
	//scTree->DrawOnTerminal(true);
   
	for (unsigned int i = 0; i < vvec_distinctClusters2.size(); ++i) {
		for (unsigned int j=0; j<vvec_distinctClusters2[i].size(); ++j) {
			SCNode * tmp = vvec_distinctClusters2[i][j];
			delete tmp;
			vvec_distinctClusters2[i][j] = NULL;
		}
	}

	for (unsigned i=0; i<vec_garbageCan.size(); ++i) {
		if (vec_garbageCan[i]) {
			SCNode* temp = vec_garbageCan[i];
			delete temp;
			vec_garbageCan[i] = NULL;
		}
	}
  
	if (scTree->root) {
		SCNode * tmp = scTree->root;
		delete tmp;
		scTree->root = NULL;
	}
	
	scTree->nodelist.clear();
	scTree->parentlist2.clear();
	delete scTree;

	return mytree;
		
}

void load_data_from_trz_file(string file, BipartitionTable &Tab){  

	string mycount, ids, bipartition, line_type, str, taxa, treeline, bitstring, branches_int, branches_frac;
	unsigned int bipart_count, nbipart, ndup_lines, NumTrees;
	vector<bool> check;
	vector<unsigned int> true_ids;
	vector <bool> is_dup;
	bipart_count = 0;
	unsigned int num_unique;
 
	//GRB
	//cout << std::numeric_limits<unsigned long>::digits << endl;
 
 
	ifstream fin(file.c_str(), ios::binary);
	if (!fin) {
		cerr << "cannot open file!\n";
		exit(2);
	}
  
	//read in taxa
	getline(fin, str);  
	int pos = str.find_first_of(" ");
	line_type = str.substr(0, pos);
	if (line_type != "TAXA"){
		cerr << "Error! No taxa labels identified for file! Exiting..\n";
		exit(2);
	}
	str = str.substr(pos+1);
	pos = str.find_first_of("\n");
	str = str.substr(0, pos);
	//fill LM
	stringstream ss(str);
	while(getline(ss, taxa, ':')){ 
		Tab.lm.push(taxa);
	}
	ss.clear();

	//read in number of trees
	getline(fin, str);
	NumTrees = get_ntrees(str);
	Tab.NumTrees = NumTrees;
	//GRB
	for (unsigned int i = 0; i < NumTrees; i++){
		boost::dynamic_bitset<> TNT(Tab.lm.size());
		Tab.taxa_in_trees.push_back(TNT);
	}
	//cout << "table created but not filled. " << endl;
	
	//read in number of unique trees
	getline(fin, str);
	num_unique = get_unique(str, NumTrees);

	//resize data structures
	Tab.tree_branches.resize(NumTrees);
	
	Tab.inverted_index.resize(NumTrees);

	check.resize(NumTrees);
	true_ids.resize(num_unique);
	is_dup.resize(NumTrees);
	unsigned int * hold_integrals = NULL;

	//read in number of bipartitions to be read
	getline(fin, str);
	parse_and_get(str, "NBIPARTITIONS", nbipart);
	//NUMBIPART = nbipart;
	//cout << "NUMBIPART = " << NUMBIPART << endl;
	//for (unsigned int i = 0; i < NUM_TREES; i++){
	//	bool *blankbs = new bool[::NUM_TAXA];
	//	for (unsigned int j = 0; j < ::NUM_TAXA; j++){
	//		blankbs[j] = (bool)0;
	//	}
	//	::biparttable.taxa_in_trees.push_back(blankbs);
	//}

  //skip nbipart lines, get the duplicate information
  for (unsigned int i = 0; i < nbipart; i++){
    getline(fin, str);
  }
  getline(fin, str); //now read the duplicates line
  parse_and_get(str, "DUPLICATES", ndup_lines);
  //cout << "ndup_lines is: " << ndup_lines << endl;
  for (unsigned int i = 0; i < NumTrees; i++){
    is_dup[i] = 0;
  }
  unsigned int decode_size = 0;
  unsigned int * found, dec_loc, dec_val;
  found = (unsigned int *)malloc(NumTrees*sizeof(unsigned int));
  for (unsigned int i = 0; i < ndup_lines; i++){
    getline(fin, str); //read in the line
    pos = str.find_first_of("\n");
    str = str.substr(0, pos); //get rid of the newline
    decode_size = decode(str, found); //decode the line
    dec_loc = found[0];
    for (unsigned int j = 1; j < decode_size; j++){
      dec_val = found[j];
      Tab.tree_dups[dec_loc].push_back(dec_val);
      //dups[dec_loc].push_back(dec_val); //add the elements of the array to the associated dup structure
      is_dup[dec_val] = 1; //also set those locations in our bool structure to a 1
    }
  }

  //populate true_ids... The order the trees were in in the input file. (first tree = tree 1) and so on. 
  unsigned int tempj = 0;
  for (unsigned int i = 0; i < num_unique; i++){
    if (is_dup[tempj]){
      while (is_dup[tempj])
	  tempj++;
    }
    true_ids[i] = tempj;
    tempj++;
  }
 
  fin.close(); //now close the file
  //GRB there should be a function break here.... This is far too big. 
  fin.open(file.c_str(), ios::binary); //reopen
  if (!fin) {
    cerr << "cannot open file!\n";
    exit(2);
  }
  for (unsigned int i  =0; i < 4; i++) {//go back to where we need to be in order to start reading in the bipartitions
    getline(fin, str); 
  }
  
  unsigned int counter = 0;
  unsigned int count = 0;
  unsigned int bipart_loc = 0;
  unsigned int numtrees_strsize, my_count;
  string numtrees_str = std::to_string(NumTrees);
  unsigned int encode_size = 0;
  vector<unsigned int> my_set_of_ids;
  numtrees_strsize = numtrees_str.size();
  encode_size = (numtrees_strsize+1)/2;
  unsigned int maxLength = 0;


  while ( counter < nbipart) { 
	//GRB NEW Round
    boost::dynamic_bitset<> tt1(NumTrees); //tt
	Tab.treetable.push_back(tt1); //tt
	
    getline(fin, str);  
    pos = str.find_first_of(" ");
    bitstring = str.substr(0, pos); //contains bitstring
    str = str.substr(pos+1);
    pos = str.find_first_of("\n");
    treeline = str.substr(0, pos);

    //process bitstring first
    // MAYBE REMOVE
    maxLength = get_bitstring_length(bitstring);//determine length of bitstring maxLength
    //GRB the bitset
    //boost::dynamic_bitset<> bs(maxLength);
    boost::dynamic_bitset<> bs(Tab.lm.size());
    
    
    //bool *bs = new bool[maxLength]; //allocate it to be maxLength	
    decode_bitstring(bitstring, bs, maxLength);
    //TAKEOUT
    Bipartition B(bs);

    
    //next, process tree line
    //first, determine the number of TIDs in the line:
    pos = treeline.find_first_of(":");
    line_type = treeline.substr(0, pos);
    treeline = treeline.substr(pos+1);
    pos = treeline.find_first_of(":");
    mycount = treeline.substr(0, pos);
    count = atoi(mycount.c_str());
    treeline = treeline.substr(pos+1);
    pos = treeline.find_first_of(" ");
    if (pos != -1){
      if (!Tab.weighted){
		fprintf(stderr, "\nFile has branch lengths! Weighted option: ON\n\n");	
		Tab.weighted = true;
		hold_integrals = (unsigned int*)malloc(NumTrees*sizeof(unsigned int));
		for (unsigned int x = 0; x < NumTrees; x++){
			hold_integrals[x] = 0;
		}
      }
    }
    if (Tab.weighted){
      ids = treeline.substr(0, pos);
      treeline = treeline.substr(pos+1);
      pos = treeline.find_first_of(" ");
      branches_int = treeline.substr(0, pos); //integral portion of branches
      treeline = treeline.substr(pos+1);
      pos = treeline.find_first_of("\n");
      branches_frac = treeline.substr(0, pos); //fractional component of branches
    }
    else{
      pos = treeline.find_first_of("\n");
      ids = treeline.substr(0, pos);
    }
	
	//THENEW
	set<unsigned int> TreeIDs;
	vector<float> branchLengths;
	bool Inverse = false;
	bool NoBranch = true;
	
    //process tree ids
    if (count == 0){  //bipartition that every tree has
	Inverse = true;
     for (unsigned int b = 0; b < NumTrees; ++b){ 
		Tab.inverted_index[b].push_back(bipart_loc);
		for(unsigned int eachbit = 0; eachbit < bs.size(); eachbit++){
			if (bs[eachbit]){
				Tab.taxa_in_trees[b].set(eachbit);
			}
		}
		
		//GRB
		if(!Tab.weighted){
		   B.add_tree(b,1.0);
	    }
		//B.add_tree(b,1.0);
		//::biparttable.searchtable[counter].push_back(b);
	    Tab.treetable[counter][b] = 1; //tt
	    
      }
      if (Tab.weighted)
	  {
		if (branches_int != ""){
	      populate_integrals(hold_integrals, branches_int, encode_size);
		}
	    my_set_of_ids.resize(NumTrees);
	    //GRB
	    decompress_branch(hold_integrals, my_set_of_ids, B, branches_frac);
      }
      bipart_count++;      
    }
    
         
    else { //count != 0 (so we have tree ids to process)
      my_count = decode(ids, found);
      assert(my_count == count);

    if (line_type == "-") { //compressed line
    Inverse = true;
		for (unsigned int i = 0; i < count; ++i) {
			unsigned int temp = found[i];
			TreeIDs.insert(temp);
			check[temp] = true;
		}
		unsigned int true_id, sec_id;
		//cout << "we do get here though? num_unique is: " << num_unique << endl;
		for (unsigned int i = 0; i < num_unique; ++i) {
			if (check[i] == false) {
				true_id = true_ids[i];
				my_set_of_ids.push_back(true_id);
				Tab.inverted_index[true_id].push_back(bipart_loc);

				if (Tab.tree_dups[true_id].size() > 0){
				  for (unsigned int j = 0; j< Tab.tree_dups[true_id].size(); j++){
						sec_id = Tab.tree_dups[true_id][j];
						Tab.inverted_index[sec_id].push_back(bipart_loc);
						my_set_of_ids.push_back(sec_id);
					}
				}
				//cout << "end of for stat" << endl;
			}
			else
				check[i] = false;
		}
		//now, take care of of the bipartition associated with this
		sort(my_set_of_ids.begin(), my_set_of_ids.end());
		unsigned int mytotalsize = my_set_of_ids.size();
		for (unsigned int j = 0; j < mytotalsize; j++){
			//GRB
			if(!Tab.weighted){
				B.add_tree(my_set_of_ids[j],1.0);
			}
			for(unsigned int eachbit = 0; eachbit < bs.size(); eachbit++){
				if (bs[eachbit]){
					Tab.taxa_in_trees[my_set_of_ids[j]].set(eachbit);
				}
			}
			//B.add_tree(my_set_of_ids[j],1.0);
			//trees.push_back(my_set_of_ids[j]);
			//::biparttable.searchtable[counter].push_back(my_set_of_ids[j]);
	        Tab.treetable[counter][my_set_of_ids[j]] = 1; //tt
		}
		if (Tab.weighted){	  
			if (branches_int != ""){
				populate_integrals(hold_integrals, branches_int, encode_size);
			}
			//GRB
			decompress_branch(hold_integrals, my_set_of_ids, B, branches_frac);			
		}
      }
      
      
      
      else { //line is not compressed
      Inverse = false;
		unsigned int true_id, sec_id;
		for (unsigned int i = 0; i < count; ++i) { 
			unsigned int temp = found[i];
			true_id = true_ids[temp];
			Tab.inverted_index[true_id].push_back(bipart_loc);
			my_set_of_ids.push_back(true_id);
			TreeIDs.insert(true_id);
			if (Tab.tree_dups[true_id].size() > 0){
				for (unsigned int j = 0; j < Tab.tree_dups[true_id].size(); j++){
					sec_id = Tab.tree_dups[true_id][j];
					Tab.inverted_index[sec_id].push_back(bipart_loc);
					my_set_of_ids.push_back(sec_id);
					TreeIDs.insert(sec_id);
				}
			}
		}	 
		sort(my_set_of_ids.begin(), my_set_of_ids.end());	 
		unsigned int mytotalsize = my_set_of_ids.size();
		for (unsigned int j = 0; j < mytotalsize; j++){
			//GRB
			if(!Tab.weighted){
				B.add_tree(my_set_of_ids[j],1.0);
			}
			for(unsigned int eachbit = 0; eachbit < bs.size(); eachbit++){
				if (bs[eachbit]){
					Tab.taxa_in_trees[my_set_of_ids[j]].set(eachbit);
				}
			}
			//B.add_tree(my_set_of_ids[j], 1.0);
			//trees.push_back(my_set_of_ids[j]);
			//::biparttable.searchtable[counter].push_back(my_set_of_ids[j]);
	        Tab.treetable[counter][my_set_of_ids[j]] = 1; //tt
		}
		if (Tab.weighted){
			if (branches_int != ""){
				populate_integrals(hold_integrals, branches_int, encode_size);
			}
			//GRB
			decompress_branch(hold_integrals, my_set_of_ids, B, branches_frac);
		}
      }
      bipart_count++;
    } //end if count != 0
    
    //GRB
    //cout << "lets see what we've got" << endl;
    //cout << "trees = ";
    //for(vector<unsigned int>::iterator it = trees.begin(); it != trees.end(); it++){
	//	cout << *it << " "; 
	//}
	//cout << endl;
    
    //THENEW
    
    if(!Tab.weighted){
		Tab.CladeMap.insert(std::make_pair(bs, TreeSet (TreeIDs, Inverse)));
	}
	else{
		Tab.CladeMap.insert(std::make_pair(bs, TreeSet (TreeIDs, Inverse, branchLengths )));
	}
    //m_mapFoo->insert(std::make_pair(0, Foo(0)));
    Tab.BipartTable.push_back(B);
    //cout << "This many times through" << counter << endl;
    //::biparttable.bipartitions.push_back(bs);
    //cout << sizeof(biparttable) << " " << sizeof(biparttable.BipartitionTable) << endl;
    my_set_of_ids.clear();
    bipart_loc++;
    counter++;
  }
  
  
  
  //set the benchmarks
  Tab.MapBenchMarks.push_back(Tab.CladeMap.begin());
  for (unsigned int i = 1; i < Tab.lm.size()+1; i++){
	Tab.MapBenchMarks.push_back(Tab.CladeMap.end());
  }
  
  unsigned int current = 0;
  for (std::map<boost::dynamic_bitset<>, TreeSet >::iterator it=Tab.CladeMap.begin(); it!=Tab.CladeMap.end(); ++it){
	
	while(it->first.count() > current){
		current++;
		Tab.MapBenchMarks[current] = it;
	}
	
	//if(it->first.count() > current){
	//cout << "it->first.count() = "<<it->first.count() << "current = " << current << endl;
	//}
  }
  
  
  for (unsigned int i = Tab.MapBenchMarks.size(); i < Tab.lm.size() +1 ; i++){
	 Tab.MapBenchMarks.push_back(Tab.CladeMap.end());
  }
  Tab.MapBenchMarks.push_back(Tab.CladeMap.end());
  
  //cout << "benchmark size = "<<Tab.MapBenchMarks.size() << endl;
  
  
  //now that all the trees are loaded, calculate  all trivial bipartitions
  Tab.calculate_trivial_bipartitions();
  cout << "loaded all the trees" << endl;
  free(hold_integrals);
  free(found);
  if (bipart_count != counter) { 
    cerr << "ERROR! Bipartitions not processed correctly!" << endl;
    cout << "bipart_count: " << bipart_count << endl;
    cout << "counter: " << counter << endl;
    exit(1);
  }
  assert(Tab.lm.size()!=0); 
 
    for (unsigned int i = 0; i < NumTrees; i++){
		for (unsigned int j = 0; j < Tab.tree_dups[i].size(); j++){
			if (Tab.tree_dups[Tab.tree_dups[i][j]].size() == 0){
				Tab.tree_dups[Tab.tree_dups[i][j]].push_back(i);
			}
		}
		//cout << endl;
    }
   
	debugstatement("Hey! you got to the end of the parsing from file function. Nice.");
   
} 

