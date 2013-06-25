//*****************************************************************/
/*
* TreeHouse is an interface for phylogenetic trees. 
* It is based on TreeZip, developed by Suzanne Matthews and 
* HashCS developed by Seung-Jin Sul.
* 
* Contributers to this project include Grant Brammer,
* Arthur Philpott, and Mark Adamo. 
*
* (c) 2012 TreeHouse : Grant Brammer
* (c) 2010 TreeZip: Suzanne Matthews
* (c) 2009 HashRF : Seung-Jin Sul 
* (c) 2009 HashCS : Seung-Jin Sul
*
* This file is part of TreeHouse.
*
* TreeHouse is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* TreeHouse is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with TreeZip.  If not, see <http://www.gnu.org/licenses/>.
*
*/

#include <ctime>
#include <iostream>

#include <readline/readline.h>
#include <readline/history.h>

#include <libxml++/libxml++.h>
#include <libxml++/parsers/textreader.h>

//Global Vars like the hashtable.
//#include "global.h"
#include "THGlobals.h"

//The label map class. 
#include "label-map.hh"

//#include "HashTableSearch.h"

//#include "compressfunc.h"
//#include "bmerge.h"
#include "buildtree.h"

#include "UserFunctions.h"
#include "BipartitionTable.h"

// The lexer/parser
#include    "pql.h"
#include "quartet.h"
#include "distance.h"

using namespace std;

bool interactive = false;
ofstream interactive_log;

static char *line_read = (char *)NULL;

//char* cmd [] ={ "hello", "world", "hell" ,"word", "quit", " " };


void * xmalloc (int size)
{
    void *buf;
 
    buf = malloc (size);
    if (!buf) {
        fprintf (stderr, "Error: Out of memory. Exiting.'n");
        exit (1);
    }
 
    return buf;
}

char * dupstr (char* s) {
  char *r;
 
  r = (char*) xmalloc ((strlen (s) + 1));
  strcpy (r, s);
  return (r);
}


char* my_generator(const char* text, int state){
    static int list_index, len;
    //char *name;
 
    if (!state) {
        list_index = 0;
        len = strlen (text);
    }

	for (unsigned int i = list_index; i < ::functionKeys.size(); i++) {
		list_index++;
		//char *name = new char[*iter->first.length()+1];
		char name[::functionKeys[i].size()+1];
		strcpy(name, ::functionKeys[i].c_str());
		//cout << "name = " << name << endl;
		if(strncmp(name, text, len) == 0)
			return dupstr(name);
    }

    /* If no names matched, then return NULL. */
    return ((char *)NULL);
 
}

static char** my_completion( const char * text , int start,  int end){ //this is where I could say to use each matching function. 
    char **matches;
 
    matches = (char **)NULL;
 
    matches = rl_completion_matches ((char*)text, &my_generator);
    
    return (matches);
}
 


/* Read a string, and return a pointer to it.  Returns NULL on EOF. */
char * do_gets ()
{
  /* If the buffer has already been allocated, return the memory
     to the free pool. */
  if (line_read != (char *)NULL)
    {
      free (line_read);
      line_read = (char *)NULL;
    }

  /* Get a line from the user. */
  line_read = readline ("");

  /* If the line has any text in it, save it on the history. */
  if (line_read && *line_read)
    add_history (line_read);

  return (line_read);
}
/*
void differentLengthsOr(const bool * first, const unsigned int firstsize, bool * second){
	for(unsigned int i = 0; i < firstsize; i++){
		if (first[i] == 1){
			second[i] = 1;
		}
	}
}
*/
void differentLengthsOr(const bool * first, const unsigned int firstsize, bool * second){
	for(unsigned int i = 0; i < firstsize; i++){
		if (first[i] == 1){
			second[i] = 1;
		}
	}
}

void populate_integrals(unsigned int * hold_integrals, string branches_int, unsigned int encode_size){
  unsigned int integral_size = branches_int.size();
  unsigned int intpos = 0;
  unsigned int treeloc, m;
  int treeval, n;
  unsigned char v;
  treeloc = 0;
  treeval = 0;
  m = 0;
  for (unsigned int i  =0; i < NUM_TREES; i++)
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
  if (count < NUM_TREES){
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
    for (unsigned int i = 0; i < NUM_TREES; ++i){
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
    if ( (str == "H") && (!HETERO) ){
      cerr << "Warning! Detected Heterogeneous collection of trees!" << endl;
      HETERO = true;
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

void load_data_from_trz_file(string file){  
  
	string mycount, ids, bipartition, line_type, str, taxa, treeline, bitstring, branches_int, branches_frac;
	unsigned int bipart_count, ntaxa, nbipart, ndup_lines;
	vector<bool> check;
	vector<unsigned int> true_ids;
	vector <bool> is_dup;
	bipart_count = 0;
	unsigned int num_unique;
 
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
	stringstream ss(str);
	ntaxa = 0;
	while(getline(ss, taxa, ':')){ 
		::biparttable.lm.push(taxa);
		ntaxa++;
	}
	::NUM_TAXA = ntaxa;
	ss.clear();

	//read in number of trees
	getline(fin, str);
	::NUM_TREES = get_ntrees(str);
	//GRB
	for (unsigned int i = 0; i < ::NUM_TREES; i++){
		boost::dynamic_bitset<> TNT(::NUM_TAXA);
		biparttable.taxa_in_trees.push_back(TNT);
	}
	//cout << "table created but not filled. " << endl;
	
	//read in number of unique trees
	getline(fin, str);
	num_unique = get_unique(str, ::NUM_TREES);

	//resize data structures
	::biparttable.tree_branches.resize(NUM_TREES);
	
	::inverted_index.resize(NUM_TREES);

	check.resize(NUM_TREES);
	true_ids.resize(num_unique);
	is_dup.resize(NUM_TREES);
	unsigned int * hold_integrals = NULL;

	//read in number of bipartitions to be read
	getline(fin, str);
	parse_and_get(str, "NBIPARTITIONS", nbipart);
	NUMBIPART = nbipart;
	cout << "NUMBIPART = " << NUMBIPART << endl;
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
  for (unsigned int i = 0; i < NUM_TREES; i++){
    is_dup[i] = 0;
  }
  unsigned int decode_size = 0;
  unsigned int * found, dec_loc, dec_val;
  found = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  for (unsigned int i = 0; i < ndup_lines; i++){
    getline(fin, str); //read in the line
    pos = str.find_first_of("\n");
    str = str.substr(0, pos); //get rid of the newline
    decode_size = decode(str, found); //decode the line
    dec_loc = found[0];
    for (unsigned int j = 1; j < decode_size; j++){
      dec_val = found[j];
      ::tree_dups[dec_loc].push_back(dec_val);
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
  string numtrees_str = std::to_string(NUM_TREES);
  unsigned int encode_size = 0;
  vector<unsigned int> my_set_of_ids;
  numtrees_strsize = numtrees_str.size();
  encode_size = (numtrees_strsize+1)/2;
  unsigned int maxLength = 0;


  while ( counter < nbipart) { 
	//GRB NEW Round
    boost::dynamic_bitset<> tt1(NUM_TREES); //tt
	::biparttable.treetable.push_back(tt1); //tt
	
    getline(fin, str);  
    pos = str.find_first_of(" ");
    bitstring = str.substr(0, pos); //contains bitstring
    str = str.substr(pos+1);
    pos = str.find_first_of("\n");
    treeline = str.substr(0, pos);

    //process bitstring first
    maxLength = get_bitstring_length(bitstring);//determine length of bitstring maxLength
    //GRB the bitset
    boost::dynamic_bitset<> bs(maxLength);
    
    //bool *bs = new bool[maxLength]; //allocate it to be maxLength	
    decode_bitstring(bitstring, bs, maxLength);
    Bipartition B(bs);

    
    //GRB not needed in new system
    //::biparttable.length_of_bitstrings.push_back(maxLength);
    
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
      if (!WEIGHTED){
		fprintf(stderr, "\nFile has branch lengths! Weighted option: ON\n\n");	
		WEIGHTED = true;
		hold_integrals = (unsigned int*)malloc(NUM_TREES*sizeof(unsigned int));
		for (unsigned int x = 0; x < NUM_TREES; x++){
			hold_integrals[x] = 0;
		}
      }
    }
    if (WEIGHTED){
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

    //process tree ids
    if (count == 0){  //bipartition that every tree has
     for (unsigned int b = 0; b < ::NUM_TREES; ++b){ 
		::inverted_index[b].push_back(bipart_loc);
		for(unsigned int eachbit = 0; eachbit < bs.size(); eachbit++){
			if (bs[eachbit]){
				biparttable.taxa_in_trees[b].set(eachbit);
			}
		}
		
		//GRB
		if(!WEIGHTED){
		   B.add_tree(b,1.0);
	    }
		//B.add_tree(b,1.0);
		//::biparttable.searchtable[counter].push_back(b);
	    ::biparttable.treetable[counter][b] = 1; //tt
	    
      }
      if (WEIGHTED)
	  {
		if (branches_int != ""){
	      populate_integrals(hold_integrals, branches_int, encode_size);
		}
	    my_set_of_ids.resize(NUM_TREES);
	    //GRB
	    decompress_branch(hold_integrals, my_set_of_ids, B, branches_frac);
      }
      bipart_count++;      
    }      
    else { //count != 0 (so we have tree ids to process)
      my_count = decode(ids, found);
      assert(my_count == count);

    if (line_type == "-") { //compressed line
		for (unsigned int i = 0; i < count; ++i) {
			unsigned int temp = found[i];
			check[temp] = true;
		}
		unsigned int true_id, sec_id;
		//cout << "we do get here though? num_unique is: " << num_unique << endl;
		for (unsigned int i = 0; i < num_unique; ++i) {
			if (check[i] == false) {
				true_id = true_ids[i];
				my_set_of_ids.push_back(true_id);
				::inverted_index[true_id].push_back(bipart_loc);

				if (::tree_dups[true_id].size() > 0){
				  for (unsigned int j = 0; j< ::tree_dups[true_id].size(); j++){
						sec_id = ::tree_dups[true_id][j];
						::inverted_index[sec_id].push_back(bipart_loc);
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
			if(!WEIGHTED){
				B.add_tree(my_set_of_ids[j],1.0);
			}
			for(unsigned int eachbit = 0; eachbit < bs.size(); eachbit++){
				if (bs[eachbit]){
					biparttable.taxa_in_trees[my_set_of_ids[j]].set(eachbit);
				}
			}
			//B.add_tree(my_set_of_ids[j],1.0);
			//trees.push_back(my_set_of_ids[j]);
			//::biparttable.searchtable[counter].push_back(my_set_of_ids[j]);
	        ::biparttable.treetable[counter][my_set_of_ids[j]] = 1; //tt
		}
		if (WEIGHTED){	  
			if (branches_int != ""){
				populate_integrals(hold_integrals, branches_int, encode_size);
			}
			//GRB
			decompress_branch(hold_integrals, my_set_of_ids, B, branches_frac);			
		}
      }
      else { //line is not compressed
		unsigned int true_id, sec_id;
		for (unsigned int i = 0; i < count; ++i) { 
			unsigned int temp = found[i];
			true_id = true_ids[temp];
			::inverted_index[true_id].push_back(bipart_loc);
			my_set_of_ids.push_back(true_id);
			if (::tree_dups[true_id].size() > 0){
				for (unsigned int j = 0; j < ::tree_dups[true_id].size(); j++){
					sec_id = ::tree_dups[true_id][j];
					::inverted_index[sec_id].push_back(bipart_loc);
					my_set_of_ids.push_back(sec_id);
				}
			}
		}	 
		sort(my_set_of_ids.begin(), my_set_of_ids.end());	 
		unsigned int mytotalsize = my_set_of_ids.size();
		for (unsigned int j = 0; j < mytotalsize; j++){
			//GRB
			if(!WEIGHTED){
				B.add_tree(my_set_of_ids[j],1.0);
			}
			for(unsigned int eachbit = 0; eachbit < bs.size(); eachbit++){
				if (bs[eachbit]){
					biparttable.taxa_in_trees[my_set_of_ids[j]].set(eachbit);
				}
			}
			//B.add_tree(my_set_of_ids[j], 1.0);
			//trees.push_back(my_set_of_ids[j]);
			//::biparttable.searchtable[counter].push_back(my_set_of_ids[j]);
	        ::biparttable.treetable[counter][my_set_of_ids[j]] = 1; //tt
		}
		if (WEIGHTED){
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
    
    biparttable.BipartitionTable.push_back(B);
    //cout << "This many times through" << counter << endl;
    //::biparttable.bipartitions.push_back(bs);
    //cout << sizeof(biparttable) << " " << sizeof(biparttable.BipartitionTable) << endl;
    my_set_of_ids.clear();
    bipart_loc++;
    counter++;
  }
  //now that all the trees are loaded, calculate  all trivial bipartitions
  biparttable.calculate_trivial_bipartitions();
  cout << "loaded all the trees" << endl;
  free(hold_integrals);
  free(found);
  if (bipart_count != counter) { 
    cerr << "ERROR! Bipartitions not processed correctly!" << endl;
    cout << "bipart_count: " << bipart_count << endl;
    cout << "counter: " << counter << endl;
    exit(1);
  }
  assert(::biparttable.lm.size()!=0); 
 
    for (unsigned int i = 0; i < ::NUM_TREES; i++){
		for (unsigned int j = 0; j < ::tree_dups[i].size(); j++){
			if (::tree_dups[::tree_dups[i][j]].size() == 0){
				::tree_dups[::tree_dups[i][j]].push_back(i);
			}
		}
		//cout << endl;
    }
   
	debugstatement("Hey! you got to the end of the parsing from file function. Nice.");
   
} 

void load_classification(string filename) {
  xmlpp::TextReader reader(filename);
  reader.read();
  if (reader.get_name() != "taxa") {
    cerr << "Error in XML parse of '" << filename << "'. Expected <taxa> as top-level tag." << endl;
    exit(3);
  }
  reader.read();
  if (reader.get_name() == "#text")
    reader.read();
  while (reader.get_name() != "taxa") {
    if (reader.get_name() != "taxon") {
      cerr << "Error in XML parse of '" << filename << "'. Expected <taxon> statement but found <" << reader.get_name() << ">." << endl;
      exit(3);
    }
    string label = reader.get_attribute("label");
    if (label == "") {
      cerr << "Error in XML parse of '" << filename << "'. One of the given taxa has no label attribute." << endl;
      exit(3);
    }
    int ind = ::biparttable.lm.index_in_labelmap(label);
    if (ind == -1) {
      cerr << "Error in XML parse of '" << filename << "'. Label '" << label << "' does not exist in labelmap." << endl;
      exit(3);
    }
    reader.read();
    if (reader.get_name() == "#text")
      reader.read();
    while (reader.get_name() != "taxon" && reader.get_name() != "taxa") { //while we are in a taxon block
      string rank = reader.get_name();
      reader.read();
      if (reader.get_name() == rank) { //if tag contents are empty
	reader.read();
	if (reader.get_name() == "#text")
	  reader.read();
	continue;
      }
      string info = reader.read_string();
      if (rank == "d") ::taxa_info[ind]->d = info;
      else if (rank == "k") ::taxa_info[ind]->k = info;
      else if (rank == "p") ::taxa_info[ind]->p = info;
      else if (rank == "c") ::taxa_info[ind]->c = info;
      else if (rank == "o") ::taxa_info[ind]->o = info;
      else if (rank == "f") ::taxa_info[ind]->f = info;
      else if (rank == "g") ::taxa_info[ind]->g = info;
      else if (rank == "tr") {
	for (unsigned int i=0; i<info.size(); i++) {
	  string tempstr = "";
	  tempstr.push_back(info[i]);
	  ::taxa_info[ind]->traits.push_back(stoi(tempstr));
	}
      }
      else {
	cerr << "Error in XML parse of '" << filename << "'. '" << rank << "' does not denote a valid taxonomic rank or trait node. Valid taxonomic ranks are d, k, p, c, o, f, and g." << endl;
	exit(3);
      }
      reader.read();
      reader.read();
      if (reader.get_name() == "#text")
	reader.read();
    }
  }
  cout << "Classification data loaded successfully!" << endl;
}	


int ANTLR3_CDECL parse_input_file(string inputstring, bool flag ){ 

    // Now we declare the ANTLR related local variables we need.
    // Note that unless you are convinced you will never need thread safe
    // versions for your project, then you should always create such things
    // as instance variables for each invocation.
    // -------------------

    // Name of the input file. Note that we always use the abstract type pANTLR3_UINT8
    // for ASCII/8 bit strings - the runtime library garauntees that this will be
    // good on all platforms. This is a general rule - always use the ANTLR3 supplied
    // typedefs for pointers/types/etc.
    //
    pANTLR3_UINT8	    fName;

    // The ANTLR3 character input stream, which abstracts the input source such that
    // it is easy to privide inpput from different sources such as files, or 
    // memory strings.
    //
    // For an ASCII/latin-1 memory string use:
    //	    input = antlr3NewAsciiStringInPlaceStream (stringtouse, (ANTLR3_UINT64) length, NULL);
    //
    // For a UCS2 (16 bit) memory string use:
    //	    input = antlr3NewUCS2StringInPlaceStream (stringtouse, (ANTLR3_UINT64) length, NULL);
    //
    // For input from a file, see code below
    //
    // Note that this is essentially a pointer to a structure containing pointers to functions.
    // You can create your own input stream type (copy one of the existing ones) and override any
    // individual function by installing your own pointer after you have created the standard 
    // version.
    //
    pANTLR3_INPUT_STREAM    input;

    // The lexer is of course generated by ANTLR, and so the lexer type is not upper case.
    // The lexer is supplied with a pANTLR3_INPUT_STREAM from whence it consumes its
    // input and generates a token stream as output.
    //
    ppqlLexer		    lxr;

    // The token stream is produced by the ANTLR3 generated lexer. Again it is a structure based
    // API/Object, which you can customise and override methods of as you wish. a Token stream is
    // supplied to the generated parser, and you can write your own token stream and pass this in
    // if you wish.
    //
    pANTLR3_COMMON_TOKEN_STREAM	    tstream;

    // The C parser is also generated by ANTLR and accepts a token stream as explained
    // above. The token stream can be any source in fact, so long as it implements the 
    // ANTLR3_TOKEN_SOURCE interface. In this case the parser does not return anything
    // but it can of cousre speficy any kind of return type from the rule you invoke
    // when calling it.
    //
    ppqlParser				psr;

    // Create the input stream based upon the argument supplied to us on the command line
    // for this example, the input will always default to ./input if there is no explicit
    // argument.
  
	char *a = new char[inputstring.size()+1];
	a[inputstring.size()] = 0;
	memcpy(a,inputstring.c_str(),inputstring.size());

	fName	= (pANTLR3_UINT8)a;
  

    // Create the input stream using the supplied file name
    // (Use antlr3AsciiFileStreamNew for UCS2/16bit input).
    //
    if (!flag) {
      input	= antlr3AsciiFileStreamNew(fName);
    }
    else {
      input = antlr3NewAsciiStringCopyStream (fName, (ANTLR3_UINT64) inputstring.size(), NULL);
    }

    // The input will be created successfully, providing that there is enough
    // memory and the file exists etc
    //
    if ( input == NULL )
    {
	    fprintf(stderr, "Unable to open file %s\n", (char *)fName);
	    exit(ANTLR3_ERR_NOMEM);
    }

    // Our input stream is now open and all set to go, so we can create a new instance of our
    // lexer and set the lexer input to our input stream:
    //  (file | memory | ?) --> inputstream -> lexer --> tokenstream --> parser ( --> treeparser )?
    //
    lxr	    = pqlLexerNew(input);	    // CLexerNew is generated by ANTLR

    // Need to check for errors
    //
    if (lxr == NULL )
    {
	    fprintf(stderr, "Unable to cuniquereate the lexer due to malloc() failure1\n");
	    exit(ANTLR3_ERR_NOMEM);
	}

    // Our lexer is in place, so we can create the token stream from it
    // NB: Nothing happens yet other than the file has been read. We are just 
    // connecting all these things together and they will be invoked when we
    // call the parser rule. ANTLR3_SIZE_HINT can be left at the default usually
    // unless you have a very large token stream/input. Each generated lexer
    // provides a token source interface, which is the second argument to the
    // token stream creator.
    // Note tha even if you implement your own token structure, it will always
    // contain a standard common token within it and this is the pointer that
    // you pass around to everything else. A common token as a pointer within
    // it that should point to your own outer token structure.
    //
    tstream = antlr3CommonTokenStreamSourceNew(ANTLR3_SIZE_HINT, TOKENSOURCE(lxr));

    if (tstream == NULL)
    {
		fprintf(stderr, "Out of memory trying to allocate token stream\n");
		exit(ANTLR3_ERR_NOMEM);
    }

    // Finally, now that we have our lexer constructed, we can create the parser
    //
    psr	    = pqlParserNew(tstream);  // CParserNew is generated by ANTLR3

    if (tstream == NULL)
    {
		fprintf(stderr, "Out of memory trying to allocate parser\n");
		exit(ANTLR3_ERR_NOMEM);
    }

    // We are all ready to go. Though that looked complicated at first glance,
    // I am sure, you will see that in fact most of the code above is dealing
    // with errors and there isn;t really that much to do (isn;t this always the
    // case in C? ;-).
    //
    // So, we now invoke the parser. All elements of ANTLR3 generated C components
    // as well as the ANTLR C runtime library itself are pseudo objects. This means
    // that they are represented as pointers to structures, which contain any
    // instance data they need, and a set of pointers to other interfaces or
    // 'methods'. Note that in general, these few pointers we have created here are
    // the only things you will ever explicitly free() as everythign else is created
    // via factories, that alloacte memory efficiently and free() everything they use
    // automatically when you close the parser/lexer/etc.
    //
    // Note that this means only that the methods are always called via the object
    // pointer and the first argument to any method, is a pointer to the structure itself.
    // It also has the side advantage, if you are using an IDE such as VS2005 that can do it,
    // that when you type ->, you will see a list of all the methods the object supports.
    //
	debugstatement("Entering the parser");
	psr->prog(psr);	
	debugstatement("Exited parser");

    //free it all. 
    psr	    ->free  (psr);	     	psr = NULL;
    tstream ->free  (tstream);	    tstream = NULL;
    lxr	    ->free  (lxr);	    	lxr = NULL;
    input   ->close (input);	    input = NULL;
    delete[] a;
    return 0;  
}



int ANTLR3_CDECL call_parser(string inputstring, bool flag ){ 

    // Now we declare the ANTLR related local variables we need.
    // Note that unless you are convinced you will never need thread safe
    // versions for your project, then you should always create such things
    // as instance variables for each invocation.
    // -------------------

    // Name of the input file. Note that we always use the abstract type pANTLR3_UINT8
    // for ASCII/8 bit strings - the runtime library garauntees that this will be
    // good on all platforms. This is a general rule - always use the ANTLR3 supplied
    // typedefs for pointers/types/etc.
    //
    pANTLR3_UINT8	    fName;

    // The ANTLR3 character input stream, which abstracts the input source such that
    // it is easy to privide inpput from different sources such as files, or 
    // memory strings.
    //
    // For an ASCII/latin-1 memory string use:
    //	    input = antlr3NewAsciiStringInPlaceStream (stringtouse, (ANTLR3_UINT64) length, NULL);
    //
    // For a UCS2 (16 bit) memory string use:
    //	    input = antlr3NewUCS2StringInPlaceStream (stringtouse, (ANTLR3_UINT64) length, NULL);
    //
    // For input from a file, see code below
    //
    // Note that this is essentially a pointer to a structure containing pointers to functions.
    // You can create your own input stream type (copy one of the existing ones) and override any
    // individual function by installing your own pointer after you have created the standard 
    // version.
    //
    pANTLR3_INPUT_STREAM    input;

    // The lexer is of course generated by ANTLR, and so the lexer type is not upper case.
    // The lexer is supplied with a pANTLR3_INPUT_STREAM from whence it consumes its
    // input and generates a token stream as output.
    //
    ppqlLexer		    lxr;

    // The token stream is produced by the ANTLR3 generated lexer. Again it is a structure based
    // API/Object, which you can customise and override methods of as you wish. a Token stream is
    // supplied to the generated parser, and you can write your own token stream and pass this in
    // if you wish.
    //
    pANTLR3_COMMON_TOKEN_STREAM	    tstream;

    // The C parser is also generated by ANTLR and accepts a token stream as explained
    // above. The token stream can be any source in fact, so long as it implements the 
    // ANTLR3_TOKEN_SOURCE interface. In this case the parser does not return anything
    // but it can of cousre speficy any kind of return type from the rule you invoke
    // when calling it.
    //
    ppqlParser				psr;

    // Create the input stream based upon the argument supplied to us on the command line
    // for this example, the input will always default to ./input if there is no explicit
    // argument.
  
  
	char *a = new char[inputstring.size()+1];
	a[inputstring.size()] = 0;
	memcpy(a,inputstring.c_str(),inputstring.size());

	fName	= (pANTLR3_UINT8)a;
  

    // Create the input stream using the supplied file name
    // (Use antlr3AsciiFileStreamNew for UCS2/16bit input).
    //
    if (!flag) {
      input	= antlr3AsciiFileStreamNew(fName);
    }
    else {
      input = antlr3NewAsciiStringCopyStream (fName, (ANTLR3_UINT64) inputstring.size(), NULL);
    }

    // The input will be created successfully, providing that there is enough
    // memory and the file exists etc
    //
    if ( input == NULL )
    {
	    fprintf(stderr, "Unable to open file %s\n", (char *)fName);
	    exit(ANTLR3_ERR_NOMEM);
    }

    // Our input stream is now open and all set to go, so we can create a new instance of our
    // lexer and set the lexer input to our input stream:
    //  (file | memory | ?) --> inputstream -> lexer --> tokenstream --> parser ( --> treeparser )?
    //
    lxr	    = pqlLexerNew(input);	    // CLexerNew is generated by ANTLR

    // Need to check for errors
    //
    if (lxr == NULL )
    {
	    fprintf(stderr, "Unable to cuniquereate the lexer due to malloc() failure1\n");
	    exit(ANTLR3_ERR_NOMEM);
	}

    // Our lexer is in place, so we can create the token stream from it
    // NB: Nothing happens yet other than the file has been read. We are just 
    // connecting all these things together and they will be invoked when we
    // call the parser rule. ANTLR3_SIZE_HINT can be left at the default usually
    // unless you have a very large token stream/input. Each generated lexer
    // provides a token source interface, which is the second argument to the
    // token stream creator.
    // Note tha even if you implement your own token structure, it will always
    // contain a standard common token within it and this is the pointer that
    // you pass around to everything else. A common token as a pointer within
    // it that should point to your own outer token structure.
    //
    tstream = antlr3CommonTokenStreamSourceNew(ANTLR3_SIZE_HINT, TOKENSOURCE(lxr));

    if (tstream == NULL)
    {
		fprintf(stderr, "Out of memory trying to allocate token stream\n");
		exit(ANTLR3_ERR_NOMEM);
    }

    // Finally, now that we have our lexer constructed, we can create the parser
    //
    psr	    = pqlParserNew(tstream);  // CParserNew is generated by ANTLR3

    if (tstream == NULL)
    {
		fprintf(stderr, "Out of memory trying to allocate parser\n");
		exit(ANTLR3_ERR_NOMEM);
    }

    // We are all ready to go. Though that looked complicated at first glance,
    // I am sure, you will see that in fact most of the code above is dealing
    // with errors and there isn;t really that much to do (isn;t this always the
    // case in C? ;-).
    //
    // So, we now invoke the parser. All elements of ANTLR3 generated C components
    // as well as the ANTLR C runtime library itself are pseudo objects. This means
    // that they are represented as pointers to structures, which contain any
    // instance data they need, and a set of pointers to other interfaces or
    // 'methods'. Note that in general, these few pointers we have created here are
    // the only things you will ever explicitly free() as everythign else is created
    // via factories, that alloacte memory efficiently and free() everything they use
    // automatically when you close the parser/lexer/etc.
    //
    // Note that this means only that the methods are always called via the object
    // pointer and the first argument to any method, is a pointer to the structure itself.
    // It also has the side advantage, if you are using an IDE such as VS2005 that can do it,
    // that when you type ->, you will see a list of all the methods the object supports.
    //
	debugstatement("Entering the parser");
	psr->prog(psr);	
	debugstatement("Exited parser");

    //free it all. 
    psr	    ->free  (psr);	     	psr = NULL;
    tstream ->free  (tstream);	    tstream = NULL;
    lxr	    ->free  (lxr);	    	lxr = NULL;
    input   ->close (input);	    input = NULL;
    delete[] a;
    return 0;  
}




double ANTLR3_CDECL call_parser(string inputstring ){ 
	
    pANTLR3_INPUT_STREAM input = antlr3NewAsciiStringCopyStream((pANTLR3_UINT8)inputstring.c_str(), (ANTLR3_UINT64) inputstring.size(), NULL);

    ppqlLexer lxr;

    pANTLR3_COMMON_TOKEN_STREAM tstream;

    ppqlParser psr;

    if ( input == NULL ){
	    cout<<"Input string is null" <<endl; 
    }

    lxr = pqlLexerNew(input);

    if (lxr == NULL ){
	    fprintf(stderr, "Unable to cuniquereate the lexer due to malloc() failure1\n");
	    exit(ANTLR3_ERR_NOMEM);
	}

    tstream = antlr3CommonTokenStreamSourceNew(ANTLR3_SIZE_HINT, TOKENSOURCE(lxr));

    if (tstream == NULL)
    {
		fprintf(stderr, "Out of memory trying to allocate token stream\n");
		exit(ANTLR3_ERR_NOMEM);
    }

    // Finally, now that we have our lexer constructed, we can create the parser
    //
    psr	    = pqlParserNew(tstream);  // CParserNew is generated by ANTLR3

    if (tstream == NULL)
    {
		fprintf(stderr, "Out of memory trying to allocate parser\n");
		exit(ANTLR3_ERR_NOMEM);
    }

  	debugstatement("Entering the parser");
  	start_clock();
	psr->prog(psr);	
	double TimeOfRun = stop_clockbp();
	debugstatement("Exited parser");

    //free it all. 
    psr	    ->free  (psr);	     	psr = NULL;
    tstream ->free  (tstream);	    tstream = NULL;
    lxr	    ->free  (lxr);	    	lxr = NULL;
    input   ->close (input);	    input = NULL;
    return TimeOfRun;  
}

int free_things(){
	/*
	while(!::biparttable.taxa_in_trees.empty()){
		delete[] ::biparttable.taxa_in_trees.back();
		::biparttable.taxa_in_trees.pop_back();
	} 
	*/
	
	//~ while(!::list_bs.empty()){
		//~ delete[] ::list_bs.back();
		//~ ::list_bs.pop_back();
	//~ }
	
	while(!::query_results.empty()){
		delete ::query_results.back();
		::query_results.pop_back();
	}
	
	for ( std::map<std::string, pqlsymbol * >::const_iterator i = symbol_table.begin(); i != symbol_table.end(); ++i ){
		delete i->second;
		//i->second = 0; // I don't think this is strictly necessary...
	}
	
	return 1;
	
}

void init_the_constants(){
  ::NUM_TREES_INIT = ::NUM_TREES;
	for(unsigned int i = 0; i < NUM_TREES; i++){
		::all_trees.insert(i);
	        ::original_trees.insert(i);
	}

	//initialize a Taxon object in taxa_info for each taxon
	for (unsigned int i = 0; i < ::NUM_TAXA; i++) {
	  ::taxa_info.push_back(new Taxon);
	  ::taxa_info[i]->label = ::biparttable.lm.name(i);
	}

	//Make taxa names constants.  
	vector<string> taxaset;
	for (unsigned int i = 0; i < ::biparttable.lm.size(); i ++ ){
		if(isalpha(::biparttable.lm.name(i)[0])){
			symbol_table[::biparttable.lm.name(i)] = new pqlsymbol(::biparttable.lm.name(i));
			constant_table[::biparttable.lm.name(i)] = true;
			taxaset.push_back(::biparttable.lm.name(i));
		}
		else{
			string str = "t";
			str.append(::biparttable.lm.name(i));
			symbol_table[str] = new pqlsymbol(::biparttable.lm.name(i));
			constant_table[str] = true;
			taxaset.push_back(str);
		}
		
		//cout << lm.name(i) << endl;
	}
	symbol_table["true"] = new pqlsymbol(true);
	constant_table["true"] = true;
	symbol_table["false"] = new pqlsymbol(false);
	constant_table["false"] = true;
	symbol_table["taxa"] = new pqlsymbol(taxaset);
	constant_table["taxa"] = true;
	symbol_table["trees"] = new pqlsymbol(all_trees);
	constant_table["trees"] = true;
	symbol_table["original_trees"] = new pqlsymbol(original_trees);
	constant_table["original_trees"] = true;
	
	//for (unsigned int i=0; i < ::NUM_TAXA; ++i) 
	//	::shuffledTaxa.push_back(i); // 1 2 3 4 5 6 7 8 9	
}

void stripwhite (char *string){
  register int i = 0;

  while (whitespace (string[i]))
    i++;

  if (i)
    strcpy (string, string + i);

  i = strlen (string) - 1;

  while (i > 0 && whitespace (string[i]))
    i--;

  string[++i] = '\0';
}


int main(int argc, char **argv){

  srand(time(NULL));

	
	init_output();
	rl_attempted_completion_function = my_completion;
	//DEBUGMODE = true;
	DEBUGMODE = false;
		
	unsigned int QueryResultsSize = query_results.size();

	if(::DEBUGMODE){
		write_to_output("Debug Mode is on\n");
		write_to_output("Session started at ");
		write_to_output(get_time_stamp());
		write_to_output("\n");		
		
		if (::HETERO){
			write_to_output( "HETERO = true\n");
		}
		else{
			write_to_output( "HETERO = false\n");
		}
		cout << "HETERO = " << HETERO << endl;
	} 
   
	if (argc != 3 && argc != 4){
	  cerr << "usage: ./sample trzfile commandfile" << endl;
	  cerr << "   or: ./sample trzfile classificationfile commandfile" << endl;
	  cerr << "This program searches the trzfile for the trees matching the search commands in the command file" << endl;
	  return 2;
	}

	start_clock();

	// Get the name of the database
	string trzfilename = argv[1];
	//std::clock_t start1 = std::clock();
	// Unpack the trees into the hash table
	load_data_from_trz_file(trzfilename);
  
	cout << "ParsingTime = " << stop_clockbp() << endl;

	start_clock();
	init_the_constants(); // This loads the constant varibles
	init_the_functs(); // This loads the language functions into usable data structures.
	
	
	//The main input file for the pql commands 

	string inputfilename; //, interactive_input;
 
	int cmdind = 2; // where commandfile is in argv

	if (argc == 4) { // if there is a classification file...
	  load_classification(argv[2]); // load it
	  cmdind = 3; // commandfile must be argv[3]
	}
	cout << "SetupTime = " << stop_clockbp() << endl;

	// Get path of TreeHouse directory, change to it
	//string thpath = argv[0];
	//thpath.resize(thpath.size() - 9);
	//chdir(thpath.c_str());
	//if (chdir(thpath.c_str()) == -1){
	//	cout << "setting the path to the current TreeHouse directory failed." << strerror (errno) << endl;
	//}

	//start_clock();
	if (strcmp(argv[cmdind],"-i") == 0) {
		testDistance();
		//TESTSTUFF();
//		TestClust();
		//TestDist();
		interactive = true;
		interactive_log.open("logs/interactive_log.txt");
		if(interactive_log){
			debugstatement("interactive log opened");
		}

		// Get path of TreeHouse directory, change to i
		string thpath = argv[0];
		thpath.resize(thpath.size() - 19);
		if (chdir(thpath.c_str()) == -1){
			cout << "setting the path to the current TreeHouse directory failed." << endl;
		}
		
		interactive_log  << "Session opened at " << get_time_stamp() << endl;
		bool done = 0;
		while (!done){
			char *line;
			line = readline ("TH>: ");
			stripwhite (line);
			string interactive_input = line;
			interactive_input.append("\n");

			if (!line || interactive_input == "quit()\n"){
				done = 1;             /* Encountered EOF at top level. */
				interactive_log  << "Session closed at " << get_time_stamp() << endl;
				interactive_log.close();
			}
			else{
				if (*line){	
					double TimeOfRun = call_parser(interactive_input);
					if( QueryResultsSize < query_results.size() ){
						QueryResultsSize++;
						interactive_log << "Command: "<< interactive_input;
						interactive_log << "Command issued at: " << get_time_stamp();
						interactive_log << "Result:" << query_results.back()->value_to_string() << endl;
						interactive_log << "Runtime:" << TimeOfRun << endl;
						interactive_log << endl;
					}
					else{
						interactive_log << "Input: "<< interactive_input;
					}			
					add_history (line);
				}
			}
			if (line)
				free (line);
		}
	}
  
	else {		
		interactive_log.open("logs/batch_log.txt");
		interactive_log  << "Session opened at " << get_time_stamp() << endl;
		interactive_log  << "Batch file : " << argv[cmdind] << endl;
		
		ifstream ifs( argv[cmdind] );
		string batch_input;

		// Get path of TreeHouse directory, change to it
		string thpath = argv[0];
		thpath.resize(thpath.size() - 19);
		//chdir(thpath.c_str());
		if (chdir(thpath.c_str()) == -1){
			cout << "setting the path to the current TreeHouse directory failed." << endl;
		}
		
		while(getline( ifs, batch_input ) ) {
			batch_input.append("\n");
			double TimeOfRun = call_parser(batch_input);
			if( QueryResultsSize < query_results.size() ){
				QueryResultsSize++;
				interactive_log << "Command: "<< batch_input;
				interactive_log << "Command issued at: " << get_time_stamp();
				interactive_log << "Result:" << query_results.back()->value_to_string() << endl;
				interactive_log << "Runtime:" << TimeOfRun << endl;
				interactive_log << endl;
			}
			else{
				interactive_log << "Input: "<< batch_input;
			}
		}
		interactive_log  << "Session closed at " << get_time_stamp() << endl;
		interactive_log.close();
	}

	
	//cout << "ExeTime = " << stop_clockbp() << endl;

	//cout << "SetTime = " << ::SetTime << endl;


	if(DEBUGMODE){
		::biparttable.lm.printMap();
		//print_hashtable();
		
		std::cout<< "taxa_in_trees = " << ::biparttable.taxa_in_trees.size() << endl;
		//print_taxa_in_trees();
		//std::cout<< "all_branches = " << ::all_branches.size() << endl;
		//std::cout<< "bs_sizes = " << ::bs_sizes.size() << endl;
		std::cout<< "all_trees = " << ::all_trees.size() << endl;

		//std::cout<< "length_of_bitstrings = " << ::length_of_bitstrings.size() << endl;
		//std::cout<< "list_bs = " << ::list_bs.size() << endl;
		//std::cout<< "list_branches = " << ::list_branches.size() << endl;
		//std::cout<< "all_bs = " << ::all_bs.size() << endl;
		std::cout<< "inverted_index = " << ::inverted_index.size() << endl;

		//::biparttable.print_searchtable();
		::biparttable.print_biparttable();
		//for (unsigned int x = 0; x < NUM_TREES; x++){
		//cout << compute_tree(::lm, ::all_bs[x], ::all_branches[x], x, 0, ::bs_sizes[x]) << endl;
		//}
	}
	
  free_things();

  return 0;

}
