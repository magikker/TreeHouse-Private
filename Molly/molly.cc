//*****************************************************************/
/*This is HashSim, a fast topological similarity matrix algorithm. It is based on 
Fast HashRF, which is in turn based on HashRF, by SeungJin Sul.

(c) 2010 HashSim, Fast HashRF: Suzanne Matthews
(c) 2009 HashRF: SeungJin Sul
*/
/*****************************************************/

#include <sys/time.h>
#include <sys/resource.h>
#include <cassert>
#include <valarray>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <algorithm>

#include "./tclap/CmdLine.h"
#include "global.h"

// For newick parser
extern "C" {
#include <newick.h>
}

using namespace std;
unsigned int * n_biparts;
bool VERBOSE = false;
ofstream fout;

bool intcomp (unsigned int i,unsigned int j) { return (i<j); }

int max(int x, int y){ 
  if (x>y)
    return x;
  else
    return y;
} 

unsigned int * maxval_a; //for proper 'a' calculation for approximate method

void process_DB(string infilename, 
		unsigned int * row, bool * isCompact){
   tree_counter = 0; //this determines if the trees are multifurcating or binary
  struct timeval ragged_start;
  struct timeval ragged_end;
  string fileLine, line_type;
  int pos;
  ifstream fin;
  unsigned int item;
  fin.open(infilename.c_str());
  if (!fin){
    cerr << "Cannot open file for reading!" << endl;
    exit(1);
  }

  //read in and generate compact table
  bool * check = (bool*)malloc(NUM_TREES*sizeof(bool));
  if ((check == NULL) ){
    cerr << "Out of memory!" << endl;
    exit(3);
  }
  for (unsigned int i =0; i < NUM_TREES; i++)
    n_biparts[i] = 0;
  for (unsigned int i =0; i < NUM_TREES; i++)
    check[i] = 0;

  //read in ALL_COUNT
  for (unsigned int i = 0; i < 3; i++)
    getline(fin, fileLine);
  pos = fileLine.find_first_of(":");
  line_type = fileLine.substr(0, pos);
  if (line_type != "ALL_COUNT"){
    cerr << "ALL_COUNT specification not found!";
    exit(2);
  }
  fileLine = fileLine.substr(pos+1);
  pos = fileLine.find_first_of("\n");
  fileLine = fileLine.substr(0,pos);
  ALL_COUNT = atoi(fileLine.c_str());

 //read in row
  getline(fin, fileLine);
  pos = fileLine.find_first_of(":");
  line_type = fileLine.substr(0, pos);
  if (line_type != "ROW"){
    cerr << "objects specification not found!";
    exit(2);
  }
  fileLine = fileLine.substr(pos+1);
  pos = fileLine.find_first_of("\n");
  fileLine = fileLine.substr(0,pos);

  cout  <<"NUM_TREES is: " << NUM_TREES << endl;
  fprintf(stderr, "Loading in correction array..\n");
  istringstream splitarray;
  splitarray.str(fileLine);
  for (unsigned int i = 0; i < NUM_TREES; i++){
    splitarray >> item;
    row[i]  = item;
  }
  //allocate compact table, compact table lengths, inverted index, and inverted index lengths
  hashtable = (unsigned int **)malloc(NUMBIPART*sizeof(unsigned int*));
  hash_lengths = (unsigned int *)malloc(NUMBIPART*sizeof(unsigned int));
  helper = (unsigned int **)malloc(NUM_TREES*sizeof(unsigned int*));
  helper_sizes = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));

  for (unsigned int i =0; i < NUM_TREES; i++)
    helper_sizes[i] = 0;

  //now load in compact table
  getline(fin, fileLine);
  pos = fileLine.find_first_of("\n");
  fileLine = fileLine.substr(0, pos);
  if (fileLine != "BEGIN COMPACT TABLE"){
    cerr << "cannot find compact table declaration!\n";
    exit(2);
  }
  int temp;
  string num;
  gettimeofday(&ragged_start, NULL);
  fprintf(stderr, "Loading in compact table..\n");
  for (unsigned int i = 0; i < NUMBIPART; i++){
    getline(fin, fileLine);
    pos = fileLine.find_first_of(":");
    line_type = fileLine.substr(0, pos); 
    fileLine = fileLine.substr(pos+1);
    pos = fileLine.find_first_of(":");
    num = fileLine.substr(0,pos);
    temp = atoi(num.c_str());
    hash_lengths[i] = temp;
    fileLine = fileLine.substr(pos+1);
    pos = fileLine.find_first_of("\n");
    fileLine = fileLine.substr(0,pos);
    hashtable[i] = (unsigned int*)malloc(temp*sizeof(unsigned int));
    if (hashtable[i] == NULL){
      cerr << "Cannot allocate memory!" << endl;
      exit(2);
    }
    splitarray.str(fileLine);
    if (line_type == "+"){
      isCompact[i] = 0;
      for (int j = 0; j < temp; j++){
	splitarray >> item;
	n_biparts[item]++;
	helper_sizes[item]++;
	hashtable[i][j] = item;
      }
    }
    else{
      isCompact[i] = 1;
      //ALL_COUNT++;
      for (int j = 0; j < temp; j++){
	splitarray >> item;
	check[item] = 1;
	helper_sizes[item]++;
	hashtable[i][j] = item;
      }
      for (unsigned int j = 0; j < NUM_TREES; j++){
	if (check[j]== 0){
	  n_biparts[j]++;
	}
	else
	  check[j] = 0;
      }
    }
  } 

  //allocate inverted index
  unsigned int isize = 0;
  for (unsigned int i = 0; i < NUM_TREES; i++){
    isize = helper_sizes[i];
    helper[i] = (unsigned int *)malloc(isize*sizeof(unsigned int));
    if (helper[i] == NULL){
      cerr << "out of memory!" << endl;
      exit(0);
    }
    helper_sizes[i] = 0;
  }

  cout << "populating index:" << endl;
  //populate inverted index
  unsigned int hsize, ipos, hitem;
  for (unsigned int i = 0; i < NUMBIPART; i++){
    hsize = hash_lengths[i];
    for (unsigned int j =0; j < hsize; j++){
      hitem = hashtable[i][j];
      ipos = helper_sizes[hitem];
      helper[hitem][ipos] = i;
      helper_sizes[hitem]++;
    }
  }
  gettimeofday(&ragged_end, NULL);
  double ragged_time = ragged_end.tv_sec - ragged_start.tv_sec + (ragged_end.tv_usec - ragged_start.tv_usec) / 1.e6;

  fprintf(stderr, "\nCompact Table Loading Time: %g s\n", ragged_time);
  fprintf(stderr, "Features in DB: %d\n", NUMBIPART);
  fprintf(stderr, "Objects in DB: %d\n", NUM_TREES);
  //fprintf(stderr, "Compacted Lines: %d\n", ALL_COUNT);
  
//now print out numbiparts.txt if necessary (needed for multifurcation)
  free(check);
  fin.close();
}

void process_DB2(string infilename, bool * isCompact){
   tree_counter = 0; //this determines if the trees are multifurcating or binary
  struct timeval ragged_start;
  struct timeval ragged_end;
  struct timeval inverted_start;
  struct timeval inverted_end;

  string fileLine, line_type;
  int pos;
  ifstream fin;
  unsigned int item;
  istringstream splitarray;
  fin.open(infilename.c_str());
  if (!fin){
    cerr << "Cannot open file for reading!" << endl;
    exit(1);
  }

  //read in and generate compact table
  bool * check = (bool*)malloc(NUM_TREES*sizeof(bool));
  if ((check == NULL) ){
    cerr << "Out of memory!" << endl;
    exit(3);
  }
  for (unsigned int i =0; i < NUM_TREES; i++)
    n_biparts[i] = 0;
  for (unsigned int i =0; i < NUM_TREES; i++)
    check[i] = 0;

  //read in ALL_COUNT
  for (unsigned int i = 0; i < 4; i++)
    getline(fin, fileLine);

  //allocate compact table, compact table lengths, inverted index, and inverted index lengths
  hashtable = (unsigned int **)malloc(NUMBIPART*sizeof(unsigned int*));
  hash_lengths = (unsigned int *)malloc(NUMBIPART*sizeof(unsigned int));
  helper = (unsigned int **)malloc(NUM_TREES*sizeof(unsigned int*));
  helper_sizes = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));

  for (unsigned int i =0; i < NUM_TREES; i++)
    helper_sizes[i] = 0;

  //now load in compact table
  getline(fin, fileLine);
  pos = fileLine.find_first_of("\n");
  fileLine = fileLine.substr(0, pos);
  if (fileLine != "BEGIN COMPACT TABLE"){
    cerr << "cannot find compact table declaration!\n";
    exit(2);
  }

  int temp;
  string num;
  gettimeofday(&ragged_start, NULL);
  fprintf(stderr, "Loading in compact table..\n");
  for (unsigned int i = 0; i < NUMBIPART; i++){
    getline(fin, fileLine);
    pos = fileLine.find_first_of(":");
    line_type = fileLine.substr(0, pos); 
    fileLine = fileLine.substr(pos+1);
    pos = fileLine.find_first_of(":");
    num = fileLine.substr(0,pos);
    temp = atoi(num.c_str());
    hash_lengths[i] = temp;
    fileLine = fileLine.substr(pos+1);
    pos = fileLine.find_first_of("\n");
    fileLine = fileLine.substr(0,pos);
    hashtable[i] = (unsigned int*)malloc(temp*sizeof(unsigned int));
    if (hashtable[i] == NULL){
      cerr << "Cannot allocate memory!" << endl;
      exit(2);
    }

    splitarray.str(fileLine);
    if (line_type == "+"){
      isCompact[i] = 0;
      for (int j = 0; j < temp; j++){
	splitarray >> item;
	n_biparts[item]++;
	helper_sizes[item]++;
	hashtable[i][j] = item;
      }
    }
    else{
      isCompact[i] = 1;
      for (int j = 0; j < temp; j++){
	splitarray >> item;
	check[item] = 1;
	hashtable[i][j] = item;
      }
      for (unsigned int j = 0; j < NUM_TREES; j++){
	if (check[j]== 0){
	  n_biparts[j]++;
	  helper_sizes[j]++;
	}
	else
	  check[j] = 0;
      }
    }
  } 

  gettimeofday(&ragged_end, NULL);
  double ragged_time = ragged_end.tv_sec - ragged_start.tv_sec + (ragged_end.tv_usec - ragged_start.tv_usec) / 1.e6;
  fprintf(stderr, "\nCompact Table Loading Time: %g s\n", ragged_time);

  //gettimeofday(&inverted_start, NULL);
  //allocate inverted index
  unsigned int isize = 0;
  for (unsigned int i = 0; i < NUM_TREES; i++){
    isize = helper_sizes[i];
    helper[i] = (unsigned int *)malloc(isize*sizeof(unsigned int));
    helper_sizes[i] = 0;
  }
  //populate inverted index
  unsigned int hsize, ipos, hitem;
  for (unsigned int i = 0; i < NUMBIPART; i++){
    hsize = hash_lengths[i];
    if (isCompact[i]){
      ALL_COUNT++;
      for (unsigned int j =0; j < hsize; j++){
	hitem = hashtable[i][j];
	check[hitem] = 1;
      }
      for (unsigned int j = 0; j < NUM_TREES; j++){ 
	if (check[j] == 0){
	  ipos = helper_sizes[j];
	  helper[j][ipos] = i;
	  helper_sizes[j]++;
	}
	else
	  check[j] = 0;
      }
    }
    else{ //not compact, procede as usual
      for (unsigned int j =0; j < hsize; j++){
	hitem = hashtable[i][j];
	ipos = helper_sizes[hitem];
	helper[hitem][ipos] = i;
	helper_sizes[hitem]++;
      }
    }
    
    free(hashtable[i]);
  }
  free(hashtable);
  free(hash_lengths);  
  gettimeofday(&inverted_end, NULL);
  double inverted_time = inverted_end.tv_sec - inverted_start.tv_sec + (inverted_end.tv_usec - inverted_start.tv_usec) / 1.e6;
  fprintf(stderr, "\nInverted Index Loading Time: %g s\n", ragged_time);
  fprintf(stderr, "Features in DB: %d\n", NUMBIPART);
  fprintf(stderr, "Objects in DB: %d\n", NUM_TREES);
  fprintf(stderr, "Compacted Lines: %d\n", ALL_COUNT);
  
//now print out numbiparts.txt if necessary (needed for multifurcation)
  free(check);
  fin.close();
}

void process_DB3(string infilename, bool * isCompact, vector<string> &molecule_ids){
   tree_counter = 0; //this determines if the trees are multifurcating or binary
  struct timeval inverted_start;
  struct timeval inverted_end;

  string fileLine, line_type;
  int pos;
  ifstream fin;
  unsigned int item;
  istringstream splitarray;
  fin.open(infilename.c_str());
  if (!fin){
    cerr << "Cannot open file for reading!" << endl;
    exit(1);
  }
 for (unsigned int i =0; i < NUM_TREES; i++)
    n_biparts[i] = 0;


 //allocate inverted index 
  helper = (unsigned int **)malloc(NUM_TREES*sizeof(unsigned int*));
  helper_sizes = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));

  for (unsigned int i =0; i < NUM_TREES; i++)
    helper_sizes[i] = 0;


  for (unsigned int i = 0; i < NUM_TREES; i++)
    isCompact[i] = 0;

  unsigned int isize, loc;
  string molecule;
  //now load in index
  gettimeofday(&inverted_start, NULL);
  for (unsigned int i = 0; i < NUM_TREES; i++){
    getline(fin, fileLine);
    pos = fileLine.find_first_of("\t");
    molecule = fileLine.substr(0, pos);
    molecule_ids.push_back(molecule); //add id to vector
    fileLine = fileLine.substr(pos+1);
    pos = fileLine.find_first_of("\n");
    fileLine = fileLine.substr(0, pos);
    //cout << "fileLine is : " << fileLine << endl;
    for(unsigned int j = 0; j < NUMBIPART; j++){
      if(fileLine[j] == '1'){
	n_biparts[i]++;
      }
    }
    isize = n_biparts[i];
    if (isize > (float)NUMBIPART/2){
      isize = NUMBIPART - isize;
      isCompact[i] = 1;
    }
    helper[i] = (unsigned int *)malloc(isize*sizeof(unsigned int));
    if (helper[i] == NULL){
      cerr << "Cannot allocate!" << endl;
      exit(2);
    }
    if (isCompact[i]){
      for (unsigned int j = 0; j < NUMBIPART; j++){
	if (fileLine[j] == '0'){
	  loc = helper_sizes[i];
	  helper[i][loc] = j;
	  helper_sizes[i]++;
	}
      }
    }
    else{
      for(unsigned int j = 0; j < NUMBIPART; j++){
	if(fileLine[j] == '1'){
	  loc = helper_sizes[i];
	  helper[i][loc] = j;
	  helper_sizes[i]++;
	}
      }
    }

  } 

  gettimeofday(&inverted_end, NULL);
  double inverted_time = inverted_end.tv_sec - inverted_start.tv_sec + (inverted_end.tv_usec - inverted_start.tv_usec) / 1.e6;
  fprintf(stderr, "\nIndex Loading Time: %g s\n", inverted_time);
  fprintf(stderr, "Features in DB: %d\n", NUMBIPART);
  fprintf(stderr, "Objects in DB: %d\n", NUM_TREES);
  fprintf(stderr, "Compacted Lines: %d\n", ALL_COUNT);
  
//now print out numbiparts.txt if necessary (needed for multifurcation)
//free(check);
  fin.close();
}

//performs a SINGLE query (aka, processes ONE bitstring!)
double processQuery(string mybitstring, unsigned int * row, string outfilename, bool * isCompact, double threshold) {
  unsigned int * avals = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  bool * bitstring = (bool *)malloc(NUMBIPART*sizeof(bool));
  if ((bitstring == NULL) || (avals == NULL)){
    cerr << "Out of memory!" << endl;
    exit(2);
  }

  struct timeval query_start;
  struct timeval query_end;

  string id;
  int pos = mybitstring.find_first_of("\t");
  id = mybitstring.substr(0,pos);
  mybitstring = mybitstring.substr(pos+1);
  pos = mybitstring.find_first_of("\n");
  mybitstring = mybitstring.substr(0,pos);
  if (mybitstring.length() != NUMBIPART){
    cerr << "Query String of wrong size!" << endl;
    exit(4);
  }
  for (unsigned int i = 0; i < NUM_TREES; i++)
    avals[i] = 0;

 if (PRINT){
    if (!fout) {
      cerr << "cannot open file for writing!" << endl;
      exit(2);
    }
    fout << "Molecule " <<mybitstring << " is most similar to:" << endl; 
  }



  //get relevant portion of bitstring and build localX
  unsigned int localX = 0;
  unsigned int  totalBiparts = 0;
 gettimeofday(&query_start, NULL);
  for (unsigned int i =0; i < NUMBIPART; i++){
    if (mybitstring[i] == '1'){
      bitstring[i] = 1;
      totalBiparts++;
    }
    else{
      bitstring[i] = 0;
      if (isCompact[i])
	localX++;
    }
  }

  unsigned int isize, sim, loc,a,b,c;
  sim = 0;
  float val, value;
  unsigned int mismatches, maxmismatch;
  double q;

  for (unsigned int i = 0; i < NUM_TREES; i++){
    val = (float)min(n_biparts[i], totalBiparts)/max(n_biparts[i], totalBiparts);
    if (val >= threshold){
      //cout << "Fingerprint " << i << endl;
      isize = helper_sizes[i]; //number of hashlines fingerprint appears on
      sim = 0;
      //mismatches = 0;
      //cout << "totalBiparts: " << totalBiparts << endl;
      //cout << "n_biparts: " << n_biparts[i] << endl;
      //q = (double)((threshold*(n_biparts[i] + totalBiparts)) + ((threshold-1)*(row[i] + localX - ALL_COUNT)))/(threshold+1);
      //maxmismatch = isize - q;
      //cout << "isize is: " << isize << endl;
      for (unsigned int j = 0; j < isize; j++){
	loc = helper[i][j]; //
	if (bitstring[loc] && !isCompact[loc])
	  sim++;
	else if (!bitstring[loc] && isCompact[loc])
	  sim++;
	//else 
	// mismatches++;
	//if (mismatches > maxmismatch)
	//  break;
      }

      //if (mismatches > maxmismatch)
      //continue;
      a = sim + ALL_COUNT - row[i] - localX; 
      b = totalBiparts - a;
      c = n_biparts[i] -a;
      //cout << i << ": [a,b,c] = [" << a <<","<< b <<"," << c <<"]" << endl; 
      value = (float)a/(a+b+c); //tanimoto
      if (PRINT){
	if (value >= threshold)
	  fout << "Molecule " << i << ": s=" << value << endl;
      }
    }
  }

  //now we can begin actual querying!
  /*unsigned int size, temp;
  for (unsigned int i = 0; i < NUMBIPART; i++){
    if (bitstring[i]){
     if (!isCompact[i]){
	size = hash_lengths[i];
	for (unsigned int j = 0; j < size; j++){
	  temp = hashtable[i][j];
	  avals[temp]++;
	}
     }
    }
    else{
      if (isCompact[i]){
	size = hash_lengths[i];
	for (unsigned int j = 0; j < size; j++){
	  temp = hashtable[i][j];
	  avals[temp]++;
	}
      }
    }
  }

  //ofstream fout; 
 
  //calculate Tanimoto and output 

  float value;
    if (PRINT)
    fout << "Molecule " << id << " is most similar to:" << endl; 
  //cout << "totalBiparts is: " << totalBiparts << endl;
  for (unsigned int i =0;  i < NUM_TREES; i++){
    a = avals[i] + ALL_COUNT - row[i] - localX; 
    b = totalBiparts - a;
    c = n_biparts[i] -a;

    //cout << i << ": [a,b,c] = [" << a <<","<< b <<"," << c <<"]" << endl; 
    value = (float)a/(a+b+c); //tanimoto
    if (PRINT){
      if (value >= threshold)
	fout << "Molecule " << i << ": s=" << value << endl;
    }
    }*/
  gettimeofday(&query_end, NULL);
  double query_time = query_end.tv_sec - query_start.tv_sec + (query_end.tv_usec - query_start.tv_usec) / 1.e6;
  return query_time;
  //fout.close();
}

double processQuery2(string mybitstring, string outfilename, double threshold, double & ave_features, bool * isCompact, vector<string> molecule_ids) {


  bool * bitstring = (bool *)malloc(NUMBIPART*sizeof(bool));
  if (bitstring == NULL){
    cerr << "Out of memory!" << endl;
    exit(2);
  }

  struct timeval query_start;
  struct timeval query_end;
  //struct timeval t1_start, t1_end, t2_start, t2_end, t3_start, t3_end;

  ave_features = 0.0;
  string id;
  int pos = mybitstring.find_first_of("\t");
  id = mybitstring.substr(0,pos);
  mybitstring = mybitstring.substr(pos+1);
  pos = mybitstring.find_first_of("\n");
  mybitstring = mybitstring.substr(0,pos);
  if (mybitstring.length() != NUMBIPART){
    cerr << "Query String of wrong size!" << endl;
    exit(4);
  }

 if (PRINT){


    if (!fout) {
      cerr << "cannot open file for writing!" << endl;
      exit(2);
    }
  }


 /*unsigned int * reg = (unsigned int *)malloc(NUMBIPART * sizeof(unsigned int));
 unsigned int * comp = (unsigned int *)malloc(NUMBIPART * sizeof(unsigned int));
 if ((reg == NULL) || (comp == NULL)){
   cerr << "cannot allocate!" << endl;
   exit(2);
 }
 unsigned int sizereg = 0;
 unsigned int sizecomp = 0;*/

 gettimeofday(&query_start, NULL);


  //get relevant portion of bitstring and build localX
  unsigned int  totalBiparts = 0;
  //gettimeofday(&t1_start, NULL);
  for (unsigned int i =0; i < NUMBIPART; i++){
    if (mybitstring[i] == '1'){
      bitstring[i] = 1;
      totalBiparts++;
    }
    else{
      bitstring[i] = 0;
    }
  }
  //gettimeofday(&t1_end, NULL);

  unsigned int isize, loc,a,b,c, count, ave_count;
  float val, value, tmp;
  a = 0;
  ave_count  = 0;
  unsigned int mismatches, maxmismatch;
  //bool compact= false;
  //int mismatches;
  //double t2_time, temp_time;
  //t2_time = 0;
  //if (PRINT)
    //fout << "Molecule " << mybitstring << " is most similar to:" << endl; 
  // cout << "mybitstring is: " << mybitstring << endl;

  //gettimeofday(&t3_start, NULL);
  for (unsigned int i = 0; i < NUM_TREES; i++){

    val = (float)min(n_biparts[i], totalBiparts)/max(n_biparts[i], totalBiparts);
    if (val >= threshold){//if this is a fingerprint we should care about
      //gettimeofday(&t2_start, NULL);
      //cout << "We will check fingerprint " << i;
      isize = helper_sizes[i]; //number of hashlines fingerprint appears on
      a = 0;
      //compact = isCompact[i];
      if (totalBiparts > n_biparts[i]){
	tmp = (float)threshold*totalBiparts;
	maxmismatch = totalBiparts - tmp;
      }
      else{
	tmp = (float)threshold*n_biparts[i];
	maxmismatch =n_biparts[i] - tmp;	
      }
      //tmp = (float)threshold*n_biparts[i];
      //maxmismatch = n_biparts[i] - tmp;

      //cout << "tmp is: " << tmp << endl;  
      mismatches = 0;
      count = 0;
      if (!isCompact[i]){
	//maxmismatch = totalBiparts - tmp;
	//cout << "( which is not compact )" << endl;
	//cout << "maxmismatch is: " << maxmismatch << endl;
	//cout << "threshold is: " << threshold << endl;
	//cout << "n_biparts is: " << n_biparts[i] << endl;
	//cout << "totalBiparts is: " << totalBiparts << endl;
	
	for (unsigned int j = 0; j < isize; j++){
	  loc = helper[i][j];
	  //cout << "bitstring[" << loc << "] is: " << bitstring[loc] << endl;
	  if (bitstring[loc]){
	    //cout << "match!" << endl;
	    count++;
	    a++;
	  }
	  else {
	    mismatches++;
	    //cout << "mismatch! Now we have " << mismatches << " mismatches." << endl;
	  }
	  //if ((isize - j) < (maxmismatch - mismatches))
	  //  break;	 
	  if (mismatches > maxmismatch)
	    break;
	}
      }
      else{
	//cout << "( which is compact )" << endl;
	//cout << "maxmismatch is: " << maxmismatch << endl;
	//cout << "threshold is: " << threshold << endl;
	//cout << "n_biparts is: " << n_biparts[i] << endl;
	//cout << "totalBiparts is: " << totalBiparts << endl;
	for (unsigned int j = 0; j < isize; j++){
	  loc = helper[i][j];
	  count++; 
	  //cout << "bitstring[" << loc << "] is: " << bitstring[loc] << endl;
	  if (!bitstring[loc]){
	    //cout << "match!" << endl;
	    a++; //this is actually d
	  }
	  else {
	    mismatches++;
	    //cout << "mismatch! Now we have " << mismatches << " mismatches." << endl;
	  }
	  //if ((isize - j) < (maxmismatch - mismatches))
	  //  break;	 
	  if (mismatches > maxmismatch)
	    break;
	}
      }
      //cout << "The number of total mismatches was: " << mismatches << " (max is: " << maxmismatch << ")" << endl;
      if (mismatches > maxmismatch)
	continue;
      //cout << "we do get here!!!" << endl;
      if (isCompact[i]){
	//b = (NUMBIPART-totalBiparts) - a;	
	//c = (NUMBIPART-n_biparts[i]) -a;
	a = totalBiparts-NUMBIPART + n_biparts[i] + a;
	b = totalBiparts-a;
	c = n_biparts[i]-a;
      }
      else{
	b = totalBiparts-a;
	c = n_biparts[i]-a;
      }
      //cout << i << ": [a,b,c] = [" << a <<","<< b <<"," << c <<"]" << endl; 
      value = (float)a/(a+b+c); //tanimoto
      //cout << "count is: " << count << endl;
      ave_features += count;
      ave_count++;    
      //cout << "value is: " << value << endl;
      //cout << "threshold is: " << threshold << endl;
      if (PRINT){
	if (value >= threshold){
	  //cout << "do we get here?" << endl;
	  fout << molecule_ids[i] << " " << value << endl;
	}
      }
      //exit(0);
      //gettimeofday(&t2_end, NULL);   
      //temp_time = t2_end.tv_sec - t2_start.tv_sec + (t2_end.tv_usec - t2_start.tv_usec) / 1.e6;
      //t2_time += temp_time;
    }
  }
  //gettimeofday(&t3_end, NULL);
  //cout << "ave_features is: " << ave_features << endl;
  //cout << "ave_count is: " << ave_count << endl;
  ave_features = (double)ave_features/ave_count;
  //cout << "ave_features is: " << ave_features << endl;
  gettimeofday(&query_end, NULL);
  double query_time = query_end.tv_sec - query_start.tv_sec + (query_end.tv_usec - query_start.tv_usec) / 1.e6;
  /*double t1_time = t1_end.tv_sec - t1_start.tv_sec + (t1_end.tv_usec - t1_start.tv_usec) / 1.e6;
  //double t2_time = t2_end.tv_sec - t2_start.tv_sec + (t2_end.tv_usec - t2_start.tv_usec) / 1.e6;
  double t3_time = t3_end.tv_sec - t3_start.tv_sec + (t3_end.tv_usec - t3_start.tv_usec) / 1.e6;
  cout << "t1 time: " << t1_time << endl;
  cout << "t2 time: " << t2_time << endl;
  cout << "t3 time: " << t3_time << endl;
  cout << "query time: " << query_time << endl;
  exit(0);*/
  return query_time;
}


/**********************
 *
 * main function
 *
 *********************/


int main(int argc, char** argv){
  string outfilename = "query.out";
  string bfilename;
  string infilename;
  double query = 0.9;

  struct timeval prog_start;
  struct timeval prog_end;
  struct timeval process_start;
  struct timeval process_end;
  gettimeofday(&prog_start, NULL);
  // TCLAP
  //since we are creating the hashtable serially, the assumption in this
  //version is that the seed file is the same as the input file and that 
  //offset is set to 0
  //this can readily be changed in future version if necessary
  try {

    // Define the command line object.
    string helpMsg  = "molly <DB file> <bitstring file> -o <output> \n";
    helpMsg += "Input file: \n";
    helpMsg += "   The DB file is prepared separately.\n";
    helpMsg += "The bitstring file is a series of bitstrings. Each bitstring is a query object made up of a set of features, represented by 0s and 1s.\n"; 
    helpMsg += "It is necessary that the length of this bitstring is the same as the number of features in the table.\n";
    helpMsg += "You can specify output using the -o option. Molecular ids are outputted if they have a Tanimoto similarity greater than 0.9.\n";
    //helpMsg += " The seed file is necessary for parallel processing. A common seed file is necessary to ensure labels are properly populated\n"; 
    //helpMsg += " -d <distance measure>: 0..6. 0 is Robinson-Foulds distance, and is used as default.\n"; 
    helpMsg += "Examples: \n";
    helpMsg += " molly  foo.db bitstring.file //input is a database called foo.db, the bitstrings for querying is in bitstring.file. Default output is used, which is query.out\n";
    helpMsg += " molly  foo.db bitstring.file -o myout.txt //Same as above, but now output is myout.txt\n";

    TCLAP::CmdLine cmd(helpMsg, ' ', "1.0.0");

    TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "database name", true, "database", "Input database name"  );
    cmd.add( fnameArg );

    TCLAP::UnlabeledValueArg<string>  bnameArg( "bname", "bitstring name", true, "bitstrings", "bitstring query name"  );
    cmd.add( bnameArg );

    TCLAP::ValueArg<double>  qArg( "q", "querythresh", "query threshold", false, 0.9, "query threshold to use"  );
    cmd.add( qArg );

    TCLAP::SwitchArg pArg("p", "print", "print to file?", false);
    cmd.add( pArg );

    TCLAP::SwitchArg vArg("v", "verbose", "verbose output?", false);
    cmd.add( vArg );

    TCLAP::ValueArg<string> oArg("o", "output", "output file", false, "query.out", "output file name");
    cmd.add( oArg );


    cmd.parse( argc, argv );

    if (strcmp(oArg.getValue().c_str(), "query.out")!= 0)
      outfilename = oArg.getValue();

    infilename = fnameArg.getValue();
    bfilename = bnameArg.getValue();
    //if (qArg.getValue())
    query = qArg.getValue();
    if (pArg.getValue())
      PRINT = true;
    if (vArg.getValue())
      VERBOSE = true;
  } catch (TCLAP::ArgException &e) { // catch any exceptions
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
  cout << "PRINT is: " << PRINT << endl;
  unsigned int * row; //stores the amount to be subtracted from rows at the end
  bool * isCompact; //keeps track of which rows are compact and which are not.
  ifstream fin;
  fin.open(infilename.c_str());
  string fileLine, line_type, item;
  //read in number of features (NUMBIPART)
  getline(fin, fileLine);
  int pos = fileLine.find_first_of("\t");
  fileLine = fileLine.substr(pos+1);
  pos = fileLine.find_first_of("\n");
  fileLine = fileLine.substr(0, pos);
  NUMBIPART = fileLine.size();
  //read in number of objects (NUM_TREES)
  NUM_TREES = 1;
  while (getline(fin, fileLine)) 
    NUM_TREES++;
  fin.close();
  cout << "NUMBER OF OBJECTS: " << NUM_TREES << endl;
  cout<< "NUMBER OF FEATURES:" << NUMBIPART << endl;

  isCompact = (bool *)malloc(NUM_TREES * sizeof(bool));
  if (isCompact == NULL){
    cerr << "Out of memory!" << endl;
    exit(2);
  }
 
  //row = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  n_biparts = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  cout << "Query Threshold: " << query << endl; 
  gettimeofday(&process_start, NULL);
  //process_DB(infilename, row, isCompact);
  vector<string> molecule_ids;
  process_DB3(infilename, isCompact, molecule_ids);

  gettimeofday(&process_end, NULL);
  double process_time = process_end.tv_sec - process_start.tv_sec + (process_end.tv_usec - process_start.tv_usec) / 1.e6;
   fprintf(stderr,"\nInput Processing Time: %g s\n", process_time);
   

   //debugging
  
   /* cout << "hashtable" << endl;
   cout << "hash_lengths: " << endl;
  for (unsigned int i = 0; i < NUMBIPART; i++) { 
    cout << hash_lengths[i] << " ";
  }
  cout << endl;

  for (unsigned int i = 0; i < NUMBIPART; i++) { 
    cout << "H[" << i << "]: ";
    for (unsigned int j = 0; j < hash_lengths[i]; j++) { 
      cout << hashtable[i][j] << " ";
    }
    cout << endl;
    } 
   */
   /*
 cout << "inverted index" << endl;
 cout << "index_lengths: " << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++) { 
    cout << helper_sizes[i] << " ";
  }
  cout << endl;
 
  for (unsigned int i = 0; i < NUM_TREES; i++) { 
    cout << "I[" << i << ":";
    if (isCompact[i])
      cout << "-: ";
    else 
      cout << "+: ";
    for (unsigned int j = 0; j < helper_sizes[i]; j++) { 
      cout << helper[i][j] << " ";
    }
    cout << endl;
  } 

  //cout << "row array: " << endl;
  //for (unsigned int i = 0; i < NUM_TREES; i++) { 
  //  cout << row[i] << " ";
  //}
  //cout << endl;
  //cout << "ALL_COUNT: " << ALL_COUNT << endl;
  //cout << "NUMBER OF UNIQUE BIPARTITIONS: " << NUMBIPART << endl;
  
//end debugging
//exit(0);
*/

   ifstream fin2;
   fin2.open(bfilename.c_str());
   if (!fin2){
     cerr << "cannot open file for reading!" << endl;
     exit(2);
   }

   fout.open(outfilename.c_str());
   if (!fout){
     cerr << "cannot open file for writing!" << endl;
     exit(2);
   }
   string myBitstring;
   unsigned int count = 1;
   double qtime, average, min, max, qcount, qtotal;
   average = 0.0;
   qtotal = 0.0;
   fprintf(stderr, "Beginning Querying...\n");
   
   if (query == 0){
     struct timeval myquery_start;
     struct timeval myquery_end;
     average = 0;
     while (!fin2.eof()){
       getline(fin2, myBitstring);
       if (!fin2.eof()){
	 if (PRINT){
	   gettimeofday(&myquery_start, NULL);
	   //fout << "Molecule " << myBitstring << " is most similar to:" << endl; 
	   for (unsigned int i = 0; i < NUM_TREES; i++)
	     fout << molecule_ids[i] << endl;
	   gettimeofday(&myquery_end, NULL);
	   qtime = myquery_end.tv_sec - myquery_start.tv_sec + (myquery_end.tv_usec - myquery_start.tv_usec) / 1.e6;
	   average += qtime;
	   if (count == 1){ 
	     min = qtime;
	     max = qtime;
	   }
	   if (qtime > max)
	     max = qtime;
	   if (qtime < min)
	     min = qtime;
	   count++;
	 }
       }
     }
     average = average/count;
     if (!PRINT){
       count = 0;
       min = 0;
       max = 0;
       qtime = 0;
     }
   }
   else{
     cout << "We get here!" << endl;
     while (!fin2.eof()){
       getline(fin2, myBitstring);
       if (!fin2.eof()){
	 if (VERBOSE)
	   cout << "Query " << count << endl; 
	 //qtime = processQuery(myBitstring, row, outfilename, isCompact, query);
	 
	 qtime = processQuery2(myBitstring, outfilename, query, qcount, isCompact, molecule_ids);
	 qtotal += qcount;
	 //cout << "Average number of features checked: " << qcount << endl;
	 if (VERBOSE)
	   cout << "Time: " << qtime << "s" << endl;
	 if (count == 1){
	   min = qtime;
	   max = qtime;
	 }
	 else{
	   if (qtime < min)
	     min = qtime;
	   if (qtime > max)
	     max = qtime;
	 }
	 average += qtime;
	 count++;
       }
       //exit(0);
     }
   }
   if (PRINT)
     fout.close();
   count--;
   average = (double)average/count;
   qtotal = (double)qtotal/count;
   fprintf(stderr, "***Querying Statistics:***\n");
   fprintf(stderr, "Number of Queries: %u\n", count);
   fprintf(stderr, "Average Query Time %g s\n", average); 
   fprintf(stderr, "Min Query Time %g s\n", min);
   fprintf(stderr, "Max Query Time %g s\n\n", max);  
   fprintf(stderr, "Average number of features checked: %g\n", qtotal);

   gettimeofday(&prog_end, NULL);
   double prog_time = prog_end.tv_sec - prog_start.tv_sec + (prog_end.tv_usec - prog_start.tv_usec) / 1.e6;
   fprintf(stderr,"\nProgram Running Time: %g s\n", prog_time);
   //fprintf(stderr, "do we get here?\n");
  return 0;

}

