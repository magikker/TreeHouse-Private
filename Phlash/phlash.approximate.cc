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
#include <cmath>

#include "label-map.hh"
#include "hashfunc.hh"
#include "hash.hh"
#include "./tclap/CmdLine.h"
#include "global.h"
#include "parsing.h"

// For newick parser
extern "C" {
#include <newick.h>
}

using namespace std;
bool APPROXIMATE = false;
bool MISMATCH = false;
bool ROOTED = false;
bool FULLMATRIX = false;
bool QOPT = false;
bool TRIVIAL = false;
unsigned int NUM_TREES2;
unsigned int TOTAL_TREES;
//vector<unsigned int> n_biparts;
unsigned int * n_biparts;

typedef struct {
  bool *my_bs;
  unsigned int h1; //hashtable location
  unsigned int pos; //chain position --> usually 0
} bipartition; //for approximate method --> definition of a bipartition

typedef struct { 
  unsigned int val;
  bool add;
} checks; //for approximation method --> decides whether certain elements should be added
bool intcomp (unsigned int i,unsigned int j) { return (i<j); }

int max(int x, int y){ 
  if (x>y)
    return x;
  else
    return y;
} 

unsigned int * maxval_a; //for proper 'a' calculation for approximate method
 
void generate_approximate_matrix(unsigned int approximate_threshold){
  //allocate maxval vector
  vector< vector<unsigned int> > full_index;

  struct timeval allocate_start;
  struct timeval allocate_end;

  struct timeval first_start;
  struct timeval first_end;

  struct timeval second_start;
  struct timeval second_end;

  struct timeval third_start;
  struct timeval third_end;

  
  gettimeofday(&allocate_start, NULL);
  maxval_a = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  cout << "approximate threshold is: " << approximate_threshold << endl;

  //create full index
  full_index.resize(NUM_TREES);
  unsigned int index_pos = 0;
  for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); ++hti) {
    unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size(); 
    if (sizeVec) {
      for (unsigned int i=0; i<sizeVec; ++i) {
	unsigned int sizeTreeIdx = vvec_hashrf._hashtab2[hti][i]._vec_treeidx.size();
	for (unsigned int j =0; j < sizeTreeIdx; ++j) { //first determine the elements that exist
	  int checker = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
	  full_index[checker].push_back(index_pos);
	} 
	index_pos++;
      } //for everything in sizevec
    }
  } //end full index creation
  
  NUMBIPART = index_pos; //index_pos is now the number of bipartitions.

  fprintf(stderr, "generating approximate matching bipartition matrix...\n");
  //allocate and initialize matrix
  bool ** bipart_matrix = (bool **)malloc(NUMBIPART * sizeof(bool*)); 
  for (unsigned int i = 0; i < NUMBIPART; i++)
    bipart_matrix[i] = (bool *)malloc(NUMBIPART*sizeof(bool));
  for (unsigned int i = 0; i < NUMBIPART; i++){
    for (unsigned int j = 0; j < NUMBIPART; j++){
      bipart_matrix[i][j] = 0;
    }
  }

  gettimeofday(&allocate_end, NULL);
  double allocate_time = allocate_end.tv_sec - allocate_start.tv_sec + (allocate_end.tv_usec - allocate_start.tv_usec) / 1.e6;
  fprintf(stderr, "Allocation Bipartition Matrix Time: %g s\n", allocate_time);

  //collect bipartitions
  gettimeofday(&allocate_start, NULL);
  vector< bipartition > new_bs;
  for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); ++hti) {
    unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size(); 
    if (sizeVec) {
      for (unsigned int i=0; i<sizeVec; ++i) {
	bipartition mybipart;
	mybipart.my_bs = vvec_hashrf._hashtab2[hti][i]._bs;
	mybipart.h1 = hti;
	mybipart.pos = i;
	new_bs.push_back(mybipart);	 
      }
    }
  }
  gettimeofday(&allocate_end, NULL);
  allocate_time = allocate_end.tv_sec - allocate_start.tv_sec + (allocate_end.tv_usec - allocate_start.tv_usec) / 1.e6;
  fprintf(stderr, "Collect Bipartitions Time: %g s\n", allocate_time);
  assert(new_bs.size() == NUMBIPART);

  //print list of bipartitions:
  /*for (unsigned int i = 0; i < NUMBIPART; i++){
    cout << i << ":  ";
    for (unsigned int x = BITSETSZ-1; x > BITSETSZ-NUM_TAXA-1; x--){
      cout << (*new_bs[i].my_bs)[x];
    }
    cout << endl;
    }*/
  //populate matrix:

  gettimeofday(&allocate_start, NULL);

  unsigned int sim = 0;
  for (unsigned int i = 0; i < NUMBIPART; i++){
    bool *first = new_bs[i].my_bs; //set first bitstring
    for (unsigned int j = i+1; j < NUMBIPART; j++){
      sim = 0;
      bool *second = new_bs[j].my_bs; //set second bitstring
      for (unsigned int x = 0; x < NUM_TAXA; x++){
	if (first[x] == second[x]) 
	  sim++;
      }
      if (sim < (NUM_TAXA-sim))
	sim = NUM_TAXA - sim;
      if (sim >= approximate_threshold){
	bipart_matrix[i][j] = sim;
	bipart_matrix[j][i] = sim;
      }
    }
  }


  //use inverted index to update hash table.
  /*cout << "full index data structure: " << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++) {
    cout << "Tree " << i << ": ";
    for (unsigned int j = 0; j < full_index[i].size(); j++) {
      cout << full_index[i][j] << " ";
    }
    cout << endl;
    }*/

  gettimeofday(&allocate_end, NULL);
  allocate_time = allocate_end.tv_sec - allocate_start.tv_sec + (allocate_end.tv_usec - allocate_start.tv_usec) / 1.e6;
  fprintf(stderr, "Update Bipartition Matrix Time: %g s\n", allocate_time);
  
  //for each tree
  gettimeofday(&allocate_start, NULL);
  unsigned int mysize = 0;
  unsigned int mypos = 0;
  unsigned int myval = 0;
  unsigned int compos  = 0;
  vector< unsigned int > to_compare;
  vector< int > checks;
  unsigned int a, b;
  double first_time, second_time, third_time;
  unsigned int h1, chainpos;
  first_time = 0;
  second_time = 0;
  third_time = 0;

  for (unsigned int i = 0; i < NUM_TREES; i++){
    mysize = full_index[i].size();
    to_compare.resize(NUMBIPART-mysize);
    mypos = 0;
    myval = 0;
    myval = full_index[i][mypos];
    compos = 0;
    //populate to_compare
    gettimeofday(&first_start, NULL);
    for (unsigned int j = 0; j < NUMBIPART; j++){
      if ( j  == myval){ //if j is equal to myval
	mypos++; //increment mypos
	myval = full_index[i][mypos]; //update myval
      }
      else { //if j is not equal to myval
	to_compare[compos] = j;
	//to_compare[compos].add = 0;
	compos++;
      }
    }
    gettimeofday(&first_end, NULL);    
    double temp_time = first_end.tv_sec - first_start.tv_sec + (first_end.tv_usec - first_start.tv_usec) / 1.e6;
    first_time+= temp_time;
    /* cout << "to compare is: " << endl;
    for (unsigned int x = 0; x < to_compare.size(); x++){
      cout << to_compare[x] << " ";
    }
    cout << endl;*/
    checks.resize(to_compare.size()+1);

    gettimeofday(&second_start, NULL);
    //cout << "checks is: " <<  endl;
    for (unsigned int x =0; x < to_compare.size(); x++){
      checks[x] = x;
      //cout << checks[x] << " ";
    }
    checks[to_compare.size()] = -1;
    //cout << checks[to_compare.size()] << endl;

    int current, next;
    unsigned int x, y; 
    for (unsigned int j = 0; j < mysize; j++){
      a = full_index[i][j];
      x = 0;
      y = 1;
      while ( (x != to_compare.size()) ){
	current = checks[x];
	next = checks[y];
	if (current < 0 )
	  break;
	b = to_compare[current];
	if (bipart_matrix[a][b]){
	  gettimeofday(&third_start, NULL);
	  checks[x] = next;
	  //cout << "adding tree " << i << " to bipartition " << b << endl;
	  h1 = new_bs[b].h1;  // get b's h1 value
	  chainpos = new_bs[b].pos; //get b's chain position
	  vvec_hashrf._hashtab2[h1][chainpos]._vec_treeidx.push_back(i); //add tree i to bipartition b's tree list 
	  //full_index[i].push_back(b); //note that this tree will contain this bipartition
	  gettimeofday(&third_end, NULL);
	  temp_time = third_end.tv_sec - third_start.tv_sec + (third_end.tv_usec - third_start.tv_usec) / 1.e6;
	  third_time += temp_time;
	}
	else{
	  x++;
	  checks[x] = next;
	}
	if (next > 0)
	  y++;
      }
    }
	
    full_index[i].clear(); 
    gettimeofday(&second_end, NULL);
    temp_time = second_end.tv_sec - second_start.tv_sec + (second_end.tv_usec - second_start.tv_usec) / 1.e6;
    second_time+= temp_time;
  } //end for all trees
  full_index.clear();
  gettimeofday(&allocate_end, NULL);
  allocate_time = allocate_end.tv_sec - allocate_start.tv_sec + (allocate_end.tv_usec - allocate_start.tv_usec) / 1.e6;
  fprintf(stderr, "Create compute_a data structure Time: %g s\n", first_time);
  fprintf(stderr, "Use compute_a data structure to update data structures: %g s\n", second_time);
  fprintf(stderr, "Amount of time to update data structures: %g s\n", third_time);
  fprintf(stderr, "Traverse and Update Hashtable etc. Time: %g s\n", allocate_time);

  //deallocate bipartition matrix
  for (unsigned int i = 0; i < NUMBIPART; i++)
    free(bipart_matrix[i]);
  free(bipart_matrix);

  //deallocate integer array
  checks.clear();
  to_compare.clear();

  //update hashtable with information in full index
  /*unsigned int arrsize = 0; 
  unsigned int h1, chainpos, bp;
  for (unsigned int i = 0; i < NUM_TREES; i++){
    arrsize = full_index[i].size();
    for (unsigned int j = 0; j < arrsize; j++){
      bp = full_index[i][j];
      h1 = new_bs[bp].h1;  // get b's h1 value
      chainpos = new_bs[bp].pos; //get b's chain position
      vvec_hashrf._hashtab2[h1][chainpos]._vec_treeidx.push_back(i); //add tree i to bipartition b's tree list 
    }
    }*/

  //deallocate new_bs
  new_bs.clear();
}

void process_trees(string infilename, 
   string qfilename,
   unsigned int & hashtable_length, 
   unsigned int * row,
   unsigned int approximate_threshold) {


  unsigned int hashtable_pos = 0;
   unsigned int hashtable_length_correction = 0;
  ofstream fout; //we will use this if we need to print out the hash table
  unsigned long long M1=0;
  unsigned long long M2=0;
  FILE *fp;
  NEWICKTREE *newickTree;
  int err;
  unsigned long uBID = 0;
  tree_counter = 0; //this determines if the trees are multifurcating or binary
  struct timeval index_start;
  struct timeval index_end;
  struct timeval hash_start;
  struct timeval hash_end;
  struct timeval ragged_start;
  struct timeval ragged_end;
  struct timeval approximate_start;
  struct timeval approximate_end;
  double approximate_time = 0;

  //counter information for amazing 
  gettimeofday(&index_start, NULL);
  //allocate helper data structures
  unsigned int hashtablesize = 0;
  if (QOPT)
    hashtablesize = TOTAL_TREES;
  else
    hashtablesize = NUM_TREES;

  helper = (unsigned int **)malloc(hashtablesize*sizeof(unsigned int *));
  helper_sizes = (unsigned int *)malloc(hashtablesize*sizeof(unsigned int));

  for (unsigned int i =0; i < hashtablesize; ++i) { 
    helper_sizes[i]= 0;
  }

  gettimeofday(&index_end, NULL);
  double temp_time = index_end.tv_sec - index_start.tv_sec + (index_end.tv_usec - index_start.tv_usec) / 1.e6;
  index_time+=temp_time;


  //generate Hashtable
  gettimeofday(&hash_start, NULL);
  LabelMap lm = collect_labels(infilename);
  if (QOPT){
    NUM_TREES = TOTAL_TREES;
    initialize_hashtable(M1, M2); //initialize contents of hashtable
    NUM_TREES = TOTAL_TREES - NUM_TREES2;
  }
  else
    initialize_hashtable(M1, M2); //initialize contents of hashtable

  fp = fopen(infilename.c_str(), "r");
  if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}
  if (MC){  //do MC procedure
    for (unsigned int treeIdx=0; treeIdx<NUM_TREES; ++treeIdx) {
      newickTree = loadnewicktree2(fp, &err);
      if (!newickTree) {
	switch (err) {
	case -1:
	  printf("Out of memory\n");
	  break;
	case -2:
	  printf("parse error\n");
	  break;
	case -3:
	  printf("Can't load file\n");
	  break;
	default:
	  printf("Error %d\n", err);
	}
	exit(0);
      }
      else {
	unsigned int numBitstr=0;
	bool processed = false;
	if (ROOTED)
	  dfs_compute_hash_mc(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, NULL);
	else
	  dfs_compute_hash_unrooted_mc(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, processed, NULL);
	if (MAXVAL == 0){
	  if (ROOTED)
	    MAXVAL = numBitstr -1;
	  else
	    MAXVAL = NUM_TAXA - 3;
	}
	killnewicktree(newickTree);
      }
    }
  }
  else {
    for (unsigned int treeIdx=0; treeIdx<NUM_TREES; ++treeIdx) {
      newickTree = loadnewicktree2(fp, &err);
      if (!newickTree) {
	switch (err) {
	case -1:
	  printf("Out of memory\n");
	  break;
	case -2:
	  printf("parse error\n");
	  break;
	case -3:
	  printf("Can't load file\n");
	  break;
	default:
	  printf("Error %d\n", err);
	}
	exit(0);
      }
      else {
	unsigned int numBitstr=0;
	bool processed = false;
	if (ROOTED)
	  dfs_compute_hash(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, NULL);
	else
	  dfs_compute_hash_unrooted2(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, processed, NULL);
	if (MAXVAL == 0){
	  if (ROOTED)
	    MAXVAL = numBitstr -1;
	  else
	    MAXVAL = NUM_TAXA - 3;
	}
	killnewicktree(newickTree);
      }
    }
  }

  fclose(fp);

  if (QOPT){
    fp = fopen(qfilename.c_str(), "r");

    if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}
    if (MC){ //do MC procedure
      for (unsigned int treeIdx=NUM_TREES; treeIdx< TOTAL_TREES; ++treeIdx) {
	newickTree = loadnewicktree2(fp, &err);
	if (!newickTree) {
	  switch (err) {
	  case -1:
	    printf("Out of memory\n");
	    break;
	  case -2:
	    printf("parse error\n");
	    break;
	  case -3:
	    printf("Can't load file\n");
	    break;
	  default:
	    printf("Error %d\n", err);
	  }
	  exit(0);
	}
	else {
	  unsigned int numBitstr=0;
	  bool processed = false;
	  if (ROOTED)
	    dfs_compute_hash_mc(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, NULL);
	  else
	    dfs_compute_hash_unrooted_mc(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, processed, NULL);
	  if (MAXVAL == 0){
	    if (ROOTED)
	      MAXVAL = numBitstr -1;
	    else
	      MAXVAL = NUM_TAXA - 3;
	  }
	  killnewicktree(newickTree);
	}
      }
    
      fclose(fp);
    }
    else{
      for (unsigned int treeIdx=NUM_TREES; treeIdx< TOTAL_TREES; ++treeIdx) {
	newickTree = loadnewicktree2(fp, &err);
	if (!newickTree) {
	  switch (err) {
	  case -1:
	    printf("Out of memory\n");
	    break;
	  case -2:
	    printf("parse error\n");
	    break;
	  case -3:
	    printf("Can't load file\n");
	    break;
	  default:
	    printf("Error %d\n", err);
	  }
	  exit(0);
	}
	else {
	  unsigned int numBitstr=0;
	  bool processed = false;
	  if (ROOTED)
	    dfs_compute_hash(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, NULL);
	  else
	    dfs_compute_hash_unrooted2(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, processed, NULL);
	  if (MAXVAL == 0){
	    if (ROOTED)
	      MAXVAL = numBitstr -1;
	    else
	      MAXVAL = NUM_TAXA - 3;
	  }
	  killnewicktree(newickTree);
	}
      } 
      fclose(fp);
    }
  }
  gettimeofday(&hash_end, NULL);
  double hash_time = hash_end.tv_sec - hash_start.tv_sec + (hash_end.tv_usec - hash_start.tv_usec) / 1.e6;
  fprintf(stderr, "Hash table Creation Time: %g s\n", hash_time);

  //if approximate is requested, call that function:
  if (APPROXIMATE){
    gettimeofday(&approximate_start, NULL);
    if (MISMATCH)
      approximate_threshold = NUM_TAXA - approximate_threshold;
    else
      approximate_threshold = (NUM_TAXA*approximate_threshold)/100;
    generate_approximate_matrix(approximate_threshold);
    gettimeofday(&approximate_end, NULL);
    double temp_time = approximate_end.tv_sec - approximate_start.tv_sec + (approximate_end.tv_usec - approximate_start.tv_usec) / 1.e6;
    approximate_time += temp_time;
   }  

  if (QOPT)
    NUM_TREES = TOTAL_TREES;

  vector<bool> check(NUM_TREES, false);
  unsigned int threshold = NUM_TREES - (NUM_TREES*P)/100;
  n_biparts = (unsigned int *)malloc(NUM_TREES * sizeof(unsigned int));
  for (unsigned int i = 0; i < NUM_TREES; i++)
    n_biparts[i] = 0;

  //count up the number of unique bipartitions
  unsigned int unique_bipartitions = 0; 
  for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); ++hti) {
    unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size(); 
    if (sizeVec) {
      unique_bipartitions += sizeVec;      
      gettimeofday(&index_start, NULL);
      for (unsigned int i=0; i<sizeVec; ++i) {
	unsigned int sizeTreeIdx = vvec_hashrf._hashtab2[hti][i]._vec_treeidx.size();	  
	if( sizeTreeIdx < NUM_TREES) { //if it's less than the total number of trees	     
	  if (sizeTreeIdx > threshold) {  //if it is greater than the threshold value
	    for (unsigned int j =0; j < sizeTreeIdx; ++j) { //first determine the elements that exist
	      int checker = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
	      check[checker] = true;
	    } 
	    for (unsigned int j = 0; j < NUM_TREES; ++j) { //then print out/process the ones that don't
	      if (check[j] == false) {
		helper_sizes[j]++;
	      }	      
	      else {
		check[j] = false;
	      }
	    }
	  } //end if compress and greater than threshold
	  else { //if it less than the threshold, print/process as usual
	    if (sizeTreeIdx > 1) { 
	      for (unsigned int j=0; j<sizeTreeIdx; ++j) {	
		unsigned int temp  = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
		helper_sizes[temp]++;
	      }	      
	    } //we don't care about single element lists, except for d calculation
	  }
	}
      } //end for everything in sizeVec
      gettimeofday(&index_end, NULL);
      temp_time = index_end.tv_sec - index_start.tv_sec + (index_end.tv_usec - index_start.tv_usec) / 1.e6;
      index_time+=temp_time;  
    } //end if sizevec
  } //end index creation

  NUMBIPART = unique_bipartitions;

  //allocate helper data structures
  gettimeofday(&index_start, NULL);
  //now, allocate helper
  unsigned int newval = 0;
  for (unsigned int i = 0; i < NUM_TREES; ++i) { 
    newval = helper_sizes[i];
    helper[i] = (unsigned int*)malloc(newval*sizeof(unsigned int));
    for (unsigned int j = 0; j < newval; ++j) { 
      helper[i][j] = 0;
    }
    helper_sizes[i]=0;
   } 
  gettimeofday(&index_end, NULL);
  temp_time = index_end.tv_sec - index_start.tv_sec + (index_end.tv_usec - index_start.tv_usec) / 1.e6;
  index_time+=temp_time;
  
  hashtable = (unsigned int **)malloc(unique_bipartitions*sizeof(unsigned int *));
  hash_lengths = (unsigned int*)malloc(unique_bipartitions*sizeof(unsigned int)); 
  for( unsigned int i = 0; i < unique_bipartitions; ++i) {
    hash_lengths[i] = 0;
  }
  //stores the length of elements in hash table 
  hashtable_length = unique_bipartitions;
  
  gettimeofday(&ragged_start, NULL);
  unsigned int amazing_loc = 0;	   
  temp_time = 0;   
  for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); ++hti) {
    unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size(); 
    if (sizeVec) {
      uBID += sizeVec;
      for (unsigned int i=0; i<sizeVec; ++i) {
	unsigned int sizeTreeIdx = vvec_hashrf._hashtab2[hti][i]._vec_treeidx.size();
	tree_counter += sizeTreeIdx;
	//cout << hti << "." << vvec_hashrf._hashtab2[hti][i]._hv2 << ":" << sizeTreeIdx << ":";
	bool increment_row_array = false;
	if( sizeTreeIdx < NUM_TREES) {	    
	  if ( (sizeTreeIdx > threshold) ) { //if it is greater than the threshold value
	    ALL_COUNT++; 
	    increment_row_array = true;
	    int allocate_size = NUM_TREES - sizeTreeIdx + 1;
	    hashtable[hashtable_pos] = (unsigned int*)malloc(allocate_size*sizeof(unsigned int)); //allocate an array of this size
	    hash_lengths[hashtable_pos] = allocate_size;
	    unsigned int pos = 1; //where we add the first element
	    unsigned int lowest = NUM_TREES;
	    for (unsigned int j =0; j < sizeTreeIdx; ++j) { //first determine the elements that exist
	      int checker = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
	      //cout << checker << " ";
	      n_biparts[checker]++;
	      check[checker] = true;
	    } 
	    //cout << endl;
	    for (unsigned int j = 0; j < NUM_TREES; ++j) { //then print out/process the ones that don't
	      if (check[j] == false) {
		hashtable[hashtable_pos][pos] = j;
		gettimeofday(&index_start, NULL);
		amazing_loc = helper_sizes[j];
		helper[j][amazing_loc] = hashtable_pos;
		helper_sizes[j]++;
		gettimeofday(&index_end, NULL);
		temp_time = index_end.tv_sec - index_start.tv_sec + (index_end.tv_usec - index_start.tv_usec) / 1.e6;
		index_time+=temp_time;
		pos++;
		if(j < lowest) 
		  lowest = j;
		if(increment_row_array)
		  row[j]++;
	      }
	      else {
		check[j] = false;
	      }
	    }
	    //update lowest element, and increment hashtable_pos
	    hashtable[hashtable_pos][0] = lowest;
	    hashtable_pos++;
	  } //end if greater than threshold
	  else { //if it less than the threshold, print/process as usual
	    if (sizeTreeIdx > 1) {
	      hashtable[hashtable_pos] = (unsigned int *)malloc((sizeTreeIdx+1)*sizeof(unsigned int));
	      hash_lengths[hashtable_pos] = sizeTreeIdx+1;
	      unsigned int pos = 1; //where we begin to increment
	      unsigned int lowest = NUM_TREES;
	      for (unsigned int j=0; j<sizeTreeIdx; ++j) {	
		unsigned int temp  = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
		//cout << temp <<  " ";
		n_biparts[temp]++;
		hashtable[hashtable_pos][pos] = temp;
		gettimeofday(&index_start, NULL);
		amazing_loc = helper_sizes[temp];
		helper[temp][amazing_loc] = hashtable_pos;
		helper_sizes[temp]++;
		gettimeofday(&index_end, NULL);
		temp_time = index_end.tv_sec - index_start.tv_sec + (index_end.tv_usec - index_start.tv_usec) / 1.e6;
		index_time+=temp_time;
		pos++;
		if (temp < lowest) 
		  lowest = temp;
	      }
	      //cout << endl;
	      hashtable[hashtable_pos][0] = lowest;
	      hashtable_pos++;
	    }//we silently ignore hash lines of length one.
	    else{
	      hashtable_length_correction++;
	      assert(sizeTreeIdx !=0);
	      //cout << vvec_hashrf._hashtab2[hti][i]._vec_treeidx[0] << endl;
	    }
	  } //end else if it's less than the threshold
	}
	else { //if sizeTreeIdx == NUM_TREES
	  for (unsigned int j=0; j<sizeTreeIdx; ++j) {	
	    int temp  = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
	    //cout << temp << " ";
	    n_biparts[temp]++;
	  }
	  //cout << endl;
	  ALL_COUNT++;
	  hashtable_length_correction++;
	}
      } //for everything in sizevec
      //vvec_hashrf._hashtab2[hti].clear();  //clear the contents of the array
    }
  } //end hashtable creation
  //vvec_hashrf._hashtab2.clear(); //clear out this vector as well. 
  //do by-tree traversal to generate index

  assert(uBID == unique_bipartitions);
  hashtable_length -= hashtable_length_correction;
  gettimeofday(&ragged_end, NULL);
  double ragged_time = ragged_end.tv_sec - ragged_start.tv_sec + (ragged_end.tv_usec - ragged_start.tv_usec) / 1.e6;

  fprintf(stderr, "\nRagged Array Creation Time: %g s\n", ragged_time);
  if (APPROXIMATE)
    fprintf(stderr, "\n**Approximation Calculation Time: %g s\n", approximate_time);

  fprintf(stderr, "UNIQUE BIPARTITIONS: %d\n", unique_bipartitions);
  NUMBIPART = unique_bipartitions;
  
//now print out numbiparts.txt if necessary (needed for multifurcation)
  if (!APPROXIMATE){
    if (tree_counter < ((NUM_TAXA-3)*NUM_TREES)) { 
      //fout.open("numbiparts.txt");
      //if (!fout) {
      //	cerr << "cannot open file for writing!\n";
      //}
      fprintf(stderr, "Trees are Multifurcating\n");
      //for (unsigned int i = 0; i < NUM_TREES; ++i) { 
      //	fout << i << "=" << n_biparts[i] << endl;
      //}
      if (!MULTIFURCATING){
	MULTIFURCATING = true;
	fprintf(stderr, "\nWARNING: Detected multifurcated trees!\n Trees will be treated as such through remainder of computation\n");
      }
      //fout.close();
    }
    else { 
      fprintf(stderr, "Trees are Binary\n");
      free(n_biparts);
    }
  }
  else{
    unsigned int first = n_biparts[0];
    bool change = false;
    for (unsigned int i = 0; i < NUM_TREES; i++){
      maxval_a[i] = n_biparts[i];
      if (n_biparts[i] != first)
	change = true;
    }
    if (change)
      MULTIFURCATING = true;
    else
      MAXVAL = first;
  }
  fprintf(stderr,"\nInverted Index Creation Time: %g s\n", index_time);

  // NUM_TREES = TOTAL_TREES - NUM_TREES2;
 //exit(0);  
}

void calculate(int a, int b, int d, int distance_measure, unsigned int &value, float &value2){
  int sig1, sig2;
  float temp, temp2;
  //cout << "----" << endl;
  value = 0;
  value2 = 0;
    switch (distance_measure) {
    case 0: 
    case 24: 
      value = b; //RF
      break;
    case 1:
      value = b+b; //Hamming
      break;
    case 2: 
      value2 =(float)a/(a+2*b); //Jaccard/Tanimoto
      break;
    case 3:
      value2 = (float)a/(a+b); //Dice
      break;
    case 4:
      value2 = (float)a/NUMBIPART; //Russel/Rao
      break;
    case 5:
      value2 = (float)a/(a+4*b); //Sokal/Sneath
      break;
    case 6:
      d = NUMBIPART - (a+b+b);
      value2 = (float)(a+d-2*b)/NUMBIPART; //Hamann
      break;
    case 7: 
      d = NUMBIPART - (a+b+b);
      value2 = (float)(a*d - b*b)/(a*d + b*b); //Yule
      break;
    case 8:
      value2 = (float)a/(sqrt((a+b)*(a+b))); //Ochai/Cosine
      break;
    case 9:
      d = NUMBIPART - (a+b+b);
      sig1 = max(a,b)+max(b,d)+max(a,b)+max(b,d);
      sig2 = max(a+b,b+d) + max(a+b,b+d);
      value2 = (float)(sig1-sig2)/(2*NUMBIPART); //Anderberg
      break;
    case 10: 
      value2 = (float)(NUMBIPART*a)/((a+b)*(a+b)); //Forbes
      break;
    case 11:
      temp = (float)(0.5*(a*b + a*b) + b*b); //Mountford
      if (temp == 0)
	value2 = NAN;
      else
	value2 = (float)a/(0.5*(a*b + a*b) + b*b); //Mountford
      break;
    case 12:
      value2 = (float)(pow(a,2))/((a+b)*(a+b)); //Sorgenfrei
      break;
    case 13:
      temp = (float)(a+b)/NUMBIPART;
      value2 = (float)log10(a) - log10(NUMBIPART) - log10(temp) - log10(temp); //Gilbert/Wells
      break;
    case 14:
      temp = (float) (NUMBIPART*a);
      temp = temp - ((a+b)*(a+b));
      temp2 =(float) (NUMBIPART*a + (a+b)*(a+b)); //Tarwid
      value2 = (float)temp/temp2;
      break;
    case 15:
      value2 = (float)a/(a+b) + (float)a/(a+b); //Johnson
      break;
    case 16:
      d = NUMBIPART - (a+b+b);
      temp = (float)(a*b + 2*b*b + b*d); //Peirce
      if (temp == 0)
	value2 = NAN;
      else
	value2 = (float)(a*b + b*b)/(a*b + 2*b*b + b*d); //Peirce
      break;
    case 17:
      temp = (float)1/(a+b);
      value2 = (float)(a/2)*(temp + temp); //driver/kroeber
      break;
    case 18:
      value2 = (float)a/(sqrt((a+b)*(a+b))) - (float)max(a+b,a+b)/2; //fager/mcgowan 
      break;
    case 19:
      value2 = (float)sqrt(b+b); //euclidean
      break;
    case 20:
      d = NUMBIPART - (a+b+b);
      temp2 = (float)abs(a*d - b*b);
      temp2 -= (float)(NUMBIPART/2);
      temp = (float)pow(temp2, 2);
      temp = (float)NUMBIPART*temp;
      temp2 = (a+b)*(a+b); 
      temp2 *= (b+d)*(b+d);
      temp2 = (float)temp/temp2;
      value2 = (float)log10(temp2); //Stiles
      break;
    case 21:
      d = NUMBIPART - (a+b+b);
      temp = a+d;
      value2 = (float)temp/(2*(b+b)+a+d); //Rogers/Tanimoto
      break;
    case 22:
      d = NUMBIPART - (a+b+b);
      temp2 = (float)NUMBIPART*a;
      temp2 -= ((a+b)*(a+b));
      temp = (float)(NUMBIPART*NUMBIPART)*temp2;
      temp2 =(float)(a+b)*(a+b)*(b+d)*(b+d);
      value2 = (float)temp/temp2; //Eyraud
      break;
    case 23:
      temp = b;
      temp2 = b+a;
      value2 = (float)temp/temp2; //weighted RF
      break;
    case 25:
      value2 = (float)(2*b)/NUMBIPART; //mean manhattan
      break;
    case 26:
      d = NUMBIPART - (a+b+b);
      value2 = (float)(a+d)/NUMBIPART; //simple matching
      break;
    case 27:
      value2 = (float)(a*d - b*b)/(sqrt((NUMBIPART*(a+b)*(a+b)))); //Dennis
      break;
    case 28:
      d = NUMBIPART - (a+b+b);
      sig1 = max(a,b)+max(b,d)+max(a,b)+max(b,d);
      sig2 = max(a+b,b+d) + max(a+b,b+d);
      value2 = (float)(sig1-sig2)/(2*NUMBIPART - sig2); //Goodman/Kruskal
      break;
    case 29:
      value2 = (float)(NUMBIPART*pow((a-0.5),2))/((a+b)*(a+b)); //Fossum
      break;
    case 30:
      value2 = (float)(b+b)/(2*a + b + b); //bray/curtis
      break;
    case 31:
      d = NUMBIPART - (a+b+b);
      temp = (float)(abs(a*(b+d)));
      temp2 = (float)(abs(b*(a+b)));
      if (temp2 == 0)
	value2 = NAN;
      else
	value2 = (float)temp/temp2; //AMPLE
      break;
    case 32:
      d = NUMBIPART - (a+b+b);
      temp = (float)(sqrt(a*d) - b);
      temp2 = (float)(sqrt(a*d) + b);
      value2 = (float)temp/temp2; //Yule-W
      break;
    case 33:
      value2 = (float)(a - (b+b))/NUMBIPART; //FAITH
      break;
    case 34:
      value2 = (float)(a - (b+b))/(a+b+b); //FAITH-2
      break;
    default:
      cerr << "Invalid option!" << endl;
      exit(0);
    }
}

void calculatem(int a, int b, int c, int d, int distance_measure, float & value){ //for multifurcating
  int sig1,sig2;
  float temp, temp2;

  //cout << "distance measure is: " << distance_measure << endl;
  //cout << "b is: " << b << endl;
  //cout << "c is: " << c << endl;
  //cout << "d is: " << d << endl;
  //cout << "a is: " << a << endl;
    switch (distance_measure) {
    case 0: 
    case 24: 
      value = (float)(b+c)/2; //RF
      //cout << "value is: " << value << endl;
      break;
    case 1: 
      value =  b+c; //Hamming
      break;
    case 2: 
      value =  (float)a/(a+b+c); //Jaccard/Tanimoto
      break;
    case 3: 
      value = (float)(2*a)/(2*a + b + c); //Dice
      break;
    case 4: 
      value = (float)a/NUMBIPART; //Russel/Rao
      break;
    case 5: 
      value = (float)a/(a + 2*b + 2*c); //Sokal/Sneath
      break; 
    case 6:
      d = NUMBIPART - (a+b+c);
      value = (float)(a+d-b-c)/NUMBIPART; //Hamann
      break;
    case 7:
      d = NUMBIPART - (a+b+c);
      value = (float)(a*d - b*c)/(a*d + b*c); //Yule
      break;
    case 8:
      value = (float)a/(sqrt((a+b)*(a+c))); //Ochai/Cosine
      break;
    case 9:
      d = NUMBIPART - (a+b+c);
      sig1 = max(a,b)+max(c,d)+max(a,c)+max(b,d);
      sig2 = max(a+c,b+d) + max(a+b,c+d);     
      value = (float)(sig1-sig2)/(2*NUMBIPART); //Anderberg
      break;
    case 10: 
      value = (float)(NUMBIPART*a)/((a+b)*(a+c)); //Forbes
      break;
    case 11:
      temp = (float)(0.5*(a*b + a*c)) + b*c;
      if (temp == 0)
	value = NAN;
      else
      value = (float)a/(0.5*(a*b + a*c) + b*c); //Mountford
      break;
    case 12:
      value = (float)(pow(a,2))/((a+b)*(a+c)); //Sorgenfrei
      break;
    case 13:
      temp = (float)(a+b)/NUMBIPART;
      temp2 = (float)(a+c)/NUMBIPART;
      value = (float)log10(a) - log10(NUMBIPART) - log10(temp) - log10(temp2); //gilbert/wells
      break;
    case 14:
      temp = (float)(NUMBIPART*a);
      temp = temp - ((a+b)*(a+c));
      value = (NUMBIPART*a + (a+b)*(a+c)); //Tarwid
      break;
    case 15:
      value = (float)(a/(a+b)) + (float)(a/(a+c)); //Johnson
      break;
    case 16:
      d = NUMBIPART - (a+b+c);
      temp = (float)(a*b + 2*b*c + c*d);
      if (temp == 0)
	value = NAN;
      else
	value = (float)(a*b + b*c)/(a*b + 2*b*c + c*d); //Peirce
      break;
    case 17:
      temp = (float)1/(a+b);
      temp2 = (float)1/(a+c);
      value = (float)(a/2)*(temp + temp2); //Driver/Kroeber
      break;
    case 18:
      value = (float)a/(sqrt((a+b)*(a+c))) - (float)max(a+b,a+c)/2; //Fager/McGowan 
      break;
    case 19:
      value = (float)sqrt(b+c); //Euclid
      break;
    case 20:
      d = NUMBIPART - (a+b+c);
      temp = (float)NUMBIPART*pow((abs(a*d - b*c)-(NUMBIPART/2)),2);
      temp2 = (a+b)*(a+c);
      temp2*= (b+d)*(c+d);
      temp2 = (float)temp/temp2;
      value = (float)log10(temp2); //Stiles
      break;
    case 21:
      d = NUMBIPART - (a+b+c);
      temp = a+d;
      value = (float)temp/(2*(b+c)+a+d); //Rogers/Tanimoto
      break;
    case 22:
      d = NUMBIPART - (a+b+c);
      temp2 = (float)NUMBIPART*a;
      temp2 -= (float)((a+b)*(a+c));
      temp = (float)(NUMBIPART*NUMBIPART)*temp2;
      temp2 =(float)(a+b)*(a+c)*(b+d)*(c+d);
      value = (float)temp/temp2; //Eyraud
      break;
    case 23:
      temp = (float)(b+c);
      temp2 = (float)(2*a + b + c);
      value = (float)temp/temp2; //RF rate
      break;
    case 25:  
      value = (float)(b+c)/NUMBIPART; //Mean Manhattan
      break;
    case 26:
      d = NUMBIPART - (a+b+c);
      value = (float)(a+d)/NUMBIPART; //Simple Matching
      break;
    case 27:
      d = NUMBIPART - (a+b+c);
      value = (float)(a*d - b*c)/(sqrt((NUMBIPART*(a+b)*(a+c)))); //Dennis
      break;
    case 28:
      d = NUMBIPART - (a+b+c);
      sig1 = max(a,b)+max(c,d)+max(a,c)+max(b,d);
      sig2 = max(a+c,b+d) + max(a+b,c+d);     
      value = (float)(sig1-sig2)/(2*NUMBIPART - sig2); //Goodman/Kruskal
      break;
    case 29:
      value = (float)(NUMBIPART*pow((a-0.5),2))/((a+b)*(a+c)); //Fossum
      break;
    case 30:
      value = (float)(b+c)/(2*a + b + c); //Bray/Curtis
      break;
    case 31:
      d = NUMBIPART - (a+b+c);
      temp = (float)(abs(a*(c+d)));
      temp2 = (float)(abs(c*(a+b)));
      if (temp2 == 0)
	value = NAN;
      else
	value = (float)temp/temp2; //AMPLE
      break;
    case 32:
      d = NUMBIPART - (a+b+c);
      temp = (float)(sqrt(a*d) - sqrt(b*c));
      temp2 = (float)(sqrt(a*d) + sqrt(b*c));
      value = (float)temp/temp2; //Yule-W
      break;
    case 33:
      value = (float)(a - (b+c))/NUMBIPART; //FAITH
      break;
    case 34:
      value = (float)(a - (b+c))/(a+b+c); //FAITH-2
      break;
    default: 
      cerr << "Invalid option!" << endl;
      exit(0);
    }
}

void by_row_distance(unsigned int hash_length, unsigned int * row, unsigned int distance_measure, string outfilename) { 

  if (PRINT || PRINT_L)
    fprintf(stderr, "outfilename is: %s\n", outfilename.c_str()); 
  double traversal_time = 0;
  //struct timeval print_start;
  //struct timeval print_end;
  struct timeval traversal_start;
  struct timeval traversal_end;
  struct timeval calculate_start;
  struct timeval calculate_end;
 struct timeval print_start;
  struct timeval print_end;
  double print_temp;
  unsigned int per_10, per_counter, per_place;
  double sum = 0;
  int count = 0;
  double average = 0;

  per_10 = NUM_TREES/10;
  per_counter = 0;
  per_place = 0;

  gettimeofday(&traversal_start, NULL);
  gettimeofday(&calculate_start, NULL);  
  //intialize RFMatrix row array
  unsigned int * RFMatrixRow = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  if (RFMatrixRow == NULL){
    fprintf(stderr, "Cannot allocate RFMatrixRow: Out of Memory!\n");
    exit(3);
  }
  long num_comparisons = 0;
  
  /*unsigned int * bipart_counts = (unsigned int*)malloc(NUM_TREES*sizeof(unsigned int));
  if (bipart_counts == NULL){
    fprintf(stderr, "Cannot allocate bipart_counts: Out of Memory!\n");
    exit(3);
    }*/

  /*if (MULTIFURCATING && !APPROXIMATE) {
    ifstream fin;
    fin.open("numbiparts.txt");
    if (!fin) { 
      cerr << "Cannot open file!" << endl;
      return;
    }
    string line;
    unsigned int array_pos = 0;
    int pos = 0;
    unsigned int mycount = 0;
    while (!fin.eof()) { //read bipartition counts in from file
      fin >> line;
      pos = line.find_first_of("=");
      line = line.substr(pos+1);
      mycount = atoi(line.c_str());
      bipart_counts[array_pos] = mycount;
      ++array_pos;
    }
    }*/

  if (MULTIFURCATING && !APPROXIMATE)
    assert(n_biparts != NULL);

  unsigned int hash_loc = 0;
  unsigned int temp = 0;
  fprintf(stderr, "\nBeginning By-Tree Traversal:\n");
  ofstream fout;
  if (PRINT || PRINT_L){
    fout.open(outfilename.c_str());
    if (!fout){
      cerr << "ERROR! Cannot open file for writing!\n";
      exit(3);
    }
  }

  if (!QOPT)
    TOTAL_TREES = NUM_TREES;
  else
    NUM_TREES = TOTAL_TREES - NUM_TREES2;

 if (TRIVIAL)
    NUMBIPART += NUM_TAXA;

  for (unsigned int pos =0; pos < NUM_TREES; ++pos) { //for every row 
    //fprintf(stderr,"\npos is: %d\n", pos);
    for (unsigned int x = 0; x < TOTAL_TREES; ++x) { //initialize 
      RFMatrixRow[x] = 0;
    }

    //calculate similarity (value: 'a' from Salim et. al paper)
    for (unsigned int i = 0; i < helper_sizes[pos]; ++i) {    
      hash_loc = helper[pos][i]; //get id
      for (unsigned int j = 1; j < hash_lengths[hash_loc]; ++j) { //hash line
	temp = hashtable[hash_loc][j];
	if(temp!=pos){
	  ++num_comparisons;
	  ++RFMatrixRow[temp];
	}
      }
    }

 
  
    unsigned int maxvalue = MAXVAL; 
    int a,b,c,d;
    unsigned int value;
    float value2;
    stringstream rs;
    unsigned int startpos;
    unsigned int endpos;
    if (QOPT){
      startpos = NUM_TREES;
      endpos = TOTAL_TREES;
    }
    else{
      startpos = 0;
      endpos = NUM_TREES;
    }


    if (FULLMATRIX){
      if (MULTIFURCATING) {
	d = 0;
	for (unsigned int y = startpos; y < endpos; ++y) {
	  if (y!= pos) {  
	    a = ALL_COUNT - row[pos] - row[y] + RFMatrixRow[y];
	    if (TRIVIAL)
	      a += NUM_TAXA; 
	    if (APPROXIMATE){
	      if (TRIVIAL){
		b = maxval_a[pos] + NUM_TAXA -a; //b
		c = maxval_a[y] + NUM_TAXA - a; //c
	      }
	      else{
		b = maxval_a[pos] -a; //b
		c = maxval_a[y] - a; //c
	      }
	    }
	    else{
	      //cout << "a is: " << a << endl;
	      //cout << "bipart_count[pos]: " << bipart_counts[pos] << endl;
	      //cout << "bipart_count[y]: " << bipart_counts[y] << endl;
	      if (TRIVIAL){
		b = n_biparts[pos] + NUM_TAXA- a; //b 
		c = n_biparts[y] + NUM_TAXA - a; //c
	      }
	      else{
		b = n_biparts[pos] - a; //b 
		c = n_biparts[y] - a; //c
	      }
	      //cout << "b is: " << b << endl;
	      //cout << "c is: " << b << endl;
	    }

	  }
	  else{
	    a = maxvalue;
	    if (TRIVIAL)
	      a += NUM_TAXA;
	    b = 0;
	    c = 0;
	  }
	  calculatem(a,b,c,d,distance_measure,value2);
	  //cout << "value2 is: " << value2 << endl;
	  if (PRINT)
	    rs << value2 << " ";
	  else if (PRINT_L)
	    rs << pos << "." << y << " " << value2 << "\n";
	  else{
	    if (value2 != NAN){
	      sum += value2;
	      count++;
	    }
	  }
	}
	if (PRINT || PRINT_L){
	  gettimeofday(&print_start, NULL);
	  fout << rs.str();
	  if (PRINT)
	    fout << endl;
	  gettimeofday(&print_end, NULL);
	  print_temp = print_end.tv_sec - print_start.tv_sec + (print_end.tv_usec - print_start.tv_usec) / 1.e6; 
	  print_time += print_temp;
	}
	rs.str("");
      } //end multifurcating
      else {
	d = 0;
	for (unsigned int y = startpos; y < endpos; ++y) { 
	  if (y != pos) { 
	    a = ALL_COUNT - row[pos] - row[y] + RFMatrixRow[y];
	    if (TRIVIAL){
	      a += NUM_TAXA; 
	      b = maxvalue + NUM_TAXA - a;
	    }
	    else{
	      b = maxvalue - a;
	    }
	  }
	  else{
	    a = maxvalue;
	    if (TRIVIAL)
	      a+= NUM_TAXA;
	    b = 0;
	  }
	  calculate(a,b,d,distance_measure, value, value2);
	  if (!PRINT_L){
	    if (PRINT){
	      if ((distance_measure < 2) || (distance_measure == 24))
		rs << value << " ";
	      else
		rs << value2 << " ";
	    }
	    else{
	   
	      if ((distance_measure < 2) || (distance_measure == 24)){
		if (value != NAN){
		  sum += value;
		  count++;
		}
	      }
	      else{
		if (value2 !=  NAN){
		  sum += value2;
		  count++;
		}
	      }

	      //cout << "we get here!" << endl;
	      //cout << "sum is: " << sum << endl;
	      //cout << "count is: " << count << endl;
	      //exit(0);
	    }
	  }
	  else{
	    if ((distance_measure < 2) || (distance_measure == 24))
	      rs << pos << "." << y << " " << value << "\n";
	    else
	      rs << pos << "." << y << " " << value2 << "\n";
	  }
	}
	if (PRINT || PRINT_L){
	  gettimeofday(&print_start, NULL);
	  fout << rs.str();
	  if (PRINT)
	    fout << endl;
	  gettimeofday(&print_end, NULL);
	  print_temp = print_end.tv_sec - print_start.tv_sec + (print_end.tv_usec - print_start.tv_usec) / 1.e6; 
	  print_time += print_temp;
	}
	rs.str("");
      }
    }
    else{ //half matrix
      if (MULTIFURCATING) {
	d = 0;
	for (unsigned int y = startpos; y < startpos+pos; ++y) {
	  a = ALL_COUNT - row[pos] - row[y] + RFMatrixRow[y];
	  if (TRIVIAL)
	    a += NUM_TAXA; 
	  if (APPROXIMATE){
	    if (TRIVIAL){
	      b = maxval_a[pos] + NUM_TAXA -a; //b
	      c = maxval_a[y] + NUM_TAXA - a; //c
	    }
	    else{
	      b = maxval_a[pos] -a; //b
	      c = maxval_a[y] - a; //c
	    }
	  }
	  else{
	    if (TRIVIAL){
	      b = n_biparts[pos] + NUM_TAXA - a; //b 
	      c = n_biparts[y] + NUM_TAXA - a; //c
	    }
	    else{
	      b = n_biparts[pos] - a; //b 
	      c = n_biparts[y] - a; //c
	    }
	  }
	  calculatem(a,b,c,d,distance_measure, value2);
	  if (PRINT)
	    rs << value2 << " ";
	  else if (PRINT_L)
	    rs << pos << "." << y << " " << value2 << "\n";
	  else {
	    if (value2 != NAN){
	      count++;
	      sum += value2;
	    }
	  }
	}
	if (PRINT || PRINT_L){
	  gettimeofday(&print_start, NULL);
	  fout << rs.str();
	  if (PRINT)
	    fout << endl;
	  gettimeofday(&print_end, NULL);
	  print_temp = print_end.tv_sec - print_start.tv_sec + (print_end.tv_usec - print_start.tv_usec) / 1.e6; 
	  print_time += print_temp;
	}
	rs.str("");
      } //end multifurcating
      else {
	d = 0;
	for (unsigned int y = startpos; y < startpos+pos; ++y) { 
	  a = ALL_COUNT - row[pos] - row[y] + RFMatrixRow[y];
	  if (TRIVIAL){
	    a += NUM_TAXA; 
	    b = maxvalue + NUM_TAXA - a;
	  }
	  else{
	    b = maxvalue - a;
	  }
	  calculate(a,b,d,distance_measure, value, value2);
	  if (!PRINT_L){
	    if (PRINT){
	      if ( (distance_measure < 2) || (distance_measure == 24))
		rs << value << " ";
	      else
		rs << value2 << " ";
	    }
	    else{
	      if ( (distance_measure < 2) || (distance_measure == 24)){
		if (value != NAN){
		  sum += value;
		  count++;
		}
	      }
	      else{
		if (value2 != NAN){
		  sum += value2;
		  count++;
		}
	      }
	    }
	  }
	  else{
	    if ((distance_measure < 2) || (distance_measure == 24))
	      rs << pos << "." << y << " " << value << "\n";
	    else
	      rs << pos << "." << y << " " << value2 << "\n";
	  }
	}
	if (PRINT || PRINT_L){
	  gettimeofday(&print_start, NULL);
	  fout << rs.str();
	  if (PRINT)
	    fout << endl;
	  gettimeofday(&print_end, NULL);
	  print_temp = print_end.tv_sec - print_start.tv_sec + (print_end.tv_usec - print_start.tv_usec) / 1.e6; 
	  print_time += print_temp;
	}
	rs.str("");
      }
    }

    //cout << "we get here?" << endl;
    //calculate printing times
    if (NUM_TREES > 100){
      if (pos % per_10 == 0){
	per_place = per_counter*10;
	fprintf(stderr, "%d%% Complete\n", per_place);
	per_counter++;
      }
    }
  }
  if (PRINT|| PRINT_L){
    //fout << "Number of comparisons: " <<  num_comparisons << endl;
    fout.close();
  }
  average = (double)sum/count;
  //cout << "sum is: " << sum << endl;
  //cout << "count is: " << count << endl;
  cout << average << endl;

  fprintf(stderr, "Finishing Traversal\n");  
  fprintf(stderr, "Number of comparisons: %u\n", num_comparisons);
  gettimeofday(&traversal_end, NULL);
  gettimeofday(&calculate_end, NULL);  
  traversal_time = traversal_end.tv_sec - traversal_start.tv_sec + (traversal_end.tv_usec - traversal_start.tv_usec) / 1.e6; 
  double calculate_time = calculate_end.tv_sec - calculate_start.tv_sec + (calculate_end.tv_usec - calculate_start.tv_usec) / 1.e6; 
  traversal_time -= print_time;
  calculate_time -= print_time;
  fprintf(stderr, "\nTraversal Time: %g s\n", traversal_time); 
  fprintf(stderr, "\nRF Calculation Time: %g s\n", calculate_time); 
  fprintf(stderr, "\nMatrix Print Time: %g s\n", print_time); 
  free(RFMatrixRow);
  if (TRIVIAL)
    NUMBIPART -= NUM_TAXA;
  //free(bipart_counts);
}


/**********************
 *
 * main function
 *
 *********************/


int main(int argc, char** argv)
{
  string outfilename = "phlash.out";
  string qfilename;
  string infilename;
  int distance_option = 0;
  
  int approximate_threshold = 0;
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
    string helpMsg  = "phlash <input file> <number of trees> -d[distance measure] -s[seed value] -m\n";

    helpMsg += "Input file: \n";
    helpMsg += "   The current version of Phlash only supports the Newick format.\n";

    helpMsg += "Example of Newick tree: \n";
    helpMsg += "   (('Chimp':0.052625,'Human':0.042375):0.007875,'Gorilla':0.060125,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
    helpMsg += "   ('Chimp':0.052625,('Human':0.042375,'Gorilla':0.060125):0.007875,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
    helpMsg += "Specify seed value: \n";   
    helpMsg += "   -s <value>, specify seed value (default: 1000) \n";
    helpMsg += "   -m specify whether to print out the matrix. Default: matrix is computed, not printed\n";

    //helpMsg += " The seed file is necessary for parallel processing. A common seed file is necessary to ensure labels are properly populated\n"; 
    helpMsg += " -d <distance measure>: 0..6. 0 is Robinson-Foulds distance, and is used as default.\n"; 
    helpMsg += " -a <approximate percentage>: 0..100. Do approximate bipartition matching, match this percentage of taxa. If you want to specify mismatches, be sure to also turn on the -t flag\n\n";
    helpMsg += "****List of implemented Distance Measures****\n";
    helpMsg += " 1: Hamming             12: Sorgeinfrie     23: RF Rate/Lance-Williams\n";
    helpMsg += " 2: Jaccard-Tanimoto    13: Gilbert-Wells   24: RF\n";
    helpMsg += " 3: Dice                14: Tarwid          25: Mean Manhattan\n";
    helpMsg += " 4: Russel-Rao          15: Johnson         26: Simple Matching\n";
    helpMsg += " 5: Sokal-Sneath        16: Peirce          27: Dennis\n";
    helpMsg += " 6: Hammann             17: Driver-Kroeber  28: Goodman-Kruskal\n";
    helpMsg += " 7: Yule                18: Fager-McGowan   29: Fossum\n";
    helpMsg += " 8: Ochai-Cosine        19: Euclidean       30: Bray-Curtis\n";
    helpMsg += " 9: Anderberg           20: Stiles          31: Ample\n";
    helpMsg += " 10: Forbes             21: Rogers-Tanimoto 32: Yule-W\n";
    helpMsg += " 11: Mountford          22: Eyraud          33: Faith    \n";
    helpMsg += " 34. Faith-2\n";
    helpMsg += "*************************************************\n\n";
    helpMsg += "Examples: \n";
    helpMsg += " phlash  foo.tre 1000 -d 0 //input is foo.tre, calculating an RF matrix\n";
    helpMsg += " phlash foo.tre 1000 -d 0 -m //input is foo.tre, calculating and printing RF matrix\n";
    helpMsg += "phlash foo.tre 1000 -d 2 -s 17 //input is foo.tre, caculating based on Jaccard-Tanimoto Similarity, seed set to 17 (matrix not printed)\n";
    helpMsg += "phlash foo.tre 1000 -d 1 -a 80 //input is foo.tre, caculating based on Hamming Similarity, approximate bipartition matching of 80%\n";
    helpMsg += "phlash foo.tre 1000 -d 1 -a 2 -t  //input is foo.tre, caculating based on Hamming Similarity, approximate bipartition based on 2 taxa mismatches (NUM_TAXA-2)\n";    
    TCLAP::CmdLine cmd(helpMsg, ' ', "1.0.0");

    TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "file name", true, "intree", "Input tree file name"  );
    cmd.add( fnameArg );

    TCLAP::ValueArg<string>  qnameArg( "q", "qfile", "q file", false, "", "optional second tree file name"  );
    cmd.add( qnameArg );

    TCLAP::UnlabeledValueArg<int>  numtreeArg( "numtree", "number of trees", true, 2, "Number of trees"  );
    cmd.add( numtreeArg );

    TCLAP::ValueArg<int>  aArg( "a", "approximate", "approximate distance", false, 0, "approximate distance threshold"  );
    cmd.add( aArg );

    TCLAP::ValueArg<unsigned int> cArg("c", "cvalue", "c value", false, 1000, "c value");
    cmd.add( cArg );

    TCLAP::ValueArg<string> oArg("o", "output", "output file", false, "phlash.out", "output file name");
    cmd.add( oArg );

    TCLAP::SwitchArg mArg("m", "matrix", "print matrix?", false);
    cmd.add( mArg );

    TCLAP::SwitchArg fArg("f", "fullmatrix", "calculate full matrix?", false);
    cmd.add( fArg );

    TCLAP::SwitchArg tArg("t", "taxamismatch", "use number of taxa mismatches for computation instead", false);
    cmd.add( tArg );

    TCLAP::SwitchArg lArg("l", "list", "print matrix as list?", false);
    cmd.add( lArg );
    
    //r changed from "rate" to "rooted" -- if we ever want to use rate, use "n" for normalized
    TCLAP::SwitchArg rArg("r", "rooted", "treat the input trees as rooted?", false);
    cmd.add( rArg );

    TCLAP::SwitchArg nArg("n", "nocheck", "no check against bitstrings (montecarlo)", false);
    cmd.add( nArg );

    TCLAP::ValueArg<int> seedArg("s", "seedvalue", "user specified seed value", false, 1000, "user specified seed value");
    cmd.add( seedArg );

    TCLAP::ValueArg<int> dArg("d", "distance", "user specified distance measure", true, 0, "user specified distance measure");
    cmd.add( dArg );

    TCLAP::SwitchArg bArg("b", "bipart", "use trivial bipartitions?", false);
    cmd.add( bArg );



    cmd.parse( argc, argv );

    NUM_TREES = numtreeArg.getValue();
    //NUM_TAXA = numtaxaArg.getValue();
    int toprint = mArg.getValue();
    int toprint_list = lArg.getValue();
    if (aArg.getValue()){
      //cout << "we get here?? aArg is:" << aArg.getValue() << endl;
      APPROXIMATE = true;
      approximate_threshold  =aArg.getValue();
    }
    if (tArg.getValue()){
      MISMATCH = true;
    }

    if (rArg.getValue())
      ROOTED = true;

    if (nArg.getValue())
      MC = true;

    if (fArg.getValue())
      FULLMATRIX = true;

    if (bArg.getValue())
      TRIVIAL = true;

    if (strcmp(oArg.getValue().c_str(), "phlash.out")!= 0)
      outfilename = oArg.getValue();

    //int toprinthash = bArg.getValue();
    distance_option = dArg.getValue();
    //int use_rate = rArg.getValue();
    if (toprint && toprint_list) { 
      cerr << "Cannot specify both print as list and print as matrix! Please choose one or the other" << endl;
      return 2;
    }
  
    if (toprint) 
      PRINT = true;

    if (toprint_list)
      PRINT_L = true;
    /*if (toprinthash) 
    PRINT_HASH = true;

    if (data_option) 
      DATA = true;

    if (use_rate) 
      RATE = true;
    */
    infilename = fnameArg.getValue();

    qfilename = qnameArg.getValue();
    if (qfilename != "")
      QOPT = true;

    if (NUM_TREES == 0) {
      string strFileLine;
      unsigned int ulLineCount;
      ulLineCount = 0;

      ifstream infile(infilename.c_str());
      if (infile) {
        while (getline(infile, strFileLine)) {
          ulLineCount++;
        }
      }

      fprintf(stderr, "Number of trees in the input file: %u\n", ulLineCount);
      NUM_TREES = ulLineCount;
      infile.close();
    }

    if  (QOPT){
      string strFileLine;
      unsigned int ulLineCount;
      ulLineCount = 0;

      ifstream infile(qfilename.c_str());
      if (infile) {
        while (getline(infile, strFileLine)) {
          ulLineCount++;
        }
      }
      fprintf(stderr, "Number of trees in the second input file: %u\n", ulLineCount);
      NUM_TREES2 = ulLineCount;
      infile.close();
      TOTAL_TREES = NUM_TREES + NUM_TREES2;
      fprintf(stderr, "Total number of trees: %u\n", TOTAL_TREES);
    }

    if (NUM_TREES < 2) {cerr << "Fatal error: at least two trees expected.\n"; exit(2);}

    if (cArg.getValue())
      C = cArg.getValue();
    
    if (seedArg.getValue()) { 
      NEWSEED = seedArg.getValue();
    }

    /*if (uArg.getValue()){
      UNIQUE = 1;
      AMAZING = true;
      P = 0;
      }*/
  } catch (TCLAP::ArgException &e) { // catch any exceptions
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  if (QOPT)
    cerr << "q file name is: " << qfilename << endl;
  else
    cerr << "q option not set!" << endl;

  gettimeofday(&process_start, NULL);
  unsigned int hashtable_length = 0; //stores the actual length of the hashtable 
  unsigned int * row; //stores the amount to be subtracted from rows at the end
  unsigned int rowsize;
  if (QOPT)
    rowsize = TOTAL_TREES;
  else 
    rowsize = NUM_TREES;

  row = (unsigned int *)malloc(rowsize * sizeof(unsigned int));
  //unsigned int threshold = NUM_TREES - (NUM_TREES*P)/100;
  //threshold--;
  for (unsigned int i = 0; i < rowsize; i++) { 
    row[i] = 0;
  }   
  if (MC)
    cout << "NOTE: Monte-Carlo Option on!" << endl;
  process_trees(infilename, qfilename, hashtable_length, row, approximate_threshold);
  gettimeofday(&process_end, NULL);
  
  double process_time = process_end.tv_sec - process_start.tv_sec + (process_end.tv_usec - process_start.tv_usec) / 1.e6;
  fprintf(stderr,"\nInput Processing Time: %g s\n", process_time);
   
     //debugging
   /*   
  cout << "hashtable" << endl;
  cout << "hashtable length is: " << hashtable_length << endl;
  cout << "hash_lengths: " << endl;
  for (unsigned int i = 0; i < hashtable_length; ++i) { 
    cout << hash_lengths[i] << " ";
  }
  cout << endl;

  for (unsigned int i = 0; i < hashtable_length; ++i) { 
    cout << "H[" << i << "]: next";
    for (unsigned int j = 0; j < hash_lengths[i]; ++j) { 
      cout << hashtable[i][j] << " ";
    }
    cout << endl;
  } 
  cout << "row array: " << endl;
  for (unsigned int i = 0; i < NUM_TREES; ++i) { 
    cout << row[i] << " ";
  }
  cout << endl;
  cout << "ALL_COUNT: " << ALL_COUNT << endl;
  
  cout << "NUMBER OF UNIQUE BIPARTITIONS: " << NUMBIPART << endl;
  cout << "helper data structure: " << endl;
  for (unsigned int i = 0; i < NUM_TREES; ++i) {
    cout << "Tree " << i << ": ";
    for (unsigned int j = 0; j < helper_sizes[i]; ++j) {
      cout << helper[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
  
  cout << "full index data structure: " << endl;
  for (unsigned int i = 0; i < NUM_TREES; ++i) {
    cout << "Tree " << i << ": ";
    for (unsigned int j = 0; j < full_index[i].size(); ++j) {
      cout << full_index[i][j] << " ";
    }
    cout << endl;
  }
exit(0);*/
//end debugging


  //gettimeofday(&traverse_start, NULL);
  fprintf(stderr, "distance measure to be used is:  %d\n", distance_option);
  by_row_distance(hashtable_length, row, distance_option, outfilename);
  gettimeofday(&prog_end, NULL);
  double prog_time = prog_end.tv_sec - prog_start.tv_sec + (prog_end.tv_usec - prog_start.tv_usec) / 1.e6;
  fprintf(stderr,"\nProgram Running Time: %g s\n", prog_time);
  //fprintf(stderr, "do we get here?\n");
  return 0;
}

