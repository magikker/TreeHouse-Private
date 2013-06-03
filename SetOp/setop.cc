//*****************************************************************/
/*
This is SetOp, a program that does fast hashing set operations.
It is part of the HaB|TAT suite, written by Suzanne J. Matthews

(c) 2010 SetOp: Suzanne J. Matthews 

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

#include "label-map.hh"
#include "hashfunc.hh"
#include "hash.hh"
#include "./tclap/CmdLine.h"
#include "global.h"
#include "parsing.h"
#include "buildtree.h"

// For newick parser
extern "C" {
#include <newick.h>
}


using namespace std;


/**********************
 *
 * main function
 *
 *********************/

vector<unsigned int> calculate_union(unsigned int NUM_TREES2, HashRFMap & vvec_setop){
  fprintf(stderr, "calculating union!\n");
  vector<unsigned int> union_ids;
  unsigned int union_count = 0;
  for (unsigned int hti=0; hti<vvec_setop._hashtab2.size(); ++hti) {
    unsigned int sizeVec = vvec_setop._hashtab2[hti].size(); 
    if (sizeVec) {
      for (unsigned int i=0; i<sizeVec; ++i) {
	union_count++;       
	union_ids.push_back(vvec_setop._hashtab2[hti][i]._vec_treeidx[0]);
      }
    } //end if sizevec
  } //end hashtable traversal
  sort(union_ids.begin(), union_ids.end());
  fprintf(stderr, "Found %u non-redundant trees\n", union_count);
  return union_ids;
}

vector<unsigned int> calculate_intersection(unsigned int NUM_TREES2, HashRFMap & vvec_setop){
  fprintf(stderr, "calculating intersection!\n");
  vector<unsigned int> intersection_ids;
  unsigned int intersection_count = 0;
  bool first, second;
  unsigned int temp_id;
  for (unsigned int hti=0; hti<vvec_setop._hashtab2.size(); hti++) {
    unsigned int sizeVec = vvec_setop._hashtab2[hti].size(); 
    if (sizeVec) {
      for (unsigned int i=0; i<sizeVec; i++) {
	unsigned int sizeTreeIdx = vvec_setop._hashtab2[hti][i]._vec_treeidx.size();
	first = false;
	second = false;
	for (unsigned int j = 0; j < sizeTreeIdx; j++){  
	  temp_id = vvec_setop._hashtab2[hti][i]._vec_treeidx[j];
	  if (temp_id < NUM_TREES)
	    first = true;
	  if (temp_id >= NUM_TREES)
	    second = true;
	}
	if (first && second){
	  intersection_ids.push_back(vvec_setop._hashtab2[hti][i]._vec_treeidx[0]);
	  intersection_count++;
	}
      }
    } //end if sizevec
  } //end hashtable traversal
  sort(intersection_ids.begin(), intersection_ids.end());
  fprintf(stderr, "Found %u non-redundant trees in intersection\n", intersection_count);
  return intersection_ids;
}

vector<unsigned int> calculate_set_difference(unsigned int NUM_TREES2, HashRFMap & vvec_setop){

  fprintf(stderr, "calculating set difference!\n");
  vector<unsigned int> difference_ids;
  unsigned int difference_count = 0;
  bool first, second;
  unsigned int temp_id;
  for (unsigned int hti=0; hti<vvec_setop._hashtab2.size(); hti++) {
    unsigned int sizeVec = vvec_setop._hashtab2[hti].size(); 
    if (sizeVec) {
      for (unsigned int i=0; i<sizeVec; i++) {
	unsigned int sizeTreeIdx = vvec_setop._hashtab2[hti][i]._vec_treeidx.size();
	first = false;
	second = false;
	for (unsigned int j = 0; j < sizeTreeIdx; j++){  
	  temp_id = vvec_setop._hashtab2[hti][i]._vec_treeidx[j];
	  if (temp_id < NUM_TREES)
	    first = true;
	  if (temp_id >= NUM_TREES)
	    second = true;
	}
	if ( (first && !second) ){
	  difference_ids.push_back(vvec_setop._hashtab2[hti][i]._vec_treeidx[0]);
	  difference_count++;
	}
      }
    } //end if sizevec
  } //end hashtable traversal
  sort(difference_ids.begin(), difference_ids.end());
  fprintf(stderr, "Found %u non-redundant trees in set difference\n", difference_count);
  return difference_ids;
}

int main(int argc, char** argv)
{
  string outfilename = "setop.out";
  string infilename, infilename2;
  unsigned int SET_OP = 0;
  unsigned int NUM_TREES2 = 0;
  vector<string> treefile;
  vector<string> treefile2;

  // TCLAP
  //since we are creating the hashtable serially, the assumption in this
  //version is that the seed file is the same as the input file and that 
  //offset is set to 0
  //this can readily be changed in future version if necessary
  try {

    // Define the command line object.
    string helpMsg  = "setop <input file> <input file2> -s[seed value] -o [output file] [set operation]\n";

    helpMsg += "Input file: \n";
    helpMsg += "   The current version of SetOp only supports the Newick format.\n";

    helpMsg += "Example of Newick tree: \n";
    helpMsg += "   (('Chimp':0.052625,'Human':0.042375):0.007875,'Gorilla':0.060125,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
    helpMsg += "   ('Chimp':0.052625,('Human':0.042375,'Gorilla':0.060125):0.007875,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
    helpMsg += "Specify seed value: \n";   
    helpMsg += "   -s <value>, specify seed value (default: 1000) \n";
    helpMsg += "   -o <output file>, specifies the name of the output file. Default: mybprof.out \n";

    helpMsg += "Examples: \n";
    helpMsg += "  setop foo.tre foo2.tre -u //union of foo.tre and foo2.tre\n";
    helpMsg += "  setop foo.tre foo2.tre -i //intersection of foo.tre and foo2.tre\n";
    helpMsg += "  setop foo.tre foo2.tre -d //set difference of foo.tre and foo2.tre\n";
    helpMsg += "  setop foo.tre foo2.tre -u -o afile.txt -s 17 //input is foo.tre and foo2.tre, seed is set to 17. Find the union and place the output in afile.txt\n";
    
    TCLAP::CmdLine cmd(helpMsg, ' ', "1.0.0");

    TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "file name", true, "intree", "input file"  );
    cmd.add( fnameArg );

    TCLAP::UnlabeledValueArg<string>  fnameArg2( "name2", "file name2", true, "intree2", "input file2"  );
    cmd.add( fnameArg2 );

    TCLAP::SwitchArg uArg("u", "union", "calculate union?", false);
    cmd.add( uArg );

    TCLAP::SwitchArg iArg("i", "intersection", "calculate intersection?", false);
    cmd.add( iArg );
    TCLAP::SwitchArg dArg("d", "difference", "calculate set difference?", false);
    cmd.add( dArg );

    TCLAP::ValueArg<string>  outputArg( "o", "output", "output file name", false, "setop.out", "user specified output file name"  );
    cmd.add( outputArg );

    TCLAP::ValueArg<int> seedArg("s", "seedvalue", "user specified seed value", false, 1000, "user specified seed value");
    cmd.add( seedArg );

    cmd.parse( argc, argv );

    NUM_TAXA = 0;
    NUM_TREES = 0;

    infilename = fnameArg.getValue();
    infilename2 = fnameArg2.getValue();

    outfilename = outputArg.getValue();

    if(uArg.getValue())
      SET_OP = 1;

    if (iArg.getValue())
      SET_OP = 2;

    if (dArg.getValue())
      SET_OP = 3;

    if (seedArg.getValue()) { 
      NEWSEED = seedArg.getValue();
    }

    if (NUM_TREES == 0) {
      string strFileLine;
      unsigned int ulLineCount;
      ulLineCount = 0;
      
      ifstream infile(infilename.c_str());
      if (infile) {
	while (getline(infile, strFileLine)) {
	  treefile.push_back(strFileLine);
	  ulLineCount++;
      }
      }
      fprintf(stderr, "Number of trees in the input file: %u\n", ulLineCount);
      NUM_TREES = ulLineCount;
      
      infile.close();
    }
    
    if (NUM_TREES2 == 0) {
      string strFileLine;
      unsigned int ulLineCount;
      ulLineCount = 0;

      ifstream infile(infilename2.c_str());
      if (infile) {
	while (getline(infile, strFileLine)) {
	  treefile2.push_back(strFileLine);
	  ulLineCount++;
	}
      }
      fprintf(stderr, "Number of trees in the second input file: %u\n", ulLineCount);
      NUM_TREES2 = ulLineCount;
      
      infile.close();
    }
    
  } catch (TCLAP::ArgException &e) { // catch any exceptions
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  if (!SET_OP){
    cerr << "No set operation specified! Exiting..." << endl;
    return 1;
  }
  fprintf(stderr, "outfilename is: %s\n", outfilename.c_str());
  //deal with null cases
  if (NUM_TREES == 0){
    ofstream fout;
    fout.open(outfilename.c_str());

    if (SET_OP == 1){ //union
      for (unsigned int i  = 0; i < NUM_TREES2; i++){
	fout << treefile2[i] << endl;
      }
    }
    fout.close();
    printf("= Total number of trees: %u\n ",NUM_TREES2);
    return 0;
  }
  else if (NUM_TREES2 == 0){
    ofstream fout;
    fout.open(outfilename.c_str());
    if ( (SET_OP == 1) || (SET_OP == 3) ){
       for (unsigned int i  = 0; i < NUM_TREES; i++){
	fout << treefile[i] << endl;
      }
    }
    fout.close();
    printf("= Total number of trees: %u\n ",NUM_TREES);
    return 0;
  }

  //process input file
  unsigned long long M1=0;
  unsigned long long M2=0;
  FILE *fp;
  NEWICKTREE *newickTree;
  int err;
  string tree; 
  struct timeval index_start;
  struct timeval index_end;
  struct timeval hash_start;
  struct timeval hash_end;
  struct timeval unique_start;
  struct timeval unique_end;


  gettimeofday(&hash_start, NULL);
  LabelMap lm = collect_labels(infilename);
  //NUM_TREES += NUM_TREES2;
  unsigned int TOTAL_TREES = NUM_TREES + NUM_TREES2;
  printf("= Total number of trees: %u\n ",TOTAL_TREES);
  initialize_hashtable(M1, M2); //initialize contents of hashtable
  //cout << "M1 is: " << M1 << endl;
  //cout << "M2 is: " << M2 << endl;

  //NUM_TREES -= NUM_TREES2;
  //process first file
  fp = fopen(infilename.c_str(), "r");
  if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}
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
    }
    else {
      unsigned int numBitstr=0;
      //bool processed = false;
      //bool collect = false;
      dfs_compute_hash(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, NULL);
      killnewicktree(newickTree);
    }
  }
  fclose(fp);


  //next, process the second file
  fp = fopen(infilename2.c_str(), "r");
  if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}
  for (unsigned int treeIdx=NUM_TREES; treeIdx < TOTAL_TREES; ++treeIdx) {
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
    }
    else {
      unsigned int numBitstr=0;
      //bool add = true;
      //bool processed = false;
      //bool collect = false;
      dfs_compute_hash(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, NULL);
      killnewicktree(newickTree);
    }
  }
  fclose(fp);

  gettimeofday(&hash_end, NULL);
  double hash_time = hash_end.tv_sec - hash_start.tv_sec + (hash_end.tv_usec - hash_start.tv_usec) / 1.e6;
  fprintf(stderr, "Hash table Creation Time: %g s\n", hash_time);

  unsigned int unique_bipartitions = 0;
  unsigned int max_bipart = NUM_TAXA - 3;
  gettimeofday(&index_start, NULL);
  unsigned int ** inverted_index = (unsigned int **)malloc(TOTAL_TREES*sizeof(unsigned int *));

  unsigned int * inverted_index_sizes = (unsigned int *)malloc(TOTAL_TREES*sizeof(unsigned int));

  for (unsigned int i = 0; i  < TOTAL_TREES; ++i){
    inverted_index[i] = (unsigned int *)malloc(max_bipart * sizeof(unsigned int));
    inverted_index_sizes[i] = 0;
  }


  unsigned int tree_id = 0;
  unsigned int loc = 0;
  unsigned int bipart_count = 0;
  for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); hti++) {
    unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size(); 
    if (sizeVec) {
      unique_bipartitions += sizeVec;
      for (unsigned int i = 0; i < sizeVec; i++){
	unsigned int sizeTreeIdx = vvec_hashrf._hashtab2[hti][i]._vec_treeidx.size();
	for (unsigned int j = 0; j < sizeTreeIdx; j++){
	  tree_id = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
	  loc = inverted_index_sizes[tree_id];
	  inverted_index[tree_id][loc] = bipart_count;
	  inverted_index_sizes[tree_id]++;
	}
	delete vvec_hashrf._hashtab2[hti][i]._bs;
	bipart_count++;
      }
    } //end if sizevec
  }

  /*for(unsigned int i = 0; i < TOTAL_TREES; i++){
    loc = inverted_index_sizes[i];
    for (unsigned int j = 0; j < loc; j++){
      cout << inverted_index[i][j] << " ";
    }
    cout << endl;
    }*/
  gettimeofday(&index_end, NULL);
  double index_time = index_end.tv_sec - index_start.tv_sec + (index_end.tv_usec - index_start.tv_usec) / 1.e6;
  fprintf(stderr, "Inverted Index Creation Time: %g s\n", index_time);

  fprintf(stderr, "%u unique bipartitions found\n", unique_bipartitions);
  assert(bipart_count == unique_bipartitions);
  //vvec_hashrf.hashrfmap_clear(); <-- need to fix this function!
  fprintf(stderr, "rehashing trees...\n");
   
  gettimeofday(&unique_start, NULL);
  HashRFMap vvec_setop;

  vvec_setop.uhashfunc_init_unique(TOTAL_TREES, NUM_TAXA, unique_bipartitions, C, NEWSEED);
  
  //new m1 and m2 values
  unsigned long long sm1 = vvec_setop._HF.getM1();
  unsigned long long sm2 = vvec_setop._HF.getM2();

  vvec_setop._hashtab2.resize(sm1); //resize hashtable
  // cout << "sm1 is: " << sm1 << endl;
  //cout << "sm2 is: " << sm2 << endl;

  for (unsigned int i = 0; i < TOTAL_TREES; i++){ 
    unsigned int pos = 0;
    loc = inverted_index_sizes[i];
    unsigned long long temp1 = 0;
    unsigned long long temp2 = 0;
    bool * cbs = new bool[unique_bipartitions];
    for (unsigned int j = 0; j < unique_bipartitions; j++) {
      cbs[j] = 0;
    }

    for (unsigned int j = 0; j < loc; j++){
      pos = inverted_index[i][j];
      cbs[pos] = 1;
      temp1+= vvec_setop._HF.getA1(pos);
      temp2+= vvec_setop._HF.getA2(pos);
    }

    temp1 = temp1 % sm1;
    temp2 = temp2 % sm2;
    vvec_setop.hashing_bs_unique(i, unique_bipartitions, temp1, temp2, cbs);
  }
  gettimeofday(&unique_end, NULL);
  double unique_time = unique_end.tv_sec - unique_start.tv_sec + (unique_end.tv_usec - unique_start.tv_usec) / 1.e6;
fprintf(stderr, "New hashtable creation time: %g s\n", unique_time);

  vector<unsigned int> my_ids;
  switch(SET_OP){
  case 1:
    my_ids = calculate_union(NUM_TREES2, vvec_setop);
    break;
  case 2:
    my_ids = calculate_intersection(NUM_TREES2, vvec_setop);
    break;
  case 3:
    my_ids = calculate_set_difference(NUM_TREES2, vvec_setop);
    break;
  default:
    cerr << "Error! No set operation specified!" << endl;
  }
 
  unsigned int id; 
  ofstream fout;
  fout.open(outfilename.c_str());
  if (!fout){
    cerr << "Cannot open file for writing!" << endl;
    return 2;
  }
  for (unsigned int i =0; i < my_ids.size(); i++){
    id = my_ids[i];
    if (id < NUM_TREES){
      fout << treefile[id] << endl;
    }
    else if (id >= NUM_TREES){
      id -= NUM_TREES;
      fout << treefile2[id] << endl;
    }
    else{
      cerr << "ERROR! Invalid Id!  " << id << endl;
      return 3;
    }
  }
  fout.close();
  return 0;
  
}
