//*****************************************************************/
/*
This is HashCS, a program that calculates consensus trees
The original HashCS code is written by SeungJin Sul
This implementation was created to be part of the HaB|TAT suite,
written by Suzanne J. Matthews.

(c) 2010 HaB|TAT: Suzanne J. Matthews
(c) 2009 HashCS: SeungJin Sul
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
#include "hashmap.hh"
#include "./tclap/CmdLine.h"
#include "global.h"
#include "parsing.h"
#include "buildtree.h"

// For newick parser
extern "C" {
#include <newick.h>
}


using namespace std;

bool sort_bitstrings(const pair<unsigned, unsigned > & a, const pair<unsigned, unsigned> &b) { 
  return a.first > b.first;
}

/**********************
 *
 * DFS for Majority
 *
 *********************/
bool *
dfs_hashcs_MJ2(
	       NEWICKNODE *startNode,
	       LabelMap &lm,
	       HashMap &vvec_hashcs,
	       unsigned int treeIdx,
	       unsigned long long m1,
	       unsigned long long m2,
	       vector<bool *> &vec_bs,
	       unsigned &numBPCollected, 
	       unsigned int RES_NUMTREES
)
{
    if (treeIdx < RES_NUMTREES) {  // treeIdx < NUM_TREES/2, if r=0.5 ==> do not collect bitstrings
      if (startNode->Nchildren == 0) { // leaf node
	string temp(startNode->label);
	unsigned int idx = lm[temp];
	
	startNode->hv1 = vvec_hashcs._HF.getA1(idx);
	startNode->hv2 = vvec_hashcs._HF.getA2(idx);
	
	return NULL; // Just return NULL
      }
      else { // non-leaf node
	for (int i=0; i < startNode->Nchildren; i++) {
	  dfs_hashcs_MJ2(startNode->child[i], lm, vvec_hashcs, treeIdx, m1, m2, vec_bs, numBPCollected, RES_NUMTREES);
	}
	
	++numBPCollected;

	unsigned long long temp1=0;
	unsigned long long temp2=0;
	
	for (int i=0; i<startNode->Nchildren; ++i) {
	  temp1 += startNode->child[i]->hv1;
	  temp2 += startNode->child[i]->hv2;
	}
	startNode->hv1 = temp1 % m1;
	startNode->hv2 = temp2 % m2;

	// Just populate the hash table with hv1 and hv2
	if (numBPCollected < NUM_TAXA-2) {
	  vvec_hashcs.hashing_bs_MJ(startNode->hv1, startNode->hv2, RES_NUMTREES);
	}
 
	return NULL; // Just return NULL
      }
    }
    else { // treeIdx >= NUM_TREES/2 ==> collect bipartitions
      if (startNode->Nchildren == 0) { // leaf node

	// Get the location of taxon label in lm
	// which is used to set "1" in bitstring.
	string temp(startNode->label);
	unsigned idx = lm[temp];
	
	startNode->hv1 = vvec_hashcs._HF.getA1(idx);
	startNode->hv2 = vvec_hashcs._HF.getA2(idx);

	// Here we make an actual bitstring
	bool *bs = new bool[NUM_TAXA];           
	for (unsigned int i = 0; i < NUM_TAXA; i++)
	  bs[i] = 0;
	bs[idx] = 1;
	
	int numOnes = 0;
	for (unsigned int i =0; i < NUM_TAXA; i++)
	  numOnes += bs[i];
	if (numOnes > 1 || numOnes == 0) {
	  cout << "numOnes = " << numOnes << endl;
	  cerr << "ERROR: bitstring creation error\n";
	  exit(0);
	}
           		       
	return bs; // This is actual bistring
      }
      else  {
	// non-leaf node	
	bool *ebs[startNode->Nchildren];
	for (int i=0; i<startNode->Nchildren; i++) {
	  ebs[i] = dfs_hashcs_MJ2(startNode->child[i], lm, vvec_hashcs, treeIdx, m1, m2, vec_bs, numBPCollected, RES_NUMTREES);
	}

	// At this point, we find a bipartition.
	// Thus, OR the bitstrings and make a bit string for the bipartition
	bool *bs = new bool[NUM_TAXA];
	for (unsigned int i =0; i < NUM_TAXA; i++)
	  bs[i] = 0;
	for (int i=0; i<startNode->Nchildren; i++) {
	  if (ebs[i]) {
	    for (unsigned int j = 0; j < NUM_TAXA; j++)
	      bs[j] |= ebs[i][j]; // This bit operation creates bitstring for a internal node.
	    delete [] ebs[i];
	    ebs[i] = NULL;
	  }
	  else {
	    cout << "ERROR: null bitstring\n";
	    exit(0);
	  }
	}
						
	int numOnes = 0;
	for (unsigned int i = 0; i < NUM_TAXA; i++)
	  numOnes += bs[i];
	if (numOnes < 1) {
	  cout << "numOnes = " << numOnes << endl;
	  cerr << "ERROR: bitstring OR error\n";
	  exit(0);
	}	
	++numBPCollected;
	
	unsigned long long temp1=0;
	unsigned long long temp2=0;
	
	for (int i=0; i<startNode->Nchildren; i++) {
	  temp1 += startNode->child[i]->hv1;
	  temp2 += startNode->child[i]->hv2;
	}
	startNode->hv1 = temp1 % m1;
	startNode->hv2 = temp2 % m2;

	// Here we need to store the actial bitsting collected.
	// *bs is the bitstring. In hashing_bs_MJ, the bitstring is checked if it is
	// majority bipartition or not. If yes, the bipartition is stored in vec_bs.
	// vec_bs is later used for constructing majority consensus tree.
	if (numBPCollected < NUM_TAXA-2) {
	  vvec_hashcs.hashing_bs_MJ(bs, NUM_TAXA, startNode->hv1, startNode->hv2, vec_bs, RES_NUMTREES);
	}			
	return bs;
      }
    }
}
/**********************
 *
 * DFS for Strict
 *
 *********************/
bool *
dfs_hashcs_SC2(
    NEWICKNODE* startNode,
    LabelMap &lm,
    HashMap &vvec_hashcs,
    unsigned int treeIdx,
    unsigned long long m1,
    unsigned long long m2,
    vector<bool *> & vec_bs,
     unsigned &numBPCollected)
{
    if (treeIdx == 0) { // for the first tree ==> collect bipartitions
        // If the node is leaf node, just set the place of the taxon name in the bit string to '1'
        // and push the bit string into stack
        if (startNode->Nchildren == 0) { // leaf node
            // Get the location of taxon label in lm
            // which is used to set "1" in bitstring.
            string temp(startNode->label);
	    //cout << "At Leaf node: " << temp << endl;
            unsigned int idx = lm[temp];

            // Implicit BPs /////////////////////////
            // Set the hash values for each leaf node.
	    unsigned long long temp1, temp2;
            temp1 = vvec_hashcs._HF.getA1(idx);
	    temp2 = vvec_hashcs._HF.getA2(idx);
	    temp1 = temp1 % m1;
	    temp2 = temp2 % m2;
	    startNode->hv1 = temp1;
	    startNode->hv2 = temp2;
            // Here we make an actual bitstring
	    bool * bs = new bool[NUM_TAXA];
	    for (unsigned int i =0; i < NUM_TAXA; i++)
	      bs[i] = 0;
            bs[idx] = 1;
            
	    //cout << "Bitstring is: " << endl;
	    //for (unsigned int i = 0; i < NUM_TAXA;i++)
	    //cout << bs[i];
	    //cout << endl;

            int numOnes = 0;
            for (unsigned int i =0; i < NUM_TAXA; i++)
	      numOnes += bs[i];

            if (numOnes > 1 || numOnes == 0) {
	      cout << "numOnes = " << numOnes << endl;
	      cerr << "ERROR: bitstring creation error\n";
	      exit(0);
            }
	    
            return bs; // return the bitstring for this leaf node (=taxon)
        }
        else { // non-leaf node
	  //cout << "At internal node!" << endl;
	  bool * ebs[startNode->Nchildren];
	  for (int i=0; i< startNode->Nchildren; i++) {
	    ebs[i] = dfs_hashcs_SC2(startNode->child[i], lm, vvec_hashcs, treeIdx, m1, m2, vec_bs, numBPCollected);
	  }

	  // At this point, we find a bipartition. Thus, OR the bitstrings and make a bit string for the bipartition
	  bool * bs = new bool[NUM_TAXA];
	  for (unsigned int i =0; i < NUM_TAXA; i++)
	    bs[i] = 0;

	  for (int i=0; i<startNode->Nchildren; i++) {
	      if (ebs[i]) {
		for (unsigned int j = 0; j < NUM_TAXA; j++) 
		  bs[j] |= ebs[i][j];
		delete [] ebs[i];
		ebs[i] = NULL;
	      }
	      else {
		cout << "ERROR: null bitstring\n";
		exit(0);
	      }
	  }
						
	  int numOnes = 0;
	  for (unsigned int i = 0; i < NUM_TAXA; i++)
	    numOnes += bs[i];
	  if (numOnes < 1) {
	    cout << "numOnes = " << numOnes << endl;
	    cerr << "ERROR: bitstring OR error\n";
	    exit(0);
	  }
	  
	  ++numBPCollected;
	  
	  // Implicit BPs ////////////
	  // After an internal node is found, compute the hv1 and hv2
	  unsigned long long temp1=0;
	  unsigned long long temp2=0;
	  
	  for (int i=0; i<startNode->Nchildren; ++i) {
	    temp1 += startNode->child[i]->hv1;
	    temp2 += startNode->child[i]->hv2;
	  }
	  startNode->hv1 = temp1 % m1;
	  startNode->hv2 = temp2 % m2;
	  
	  // Store bit strings in hash table
	  // Here we need to store the actial bitsting collected.
	  // *bs is the bitstring. In hashing_bs_SC, the bitstring is checked if it is
	  // strict bipartition or not. If yes, the bipartition is stored in vec_bs.
	  // vec_bs is later used for constructing strict consensus tree.
	  // vvec_hashcs.hashing_bs_SC(*bs, treeIdx, NUM_TAXA, NUM_TREES, startNode->hv1, startNode->hv2, vec_bs);
	  if (numBPCollected < NUM_TAXA-2) {
	    vvec_hashcs.hashing_bs_SC(bs, NUM_TAXA, startNode->hv1, startNode->hv2);
	  } 
	  
	  return bs;
        }
    }
    else { // for the other trees
      if (startNode->Nchildren == 0) { // leaf node
	string temp(startNode->label);
	unsigned int idx = lm[temp];
	
	// Implicit BPs /////////////////////////
	// Set the hash values for each leaf node.
	startNode->hv1 = vvec_hashcs._HF.getA1(idx);
	startNode->hv2 = vvec_hashcs._HF.getA2(idx);	
	return NULL;
      }
      else { // non-leaf node
	for (int i=0; i<startNode->Nchildren; ++i) {
	  dfs_hashcs_SC2(startNode->child[i], lm, vvec_hashcs, treeIdx, m1, m2, vec_bs, numBPCollected);
	}
	
	++numBPCollected;
						
	// Implicit BPs ////////////
	// After an internal node is found, compute the hv1 and hv2
	unsigned long long temp1=0;
	unsigned long long temp2=0;
	
	for (int i=0; i<startNode->Nchildren; ++i) {
	  temp1 += startNode->child[i]->hv1;
	  temp2 += startNode->child[i]->hv2;
	}
	startNode->hv1 = temp1 % m1;
	startNode->hv2 = temp2 % m2;
	
	// Store bit strings in hash table
	//vvec_hashcs.hashing_bs_SC(*bs, treeIdx, NUM_TAXA, NUM_TREES, startNode->hv1, startNode->hv2, vec_bs);
	if (numBPCollected < NUM_TAXA-2) {
	  vvec_hashcs.hashing_bs_SC(treeIdx, NUM_TAXA, NUM_TREES, startNode->hv1, startNode->hv2, vec_bs);
	}
	
	return NULL;
      }
    }
}

/**********************
 *
 * main function
 *
 *********************/
int main(int argc, char** argv)
{
  string infilename;
  string outfilename = "consensus.out";
  bool f_mj = false;
  float  RES_RATE  = 0.5;  // Default resolution rate  
  unsigned int  RES_NUMTREES    = 0;  // boundary number of trees

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
    helpMsg += "  bprof foo.tre foo2.tre -u -o afile.txt -s 17 //input is foo.tre and foo2.tre, seed is set to 17. Find the union and place the output in afile.txt\n";
    
    TCLAP::CmdLine cmd(helpMsg, ' ', "1.5.0");

    TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "file name", true, "intree", "input file"  );
    cmd.add( fnameArg );

    TCLAP::SwitchArg jArg("j", "majority", "compute majority consensus tree", false);
    cmd.add( jArg );

    TCLAP::ValueArg<int> rArg("r", "rvalue", "rate value", false, 50, "rate value");
    cmd.add( rArg );

    TCLAP::ValueArg<string>  outputArg( "o", "output", "output file name", false, "consensus.out", "user specified output file name"  );
    cmd.add( outputArg );

    TCLAP::ValueArg<int> seedArg("s", "seedvalue", "user specified seed value", false, 1000, "user specified seed value");
    cmd.add( seedArg );

    cmd.parse( argc, argv );

    NUM_TAXA = 0;
    NUM_TREES = 0;

    infilename = fnameArg.getValue();

    if (strcmp(outputArg.getValue().c_str(), "consensus.out")){
      outfilename = outputArg.getValue();
    }

    if (seedArg.getValue()) { 
      NEWSEED = seedArg.getValue();
    }

    if (rArg.getValue()) {
      if (rArg.getValue() < 0 ||  rArg.getValue() > 100) {
	cout << "ERROR: 0 < rate <= 100\n";
	exit(0);
      }
    }

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

    if( jArg.getValue())
      f_mj = true;
    

    if (f_mj) {
      RES_RATE = float(rArg.getValue())/ 100.0;
      RES_NUMTREES = NUM_TREES * RES_RATE;
      if (RES_NUMTREES == NUM_TREES)
	f_mj = false;
    }
  } catch (TCLAP::ArgException &e) { // catch any exceptions
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  //process input file
  HashMap vvec_hashcs;
  unsigned long long M1=0;
  unsigned long long M2=0;
  FILE *fp;
  NEWICKTREE *newickTree;
  int err;
  string tree; 
 
  LabelMap lm = collect_labels(infilename);
  if (NEWSEED != 1000)
    vvec_hashcs.UHashfunc_init(NUM_TREES, NUM_TAXA, C, NEWSEED);
  else
    vvec_hashcs.UHashfunc_init(NUM_TREES, NUM_TAXA, C); 
  M1 = vvec_hashcs._HF.getM1();
  M2 = vvec_hashcs._HF.getM2();
    
  // Increase the size of hash table to M1
  vvec_hashcs._hashtab_hashCS.resize(M1); 

  // to store collected strict bipartitions
  // NOTE: strict/majority bipartitions are not stored in hash table.
  vector<bool *> vec_bs; 
  assert(vec_bs.size() == 0);    


  //process first file
  fp = fopen(infilename.c_str(), "r");
  if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}

  if (f_mj) { // MAJORITY
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
	unsigned int numBPCollected=0;
	dfs_hashcs_MJ2(newickTree->root, lm, vvec_hashcs, treeIdx, M1, M2, vec_bs, numBPCollected, RES_NUMTREES);                
	killnewicktree(newickTree);
      }
    }
  }
  else { // STRICT
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
	unsigned int numBPCollected=0;
	dfs_hashcs_SC2(newickTree->root, lm, vvec_hashcs, treeIdx, M1, M2, vec_bs, numBPCollected);
	killnewicktree(newickTree);
      }
    }
  }

  vector<bool *> list_bs;
  bool *bs2 = new bool[NUM_TAXA];
  for (unsigned int i =0; i < NUM_TAXA; i++)
    bs2[i] = 1;
  list_bs.push_back(bs2);
  vector<float> list_branches; //empty for now
  //sort by number of ones
  cout << "vec_bs.size() is: " << vec_bs.size() << endl;
  vector< pair<unsigned,  unsigned > > ordered_bitstrings;
  for (unsigned i = 0; i < vec_bs.size(); i++){
    unsigned num_ones = 0;
    bool * temp = vec_bs[i];
    for (unsigned int j=0; j<NUM_TAXA; j++) {
      num_ones += temp[j];
    }
    ordered_bitstrings.push_back(make_pair(num_ones, i));
  }
  sort(ordered_bitstrings.begin(), ordered_bitstrings.end(), sort_bitstrings);
  vector<pair<unsigned, unsigned> >::iterator o_itr = ordered_bitstrings.begin();
  unsigned loc;
  while (o_itr != ordered_bitstrings.end()){
    loc = o_itr->second;
    bool *temp = vec_bs[loc];
    list_bs.push_back(temp);
    o_itr++;
  }
  string consensus = compute_tree(lm, list_bs, list_branches, 0, false);

  //  unsigned int unique_bipartitions = 0;
  // unsigned int max_bipart = NUM_TAXA - 3;

  //fprintf(stderr, "%u unique bipartitions found\n", unique_bipartitions);
  unsigned int my_rate = RES_RATE*100;
  if (!f_mj)
    my_rate = 100;
  fprintf(stderr, "Resolution rate: %u%\n", my_rate);
  ofstream fout;
  fout.open(outfilename.c_str());
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    return 2;
  }
  fout << consensus << endl;
  fout.close();
 
  return 0;
  
}
