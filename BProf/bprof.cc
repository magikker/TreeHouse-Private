//*****************************************************************/
/*
This is BProf, a program that profiles bipartitions in a tree.
It is part of the HaB|TAT suite, written by Suzanne J. Matthews

(c) 2010 BProf: Suzanne J. Matthews 

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

bool reverse_string (string i,string j) { return i.compare(j); }


/**********************
 *
 * main function
 *
 *********************/
int main(int argc, char** argv)
{
  string outfilename = "mybprof.out";
  string outfilename2 = "bitstringinfo.out";
  string infilename;
  bool REDUCED = 0;
  //struct timeval prog_start;
  //struct timeval prog_end;
  //struct timeval process_start;
  //struct timeval process_end;
  //gettimeofday(&prog_start, NULL);
  // TCLAP
  //since we are creating the hashtable serially, the assumption in this
  //version is that the seed file is the same as the input file and that 
  //offset is set to 0
  //this can readily be changed in future version if necessary
  try {

    // Define the command line object.
    string helpMsg  = "bprof <input file> <number of taxa> -s[seed value] -o [output file]\n";

    helpMsg += "Input file: \n";
    helpMsg += "   The current version of BProf only supports the Newick format.\n";

    helpMsg += "Example of Newick tree: \n";
    helpMsg += "   (('Chimp':0.052625,'Human':0.042375):0.007875,'Gorilla':0.060125,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
    helpMsg += "   ('Chimp':0.052625,('Human':0.042375,'Gorilla':0.060125):0.007875,\n";
    helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
    helpMsg += "Specify seed value: \n";   
    helpMsg += "   -s <value>, specify seed value (default: 1000) \n";
    helpMsg += "   -o <output file>, specifies the name of the output file. Default: mybprof.out \n";

    helpMsg += "Examples: \n";
    helpMsg += "  bprof foo.tre //input is foo.tre, output located in mybprof.out\n";
    helpMsg += "  bprof foo.tre -s 134//input is foo.tre, and the program is seeded with value 134. Output is located in mybprof.out\n";
    helpMsg += "  bprof foo.tre -o afile.txt -s 17 //input is foo.tre, seed is set to 17, and the output file should be afile.txt\n";
    
    TCLAP::CmdLine cmd(helpMsg, ' ', "1.0.0");

    TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "file name", true, "intree", "input file"  );
    cmd.add( fnameArg );

    //   TCLAP::ValueArg<unsigned int> cArg("c", "cvalue", "c value", false, 1000, "c value");
    //cmd.add( cArg );

    TCLAP::SwitchArg bArg("b", "hashtable", "print hashtable?", false);
    cmd.add( bArg );

    TCLAP::SwitchArg rArg("r", "reduced", "print reduced matrix?", false);
    cmd.add( rArg );

    TCLAP::ValueArg<string>  outputArg( "o", "output", "output file name", false, "mybprof.out", "user specified output file name"  );
    cmd.add( outputArg );

    TCLAP::ValueArg<int> seedArg("s", "seedvalue", "user specified seed value", false, 1000, "user specified seed value");
    cmd.add( seedArg );

    cmd.parse( argc, argv );

    NUM_TAXA = 0;
    NUM_TREES = 0;
    PRINT_HASH = 0;
    infilename = fnameArg.getValue();
    
    if (strcmp(outputArg.getValue().c_str(), "mybprof.out")){
      outfilename = outputArg.getValue();
    }

    if(bArg.getValue())
      PRINT_HASH = true;

    if (rArg.getValue())
      REDUCED = true;

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

    /*if (NUM_TREES < 2) {
      cerr << "Fatal error: at least two trees expected.\n"; 
      return 2;
      }*/

    //if (cArg.getValue())
    //  C = cArg.getValue();
    
    if (seedArg.getValue()) { 
      NEWSEED = seedArg.getValue();
    }

  } catch (TCLAP::ArgException &e) { // catch any exceptions
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  //process input file
  unsigned long long M1=0;
  unsigned long long M2=0;
  FILE *fp;
  NEWICKTREE *newickTree;
  int err;
  string tree; 
 
  vector<string> my_bitstrings;
  //char * bitstringstr = (char *)malloc((NUM_TAXA+1)*sizeof(char));
  //unsigned int x = 0;
  //char * line = (char*)malloc(200000*sizeof(char));
  //vector<pair<unsigned, unsigned> > mmap_cluster,
  //  vector<vector<SCNode*> > vvec_distinctClusters2,
  //collect labels, and compute hashtable (collecting bitstrings)
  LabelMap lm = collect_labels(infilename);
  initialize_hashtable(M1, M2); //initialize contents of hashtable
  fp = fopen(infilename.c_str(), "r");
  if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}
  //cout << endl;
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
      bool processed = false;
      dfs_compute_hash_unrooted2(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, processed, NULL);
      killnewicktree(newickTree);
    }
  }
  fclose(fp);

  unsigned int unique_bipartitions = 0;

  for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); ++hti) {
    unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size(); 
    if (sizeVec) {
      unique_bipartitions += sizeVec;
    } //end if sizevec
  }

  fprintf(stderr, "%u unique bipartitions found\n", unique_bipartitions);

  unsigned int tree_counter = 0;

  //allocate t x k matrix 
  bool ** tree_table = (bool **)malloc(NUM_TREES*sizeof(bool *));
  if (tree_table == NULL){
    cerr << "cannot allocate tree table!" << endl;
    return 2;
  }

  for (unsigned int i = 0; i < NUM_TREES; ++i) { 
    tree_table[i] = (bool*)malloc(unique_bipartitions*sizeof(bool)); 
    if (tree_table[i] == NULL) { 
      cerr << "cannot allocate row " << i << "of tree table!" << endl;
      return 2;
    }
  }
  
  //initialize everything in t x k table to 0
  for (unsigned int i = 0; i < NUM_TREES; ++i) {
    for (unsigned int j = 0; j < unique_bipartitions; ++j) { 
      tree_table[i][j] = 0;
    }
  }


  //first get labels
  ofstream fout;
   fout.open(outfilename2.c_str());
   if (!fout){
     cerr << "cannot open file bitstringinfo.out for writing!" << endl;
     return 3;
   }
   for (unsigned int i = 0; i < NUM_TAXA; ++i) { 
     fout << lm.name(i) << " ";
   }
   fout << endl;

  //traverse to fill out t x k matrix
  unsigned int uBID = 0;
  unsigned int bipart_count = 0;
  for (unsigned int hti=0; hti < vvec_hashrf._hashtab2.size(); ++hti) {
    unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size(); 
    if (sizeVec) {
      uBID += sizeVec;
      for (unsigned int i=0; i<sizeVec; ++i) {
	unsigned int sizeTreeIdx = vvec_hashrf._hashtab2[hti][i]._vec_treeidx.size();
	tree_counter += sizeTreeIdx;

	for (unsigned int j=0; j<NUM_TAXA; j++) {
	  fout << vvec_hashrf._hashtab2[hti][i]._bs[j];
	  if (PRINT_HASH)
	    cout << vvec_hashrf._hashtab2[hti][i]._bs[j];
	}
	fout << endl;
	if (PRINT_HASH)
	  cout << ":";
	for (unsigned int j = 0; j < sizeTreeIdx; ++j) {	
	  unsigned int temp  = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
	  if (PRINT_HASH)
	  cout << temp << " ";
	  tree_table[temp][bipart_count] = true;
	}
	if (PRINT_HASH)
	  cout << endl;
	bipart_count++;
      }//for everything in sizeVec
    } // if (sizeVec)
  } //end hashtable traversal
  fout.close();
  assert(unique_bipartitions == uBID);

  //gettimeofday(&process_end, NULL);

  //double process_time = process_end.tv_sec - process_start.tv_sec + (process_end.tv_usec - process_start.tv_usec) / 1.e6;
  //fprintf(stderr,"\nInput Processing Time: %g s\n", process_time);

  if (tree_counter < ((NUM_TAXA-3)*NUM_TREES)) { 
    fprintf(stderr, "Trees are Multifurcating\n");
  }
  else { 
    fprintf(stderr, "Trees are Binary\n");
  }

  fprintf(stderr, "Writing Results to file...");


  fout.open(outfilename.c_str());
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    return 3;
  }
  //traverse table and output to file
  /*fout << "\t";
  for (unsigned int i = 0; i < unique_bipartitions; ++i) { 
    fout << "B" << i << "\t";
  }
  fout << endl;
  */
  
 
  for (unsigned int i = 0; i < NUM_TREES; ++i) {
    //fout << "T" << i << ":\t";
    for (unsigned int j = 0; j < unique_bipartitions; ++j) { 
      fout << tree_table[i][j] << " ";
    }
    fout << endl;
  }
 
  
  fout.close();
  fprintf(stderr, "Done. Results located in file %s\n", outfilename.c_str());
 
  return 0;
  
}
