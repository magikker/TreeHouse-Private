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
#include "sample.hh"
#include "traverse.hh"

// For newick parser
extern "C" {
#include <newick.h>
}

using namespace std;

bool is_weighted(string filename){ //central assumption is that trees are outputted from phylogenetic search. So if first tree is weighted, all of them are.
  ifstream fin;
  string line;
  fin.open(filename.c_str());
  getline(fin, line);
  int pos = line.find_first_of(":");
  if (pos != -1)
    return true;
  return false;
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
  unsigned int nsamples = 0;
  float * ab_array;
  struct timeval prog_start;
  struct timeval prog_end;
  struct timeval process_start;
  struct timeval process_end;
  bool sort_sample = false;
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

    TCLAP::CmdLine cmd(helpMsg, ' ', "1.0.0");

    TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "file name", true, "intree", "Input tree file name"  );
    cmd.add( fnameArg );

    TCLAP::ValueArg<string>  qnameArg( "q", "qfile", "q file", false, "", "optional second tree file name"  );
    cmd.add( qnameArg );

    TCLAP::UnlabeledValueArg<int>  numtreeArg( "numtree", "number of trees", true, 2, "Number of trees"  );
    cmd.add( numtreeArg );

    TCLAP::ValueArg<unsigned int> cArg("c", "cvalue", "c value", false, 1000, "c value");
    cmd.add( cArg );

    TCLAP::ValueArg<string> oArg("o", "output", "output file", false, "phlash.out", "output file name");
    cmd.add( oArg );

    TCLAP::SwitchArg mArg("m", "matrix", "print matrix?", false);
    cmd.add( mArg );

    TCLAP::SwitchArg fArg("f", "fullmatrix", "calculate full matrix?", false);
    cmd.add( fArg );

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

    TCLAP::ValueArg<int> xArg("x", "sample", "number of samples", false, 0, "a number of samples for sampling procedure");
    cmd.add( xArg );

    TCLAP::SwitchArg yArg("y", "sort", "sorted sampling?", false);
    cmd.add( yArg );
    cmd.parse( argc, argv );

    NUM_TREES = numtreeArg.getValue();
    int toprint = mArg.getValue();
    int toprint_list = lArg.getValue();
    nsamples = xArg.getValue();
    if (yArg.getValue())
      sort_sample = true;
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
    cout << "nsamples is: " << nsamples << endl; 
    infilename = fnameArg.getValue();

    qfilename = qnameArg.getValue();
    if (qfilename != ""){
      QOPT = true;
      if (nsamples != 0){
	cerr << "cannot use q option with sampling option!" << endl;
	return 2;
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

  } catch (TCLAP::ArgException &e) { // catch any exceptions
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  if (nsamples != 0){
    //follow sampling procedure
    if (is_weighted(infilename))
      WEIGHTED = true;
    if (sort_sample)
      cout << "we should sort the samples!." << endl;
    setup_sampling(infilename, distance_option, outfilename, nsamples, sort_sample);
    return 0;
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

  if (is_weighted(infilename)){
    fprintf(stderr,"Trees are detected as being weighted!\n");
    ROOTED = true; //mark trees as rooted, weighted trees
    WEIGHTED = true;
    ab_array = (float *)malloc(rowsize*sizeof(float));
    for (unsigned int i = 0; i < rowsize; i++)
      ab_array[i] = 0.0;
    process_trees_weighted(infilename, qfilename, distance_option, outfilename, ab_array);
    weighted_row_distance(ab_array, distance_option, outfilename);
    return 0;
  }
  row = (unsigned int *)malloc(rowsize * sizeof(unsigned int));
  //unsigned int threshold = NUM_TREES - (NUM_TREES*P)/100;
  //threshold--;
  for (unsigned int i = 0; i < rowsize; i++) { 
    row[i] = 0;
  }   
  if (MC)
    cout << "NOTE: Monte-Carlo Option on!" << endl;
  process_trees(infilename, qfilename, hashtable_length, row);
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

