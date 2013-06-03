/*****************************************************************/
/*This is ConvertDB

(c) 2010 ConvertDB: Suzanne Matthews
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

#include "./tclap/CmdLine.h"
#include "global.h"

using namespace std;
unsigned int * n_biparts;

void process_DB(string infilename, unsigned int * row, unsigned int & hashlength_correction, bool * compacted) {
  string fileLine;
  unsigned int size0, size1, threshold, temp;
  unsigned int * vec0, *vec1;
  unsigned int hashpos = 0;
  bool compact = false;
  hashlength_correction = 0;
  size0 = 0;
  size1=0;
  threshold = NUM_TREES/2;
  vec0 = (unsigned int *)malloc(threshold*sizeof(unsigned int));
  vec1 = (unsigned int *)malloc(threshold*sizeof(unsigned int));
  ifstream fin; 
  fin.open(infilename.c_str());
  if (!fin){cerr << "Error! Cannot open file for reading!\n"; exit(1);}
  
  for (unsigned int bidx = 0; bidx < NUMBIPART; bidx++) {
    //read a line from the file
    getline(fin, fileLine); 
    for (unsigned int i = 0; i < NUM_TREES; i++){
      if (fileLine[i] == '0'){
	if (size0 < threshold){
	  vec0[size0] = i;
	  size0++;
	}
      }
      else{
	if (size1 < threshold){
	  vec1[size1] = i;
	  size1++;
	}
	else{
	  compact = true;
	}
      }
    }
    if (compact == true){
      if (size0 != 0){	
	hashtable[hashpos]=(unsigned int*)malloc(size0*sizeof(unsigned int));
	hash_lengths[hashpos] = size0;
	compacted[hashpos] = true;
	for (unsigned int i = 0; i < size0; i++){
	  temp = vec0[i];
	  hashtable[hashpos][i] = temp;
	  row[temp]++;
	}
	hashpos++;
      }
      else
	hashlength_correction++;
      ALL_COUNT++;
    }
    else{
      if (size1 != 0){
	hashtable[hashpos] = (unsigned int *)malloc(size1*sizeof(unsigned int));
	hash_lengths[hashpos] = size1;
	for (unsigned int i = 0; i < size1; i++)
	  hashtable[hashpos][i] = vec1[i];
	hashpos++;
      }
      else
	hashlength_correction++;
    }
    size0 = 0;
    size1 = 0;
    compact = false;
    //if char is 0 and size0 is less than threshold, add id to vec0[size0]; size0++;
    //if char is 1 and size1 is less than threshold, add id to vec1[size1]; size1++;
  }

  fin.close();  
}

/**********************
 *
 * main function
 *
 *********************/


int main(int argc, char** argv)
{
  string outfilename = "createDB.out";
  string qfilename;
  string infilename;
  bool * compact;
  struct timeval prog_start;
  struct timeval prog_end;
  struct timeval process_start;
  struct timeval process_end;
  struct timeval print_start;
  struct timeval print_end;
  gettimeofday(&prog_start, NULL);
  // TCLAP
  //since we are creating the hashtable serially, the assumption in this
  //version is that the seed file is the same as the input file and that 
  //offset is set to 0
  //this can readily be changed in future version if necessary
  try {

    // Define the command line object.
    string helpMsg  = "createDB <input file>\n";
    helpMsg += "Input file: \n";
    helpMsg += "   createDB takes in a transposed bitstring matrix as input.\n";
    helpMsg += "If you have a selection of bitstrings, where the rows are the objects\n";
    helpMsg += "And the columns are the features, it is first necessary to transpose it. \n";
    helpMsg += " please use the transposeDB script for preprocessing.\n";
    helpMsg += "Examples: \n";
    helpMsg += " createDB  foo.bitstring //input is foo.bitstring\n";
    helpMsg += " createDB  foo.bitstring -o out //input is foo.bitstring, outfile is out\n";
 
    TCLAP::CmdLine cmd(helpMsg, ' ', "1.0.0");
    TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "file name", true, "intree", "Input tree file name"  );
    cmd.add( fnameArg );

    TCLAP::ValueArg<string> oArg("o", "output", "output file", false, "createDB.out", "output file name");
    cmd.add( oArg );
    cmd.parse( argc, argv );

    if (strcmp(oArg.getValue().c_str(), "createDB.out")!= 0)
      outfilename = oArg.getValue();

    infilename = fnameArg.getValue();

    if (NUMBIPART == 0) {
      string strFileLine;
      unsigned int ulLineCount;
      ulLineCount = 0;
      ifstream infile(infilename.c_str());
      if (infile) {
        while (getline(infile, strFileLine)) {
	  if (ulLineCount==0){ //do something to count the number of "trees/objects".
	    NUM_TREES = strlen(strFileLine.c_str()); //subtract one to get rid of terminal char
	  }
          ulLineCount++;
        }
      }
      NUMBIPART = ulLineCount;
      infile.close();
    }
  } catch (TCLAP::ArgException &e) { // catch any exceptions
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
  ofstream fout;
  fout.open(outfilename.c_str());
  if (!fout){cerr << "Error! Cannot open file for writing!\n"; exit(2);}

  gettimeofday(&process_start, NULL);
  fprintf(stderr, "Number of features in the input file: %u\n", NUMBIPART);
  fprintf(stderr, "Number of objects in the input file: %u\n", NUM_TREES);
  unsigned int * row; //stores the amount to be subtracted from rows at the end
  row = (unsigned int *)malloc(NUM_TREES * sizeof(unsigned int));
  hash_lengths = (unsigned int *)malloc(NUMBIPART* sizeof(unsigned int));
  compact = (bool *)malloc(NUMBIPART*sizeof(bool));
  for (unsigned int i = 0; i < NUMBIPART; i++)
    compact[i] = false;
  hashtable = (unsigned int **)malloc(NUMBIPART*sizeof(unsigned int*));
  unsigned int hashtable_correction = 0;
  for (unsigned int i = 0; i < NUM_TREES; i++) { 
    row[i] = 0;
  }   
  
  process_DB(infilename, row, hashtable_correction, compact);
  gettimeofday(&process_end, NULL);   
  double process_time = process_end.tv_sec - process_start.tv_sec + (process_end.tv_usec - process_start.tv_usec) / 1.e6;
  fprintf(stderr,"\nDB Processing Time: %g s\n", process_time);

  unsigned int hashtable_length = NUMBIPART - hashtable_correction;
  fprintf(stderr, "Number of Features corrected by %u elements\n", hashtable_correction);
  gettimeofday(&print_start, NULL);
  fout << "FEATURES:" << hashtable_length << endl;
  fout << "OBJECTS:" << NUM_TREES << endl;
  fout << "ALL_COUNT:" << ALL_COUNT << endl;
  fout << "ROW:"; 
  for (unsigned int i = 0; i < NUM_TREES; ++i) 
    fout << row[i] << " ";
  fout << endl;
  fout << "BEGIN COMPACT TABLE" << endl;
   for (unsigned int i = 0; i < hashtable_length; ++i) { 
     if (compact[i])
       fout << "-:" << hash_lengths[i] << ":";
     else 
       fout << "+:" << hash_lengths[i] << ":";
     for (unsigned int j = 0; j < hash_lengths[i]; ++j)  
       fout << hashtable[i][j] << " ";
     fout << endl;
   }
   fout << "END COMPACT TABLE" << endl;
   fout.close();  
   gettimeofday(&print_end, NULL);
  double print_time = print_end.tv_sec - print_start.tv_sec + (print_end.tv_usec - print_start.tv_usec) / 1.e6;
   fprintf(stderr,"\nPrint DB Time: %g s\n", print_time);
   gettimeofday(&prog_end, NULL);
   double prog_time = prog_end.tv_sec - prog_start.tv_sec + (prog_end.tv_usec - prog_start.tv_usec) / 1.e6;
   fprintf(stderr,"\nProgram Running Time: %g s\n", prog_time);
  return 0;

}

