#include <iostream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <sys/time.h>
#include <algorithm>
#include "sample.hh"
#include "global.h"
#include "label-map.hh"
#include "parsing.h"
#include "hash.hh"
#include "traverse.hh"
#include "compute.hh"

using namespace std;
bool mysort(float i, float j){return i < j;}
bool mysort2(float i, float j){return i > j;}

void setup_sampling(string infile, unsigned int dist, string outfile, unsigned int nsamples, bool sort_sampling){
  //first figure out the sample size from the number of samples
  if (NUM_TREES %nsamples != 0){
    cerr << "Error! Cannot equally divide input file containing " << NUM_TREES << " into " << nsamples << " samples! Please check your input and re-run!" << endl;
    exit(2);
  }
  unsigned int sample_size = NUM_TREES / nsamples;
  cout << "The sample size is: " << sample_size << endl;
  if (WEIGHTED) //determine if trees are weighted
    weighted_sampling(infile, dist, outfile, nsamples, sample_size, sort_sampling);
  else
    unweighted_sampling(infile, dist, outfile, nsamples, sample_size, sort_sampling);
}

void  add_and_compute_w(unsigned int begin, unsigned int offset, float ** sample_matrix, float * ab_array, unsigned int distance_measure, float ** float_matrix, unsigned int sample_size, vector< vector<float> > & stdev_vec){  
  double traversal_time = 0;
  struct timeval traversal_start;
  struct timeval traversal_end;
  struct timeval calculate_start;
  struct timeval calculate_end;
  unsigned int per_10, per_counter, per_place;
  per_10 = NUM_TREES/10;
  per_counter = 0;
  per_place = 0;
  gettimeofday(&traversal_start, NULL);
  gettimeofday(&calculate_start, NULL);  
  TOTAL_TREES = NUM_TREES; //now we are calculating one row at a time.
  bool is_distance_option = false;
  if ( (distance_measure ==1) || (distance_measure == 23) || (distance_measure == 24) || (distance_measure == 19)) 
    is_distance_option = true;
  if (is_distance_option)
    assert(float_matrix != NULL);
  else
    assert(float_matrix == NULL);
  //intialize RFMatrix row array
  float * RFMatrixRow = (float *)malloc(TOTAL_TREES*sizeof(float));
  if (RFMatrixRow == NULL){
    fprintf(stderr, "Cannot allocate RFMatrixRow: Out of Memory!\n");
    exit(3);
  }
  long num_comparisons = 0; 
  unsigned int hash_loc = 0;
  unsigned int temp = 0;

  if (TRIVIAL)
    NUMBIPART += NUM_TAXA;

  float x_weight, y_weight;
  unsigned int col_begin = offset;
  unsigned int row_end = begin+sample_size;
  unsigned int col_end = TOTAL_TREES;
  for (unsigned int pos =begin; pos < row_end; ++pos) { //for the rows between begin and end 
    for (unsigned int x = 0; x < TOTAL_TREES; ++x) { //initialize 
      RFMatrixRow[x] = 0.0;
    }  
    //calculate similarity (value: '\sum xiyi' for measures in Ellis et. al paper)
    //we only need to do this for similarity measures -- if it is a distance measure, forget it. \sum xiyi is not used.
    if ( !is_distance_option){
      for (unsigned int i = 0; i < helper_sizes[pos]; ++i) {    
	hash_loc = windex[pos][i].t_id; //get id
	x_weight = windex[pos][i].weight; //get weight
	for (unsigned int j = 0; j < hash_lengths[hash_loc]; ++j) { //hash line
	  temp = weighted_hashtable[hash_loc][j].t_id; //get tid
	  y_weight = weighted_hashtable[hash_loc][j].weight; //get weight
	  if ( (temp >= col_begin) && (temp < col_end)){
	    y_weight*=x_weight; // xiyi
	    //temp = temp % TOTAL_TREES;
	    ++num_comparisons;
	    RFMatrixRow[temp]+= y_weight; //sum xiyi
	  }
	}
      }
    }
    else if ( (distance_measure == 1) || (distance_measure == 23) || (distance_measure == 24) ){
      float my_temp = 0.0;
      for (unsigned int i = col_begin; i < TOTAL_TREES; i++){
	//true_val = col_begin+i;
	for (unsigned int j = 0; j < NUMBIPART; j++){
	  //cout << float_matrix[pos][j] << " and " << float_matrix[true_val][j] << endl;
	  my_temp = fabs(float_matrix[pos][j] - float_matrix[i][j]);
	  //cout << "my_temp is: " << my_temp << endl;
	  RFMatrixRow[i] += my_temp;
	}
      }
    }
    else if (distance_measure == 19){
      float my_temp = 0.0;
      for (unsigned int i = col_begin; i < TOTAL_TREES; i++){
	//true_val = end+i;
	for (unsigned int j = 0; j < NUMBIPART; j++){
	  my_temp = fabs(float_matrix[pos][j] - float_matrix[i][j]);
	  my_temp*=my_temp;
	  RFMatrixRow[i] += my_temp;
	}
      }
    }
    else{
      cerr << "INVALID OPTION! SHOULD NEVER GET HERE!" << endl;
      exit(1);
    }

    float value;
    stringstream rs;
    unsigned int startpos;
    unsigned int endpos;
    float main, ab, ac;
    startpos = 0;
    endpos = NUM_TREES;
    unsigned int true_x = begin / sample_size;
    for (unsigned int y = col_begin; y < TOTAL_TREES; ++y) { 
      //true_val = y+col_begin;
      main = RFMatrixRow[y];
      ab = ab_array[pos];
      ac = ab_array[y];
      value = main;
      //cout << "main is: " << main << endl;
      //cout << "ab is: " << ab << endl;
      //cout << "ac is: " << ac << endl;
      if ( (distance_measure != 1) && (distance_measure != 19) )
	calculatew(main,ab,ac,distance_measure, value);
      //cout << "value is: " << value << endl;
      unsigned int true_y = y / sample_size;
      if (value != 0){
	sample_matrix[true_x][true_y] += value;
	stdev_vec[true_y].push_back(value);
      }
      //unsigned int true_x = pos % TOTAL_TREES;
    }
    
    if (NUM_TREES > 100){
      if (pos % per_10 == 0){
	per_place = per_counter*10;
	fprintf(stderr, "%d%% Complete\n", per_place);
	per_counter++;
      }
    }
  }

  fprintf(stderr, "Finishing Traversal\n");  
  fprintf(stderr, "Number of comparisons: %lu\n", num_comparisons);
  gettimeofday(&traversal_end, NULL);
  gettimeofday(&calculate_end, NULL);  
  traversal_time = traversal_end.tv_sec - traversal_start.tv_sec + (traversal_end.tv_usec - traversal_start.tv_usec) / 1.e6; 
  double calculate_time = calculate_end.tv_sec - calculate_start.tv_sec + (calculate_end.tv_usec - calculate_start.tv_usec) / 1.e6; 
  fprintf(stderr, "\nTraversal Time: %g s\n", traversal_time); 
  fprintf(stderr, "\nRF Calculation Time: %g s\n", calculate_time); 
  free(RFMatrixRow);
}

void fill_stdev_matrix(float ** stdev_matrix, float ** sample_matrix, vector< vector<float> > & stdev_vec, unsigned int start_loc, unsigned int denom, unsigned int matLength){
  float temp = 0.0;
  float mean = 0.0;
  for (unsigned int i = start_loc; i < matLength; i++){
    mean = sample_matrix[start_loc][i]/denom; //calculate the true mean
    unsigned int mysize = stdev_vec[i].size();
    assert(mysize == denom);
    for (unsigned int j = 0; j < denom; j++){
      temp = stdev_vec[i][j] - mean;
      temp = temp * temp;
      stdev_matrix[start_loc][i] += temp;
    }
    temp = stdev_matrix[start_loc][i]/(denom-1);
    temp = sqrt(temp); //don't forget to calculate the square root!
    stdev_matrix[start_loc][i] = temp; // final stdev value
    sample_matrix[start_loc][i] = mean; //final average value
    stdev_vec[i].clear(); //empty that location in the vector
  }
}

void weighted_sampling(string infile, unsigned int distance_option, string outfile, unsigned int nsamples, unsigned int sample_size, bool sort_sampling){
  ROOTED = true;
  WEIGHTED = true;

  //declare index and hashtable
  float * ab_array = (float*)malloc(NUM_TREES*sizeof(float));
  for (unsigned int i  = 0; i < NUM_TREES; i++){
    ab_array[i] = 0.0;
  }
  assert(QOPT == false); //make sure there is no q-option turned on
  string qfile = "";
  process_trees_weighted(infile, qfile, distance_option, outfile, ab_array); //this populates the hashtable, row data structure, and the inverted index

  float ** float_matrix = NULL; //for certain distance measures
  bool is_distance_option = false;
  if ( (distance_option ==1) || (distance_option == 23) || (distance_option == 24) || (distance_option == 19)) {
    is_distance_option = true;
    //we need to calculate |xi - yi| here. So, we need to allocate a new data structure, and deallocate the hashtable
    float_matrix = (float ** )malloc(NUM_TREES*sizeof(float*));
    if (float_matrix == NULL) {
      cerr << "error! could not allocate memory!" << endl;
      exit(1);
    }
    for (unsigned int i =0; i <NUM_TREES; i++){
      float_matrix[i] = (float*)malloc(NUMBIPART*sizeof(float));
      if (float_matrix[i] == NULL){
	cerr << "error! could not allocate memory!" << endl;
	exit(1);
      }
      for (unsigned int j = 0; j < NUMBIPART; j++){
	float_matrix[i][j] = 0.0;
      }
    }
    //now, take the hashtable and populate the above matrix, and deallocate the hashtable
    //float_matrix == full inverted index (neede for weighted distance measures)
    unsigned int tree_loc;
    float my_weight;
    for (unsigned int i = 0; i < NUMBIPART; i++){
      for (unsigned int j = 0; j < hash_lengths[i]; j++){
	tree_loc = weighted_hashtable[i][j].t_id;
	my_weight = weighted_hashtable[i][j].weight;
	float_matrix[tree_loc][i] = my_weight;
	//cout << "float_matrix[" << tree_loc << "][" << i << "] " << my_weight << endl;
      }
    }

    //deallocate hashtable and index -- not need for weighted measures
    for (unsigned int i = 0; i < NUMBIPART; i++)
    free(weighted_hashtable[i]);
    free(weighted_hashtable);
    for (unsigned int i = 0; i < TOTAL_TREES; i++)
      free(windex[i]);
    free(windex);
  }
 
  unsigned int N_PARTS = 1000; //1000
  //average matrix
  float ** sample_matrix = (float **)malloc(N_PARTS*sizeof(float *));
  for (unsigned int i = 0; i < N_PARTS; i++){
    sample_matrix[i] = (float *)malloc(N_PARTS * sizeof(float));
    for (unsigned int j  =0; j < N_PARTS; j++)
      sample_matrix[i][j] = 0.0;
  }

  //stdev matrix
  float ** stdev_matrix = (float **)malloc(N_PARTS*sizeof(float *));
  for (unsigned int i = 0; i < N_PARTS; i++){
    stdev_matrix[i] = (float *)malloc(N_PARTS * sizeof(float));
    for (unsigned int j  =0; j < N_PARTS; j++)
      stdev_matrix[i][j] = 0.0;
  }

  unsigned int mysample_size = NUM_TREES/N_PARTS;
  unsigned int count = 0;
  unsigned int beg, end;
  beg = 0;
  end = mysample_size;
  vector< vector<float> > stdev_vec;
  unsigned int denom = mysample_size*mysample_size;
  stdev_vec.resize(N_PARTS); //set the length to N_PARTS
  for (unsigned int i = 0; i < NUM_TREES; i+=mysample_size){
    cout << "row " << count << endl;
    add_and_compute_w(beg, i, sample_matrix, ab_array, distance_option, float_matrix, mysample_size, stdev_vec);
    fill_stdev_matrix(stdev_matrix, sample_matrix, stdev_vec, count, denom, N_PARTS);
    beg += mysample_size;
    count++;
  }

  //float real_val = 0.0;
  ofstream fout;
  string stdev_out = outfile+".stdev";
  outfile += ".mean";
  fout.open(outfile.c_str());
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }

  for (unsigned int x = 0; x < N_PARTS; x++)
    fout << "c" << x << " ";
  fout << endl;
  for (unsigned int x = 0; x < N_PARTS; x++){
    fout << "c" << x << " ";
    for (unsigned int y  = 0; y < N_PARTS; y++){
      //real_val = sample_matrix[x][y] / nsamples;
      if (x <= y)
	fout << sample_matrix[x][y] << " ";
      else
	fout << sample_matrix[y][x] << " ";
    }
    fout << endl;
  }
  fout.close();

  fout.open(stdev_out.c_str());
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  for (unsigned int x = 0; x < N_PARTS; x++)
    fout << "c" << x << " ";
  fout << endl;
  for (unsigned int x = 0; x < N_PARTS; x++){
    fout << "c" << x << " ";
    for (unsigned int y  = 0; y < N_PARTS; y++){
      //real_val = sample_matrix[x][y] / nsamples;
      if (x <= y)
	fout << stdev_matrix[x][y] << " ";
      else
	fout << stdev_matrix[y][x] << " ";
    }
    fout << endl;
  }
  fout.close();

  for (unsigned int i  = 0; i < N_PARTS; i++){
    free(sample_matrix[i]);
    free(stdev_matrix[i]);
  }
  free(sample_matrix);
  free(stdev_matrix);
  if (is_distance_option){
    for (unsigned int i = 0; i < NUM_TREES; i++)
      free(float_matrix[i]);
    free(float_matrix);
  }
}

void add_and_compute_uw(unsigned int begin, unsigned int offset, float ** sample_matrix, unsigned int * row, unsigned int hash_length, unsigned int dist, unsigned int sample_size, vector< vector<float> > & stdev_vec){
  double traversal_time = 0;
  struct timeval traversal_start;
  struct timeval traversal_end;
  struct timeval calculate_start;
  struct timeval calculate_end;
  unsigned int per_10, per_counter, per_place;
  per_10 = NUM_TREES/10;
  per_counter = 0;
  per_place = 0;
  
  gettimeofday(&traversal_start, NULL);
  gettimeofday(&calculate_start, NULL);  
  long num_comparisons = 0;
  
  if (MULTIFURCATING)
    assert(n_biparts != NULL);

  unsigned int hash_loc = 0;
  unsigned int temp = 0;
  fprintf(stderr, "\nBeginning By-Tree Traversal:\n");
  ofstream fout;
  TOTAL_TREES = NUM_TREES;

  if (TRIVIAL)
    NUMBIPART += NUM_TAXA;
  unsigned int * RFMatrixRow = (unsigned int *)malloc(TOTAL_TREES*sizeof(unsigned int));
  if (RFMatrixRow == NULL){
    fprintf(stderr, "Cannot allocate RFMatrixRow: Out of Memory!\n");
    exit(3);
  }

  unsigned int row_end = begin+sample_size;
  unsigned int col_begin = offset;
  unsigned int col_end = TOTAL_TREES;

  for (unsigned int pos =begin; pos < row_end; ++pos) { //for the rows between begin and end 
    //fprintf(stderr,"\npos is: %d\n", pos);
    //calculate similarity (value: 'a' from Salim et. al paper)
    for (unsigned int x = 0; x < TOTAL_TREES; x++)
      RFMatrixRow[x]= 0;
    for (unsigned int i = 0; i < helper_sizes[pos]; ++i) {    //for all the bipartitions associated with that tree
      hash_loc = helper[pos][i]; //get id
      for (unsigned int j = 1; j < hash_lengths[hash_loc]; j++) { //hash line
	temp = hashtable[hash_loc][j]; //get the tree id
	if( (temp >= col_begin) && (temp < col_end) ){
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
    unsigned int true_x = pos / sample_size;
    if (MULTIFURCATING) {
      d = 0;
      for (unsigned int y = col_begin; y < TOTAL_TREES; y++) {
	if (y!=pos) { 	    
	  a = ALL_COUNT - row[pos] - row[y] + RFMatrixRow[y];
	  if (TRIVIAL){
	    a += NUM_TAXA; 
	    b = n_biparts[pos] + NUM_TAXA- a; //b 
	    c = n_biparts[y] + NUM_TAXA - a; //c
	  }
	  else{
	    b = n_biparts[pos] - a; //b 
	    c = n_biparts[y] - a; //c
	  }
	}
	else{
	  a = maxvalue;
	  if (TRIVIAL)
	    a += NUM_TAXA;
	  b = 0;
	  c = 0;
	}
	calculatem(a,b,c,d,dist,value2);
	unsigned int true_y = y / sample_size;
	if (value2 != 0){
	  sample_matrix[true_x][true_y] += value2;
	  stdev_vec[true_y].push_back(value2);
	}
      }
    }
    else { //binary
      d = 0;
      for (unsigned int y = col_begin; y < TOTAL_TREES; ++y) {
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
	calculate(a,b,d,dist, value, value2);
	unsigned int true_y = y / sample_size;
	//unsigned int true_x = pos % TOTAL_TREES;
	if ((dist == 1) || (dist == 24)){
	  if (value!= 0){
	    sample_matrix[true_x][true_y] += value;
	    stdev_vec[true_y].push_back(value);
	  }
	}
	else{
	  if (value2 != 0){
	    sample_matrix[true_x][true_y] += value2;
	    stdev_vec[true_y].push_back(value2);
	  }
	}
      } //for all values in the calculated row
    } //end if binary
    
    if (NUM_TREES > 100){
      if (pos % per_10 == 0){
        per_place = per_counter*10;
        fprintf(stderr, "%d%% Complete\n", per_place);
        per_counter++;
      }
    }
 }

 fprintf(stderr, "Finishing Traversal\n");  
 fprintf(stderr, "Number of comparisons: %lu\n", num_comparisons);
 gettimeofday(&traversal_end, NULL);
 gettimeofday(&calculate_end, NULL);  
  traversal_time = traversal_end.tv_sec - traversal_start.tv_sec + (traversal_end.tv_usec - traversal_start.tv_usec) / 1.e6; 
  double calculate_time = calculate_end.tv_sec - calculate_start.tv_sec + (calculate_end.tv_usec - calculate_start.tv_usec) / 1.e6; 
  fprintf(stderr, "\nTraversal Time: %g s\n", traversal_time); 
  fprintf(stderr, "\nRF Calculation Time: %g s\n", calculate_time); 
  free(RFMatrixRow);
  if (TRIVIAL)
    NUMBIPART -= NUM_TAXA;
  //free(bipart_counts);
}

void unweighted_sampling(string infile, unsigned int dist, string outfile, unsigned int nsamples, unsigned int sample_size, bool sort_sampling){
  unsigned int hashtable_length = 0;
  string qfile = "";
  unsigned int * row = (unsigned int *)malloc(NUM_TREES*sizeof(unsigned int));
  for (unsigned int i = 0; i < NUM_TREES; i++)
    row[i]=0;
  assert(QOPT == false);
  process_trees(infile, qfile, hashtable_length, row); //this populates the hashtable, row data structure, and the inverted index
  //now, allocate a sample_size x sample_size matrix
  unsigned int N_PARTS = 1000; //1000
  //average matrix
  float ** sample_matrix = (float **)malloc(N_PARTS*sizeof(float *));
  for (unsigned int i = 0; i < N_PARTS; i++){
    sample_matrix[i] = (float *)malloc(N_PARTS * sizeof(float));
    for (unsigned int j  =0; j < N_PARTS; j++)
      sample_matrix[i][j] = 0.0;
  }

  //stdev matrix
  float ** stdev_matrix = (float **)malloc(N_PARTS*sizeof(float *));
  for (unsigned int i = 0; i < N_PARTS; i++){
    stdev_matrix[i] = (float *)malloc(N_PARTS * sizeof(float));
    for (unsigned int j  =0; j < N_PARTS; j++)
      stdev_matrix[i][j] = 0.0;
  }

  unsigned int mysample_size = NUM_TREES/N_PARTS;
  unsigned int denom = mysample_size*mysample_size;
  unsigned int count = 0;
  unsigned int beg = 0;
  vector< vector<float> > stdev_vec;
  stdev_vec.resize(N_PARTS);
  cout << "mysample_size is: " << mysample_size << endl;
  for (unsigned int i = 0; i < NUM_TREES; i+=mysample_size){
    cout << "row " << count << endl;
    add_and_compute_uw(beg, i, sample_matrix, row, hashtable_length, dist, mysample_size, stdev_vec);
    fill_stdev_matrix(stdev_matrix, sample_matrix, stdev_vec, count, denom, N_PARTS);
    count++;
    beg += mysample_size;
  }

  ofstream fout;
  string stdev_out = outfile+".stdev";
  outfile += ".mean";
  fout.open(outfile.c_str());
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  for (unsigned int x  = 0; x < N_PARTS; x++)
    fout << "c" << x << " ";
  fout << endl;
  for (unsigned int x = 0; x < N_PARTS; x++){
    fout << "c" << x << " ";
    for (unsigned int y  = 0; y < N_PARTS; y++){
      if (x <= y)
	fout << sample_matrix[x][y] << " ";
      else
	fout << sample_matrix[y][x] << " ";
     }
    fout << endl;
  }
  fout.close();

  fout.open(stdev_out.c_str());
  if (!fout){
    cerr << "cannot open file for writing!" << endl;
    exit(1);
  }
  for (unsigned int x  = 0; x < N_PARTS; x++)
    fout << "c" << x << " ";
  fout << endl;
  for (unsigned int x = 0; x < N_PARTS; x++){
    fout << "c" << x << " ";
    for (unsigned int y  = 0; y < N_PARTS; y++){
      if (x <= y)
	fout << stdev_matrix[x][y] << " ";
      else
	fout << stdev_matrix[y][x] << " ";
     }
    fout << endl;
  }
  fout.close();

  for (unsigned int i  = 0; i < N_PARTS; i++){
    free(sample_matrix[i]);
    free(stdev_matrix[i]);
  }
  free(sample_matrix);
  free(stdev_matrix);
}
