#include <iostream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include "parsing.h"
#include "hash.hh"
#include "traverse.hh"
#include "compute.hh"
#include <sys/time.h>

using namespace std;
bool ROOTED = false;
bool FULLMATRIX = false;
bool QOPT = false;
bool TRIVIAL = false;
unsigned int NUM_TREES2;
unsigned int TOTAL_TREES;
unsigned int * n_biparts;
float * ab_array;
struct weighted_element ** windex;
struct weighted_element ** weighted_hashtable;

void weighted_row_distance(float * ab_array, unsigned int distance_measure, string outfilename) { 
  if (PRINT || PRINT_L)
    fprintf(stderr, "outfilename is: %s\n", outfilename.c_str()); 
  
  double traversal_time = 0;
  struct timeval traversal_start;
  struct timeval traversal_end;
  struct timeval calculate_start;
  struct timeval calculate_end;
  struct timeval print_start;
  struct timeval print_end;
  double print_temp;
  unsigned int per_10, per_counter, per_place;
  double sum = 0;
  unsigned int count = 0;
  double average = 0;

  per_10 = NUM_TREES/10;
  per_counter = 0;
  per_place = 0;

  gettimeofday(&traversal_start, NULL);
  gettimeofday(&calculate_start, NULL);  
  //intialize RFMatrix row array
  float * RFMatrixRow = (float *)malloc(NUM_TREES*sizeof(float));
  if (RFMatrixRow == NULL){
    fprintf(stderr, "Cannot allocate RFMatrixRow: Out of Memory!\n");
    exit(3);
  }
  long num_comparisons = 0; 
  unsigned int hash_loc = 0;
  unsigned int temp = 0;
  ofstream fout;
  
  if (PRINT || PRINT_L){
    fout.open(outfilename.c_str());
    if (!fout){
      cerr << "ERROR! Cannot open file for writing!\n";
      exit(3);
    }
  }
 
  //reset NUM_TREES, and set TOTAL_TREES if necessary
  if (!QOPT)
    TOTAL_TREES = NUM_TREES;
  else
    NUM_TREES = TOTAL_TREES - NUM_TREES2;

 if (TRIVIAL)
    NUMBIPART += NUM_TAXA;

 float x_weight, y_weight;
 float ** float_matrix = NULL; //for certain distance measures
 bool is_distance_option = false;
 if ( (distance_measure ==1) || (distance_measure == 23) || (distance_measure == 24) || (distance_measure == 19)) 
   is_distance_option = true;
 if ( is_distance_option){
   //we need to calculate |xi - yi| here. So, we need to allocate a new data structure, and deallocate the hashtable
   float_matrix = (float ** )malloc(TOTAL_TREES*sizeof(float*));
   if (float_matrix == NULL) {
     cerr << "error! could not allocate memory!" << endl;
     exit(1);
   }
   for (unsigned int i =0; i <TOTAL_TREES; i++){
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
   unsigned int tree_loc;
   float my_weight;
   for (unsigned int i = 0; i < NUMBIPART; i++){
     for (unsigned int j = 0; j < hash_lengths[i]; j++){
       tree_loc = weighted_hashtable[i][j].t_id;
       my_weight = weighted_hashtable[i][j].weight;
       float_matrix[tree_loc][i] = my_weight;
     }
   }

   //deallocate hashtable and index -- not need for weighted measures
   for (unsigned int i = 0; i < NUMBIPART; i++)
     free(weighted_hashtable[i]);
   for (unsigned int i = 0; i < TOTAL_TREES; i++)
     free(windex[i]);
 }

 for (unsigned int pos =0; pos < NUM_TREES; ++pos) { //for every row 
   for (unsigned int x = 0; x < TOTAL_TREES; ++x) { //initialize 
     RFMatrixRow[x] = 0.0;
   }
   
    //calculate similarity (value: '\sum xiyi' for measures in Ellis et. al paper)
    //we only need to do this for similarity measures -- if it is a distance measure, forget it. \sum xiyi is not used.
    if ( !is_distance_option){
      for (unsigned int i = 0; i < helper_sizes[pos]; ++i) {    
	hash_loc = windex[pos][i].t_id; //get id
	x_weight = windex[pos][i].weight;
	for (unsigned int j = 0; j < hash_lengths[hash_loc]; ++j) { //hash line
	  temp = weighted_hashtable[hash_loc][j].t_id;
	  y_weight = weighted_hashtable[hash_loc][j].weight;
	  y_weight*=x_weight; // xiyi
	  RFMatrixRow[temp]+= y_weight; //sum xiyi
	  ++num_comparisons;
	}
      }
    }
    else if ( (distance_measure == 1) || (distance_measure == 23) || (distance_measure == 24) ){
      for (unsigned int i = 0; i < TOTAL_TREES; i++){
	for (unsigned int j = 0; j < NUMBIPART; j++){
	  RFMatrixRow[i]+= fabs(float_matrix[pos][j] - float_matrix[i][j]);
	}
      }
    }
    else if (distance_measure == 19){
      float my_temp = 0.0;
      for (unsigned int i = 0; i < TOTAL_TREES; i++){
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
    if (QOPT){
      startpos = NUM_TREES;
      endpos = TOTAL_TREES;
    }
    else{
      startpos = 0;
      endpos = NUM_TREES;
    }

    if (FULLMATRIX){ 
      for (unsigned int y = startpos; y < endpos; ++y) { 
	main = RFMatrixRow[y];
	ab = ab_array[pos];
	ac = ab_array[y];
	value = main;
	if ( (distance_measure != 1) && (distance_measure != 19) )
	  calculatew(main,ab,ac,distance_measure, value);
	if (!PRINT_L){
	  if (PRINT){
	    rs << value << " ";
	  }
	  else{  
	    if (value != NAN){
	      sum += value;
	      count++;
	    }
	  }
	}
	else{ //PRINT_L
	  rs << pos << "." << y << " " << value << "\n";
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
    else{ //half matrix
      for (unsigned int y = startpos; y < startpos+pos; ++y) { 
	main = RFMatrixRow[y];
	ab = ab_array[pos];
	ac = ab_array[y];
	value = main;
	if ( (distance_measure != 1) && (distance_measure != 19) )
	  calculatew(main,ab,ac,distance_measure, value);
	if (!PRINT_L){
	  if (PRINT){
	      rs << value << " ";
	  }
	  else{
	    if (value != NAN){
	      sum += value;
	      count++;
	    }
	  }
	}
	else{
	  rs << pos << "." << y << " " << value << "\n";
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
    } //end else half matrix
   
    if (NUM_TREES > 100){
      if (pos % per_10 == 0){
	per_place = per_counter*10;
	fprintf(stderr, "%d%% Complete\n", per_place);
	per_counter++;
      }
    }
 }
  if (PRINT|| PRINT_L){
    fout.close();
  }
  average = (double)sum/count;
  cout << average << endl;

  fprintf(stderr, "Finishing Traversal\n");  
  fprintf(stderr, "Number of comparisons: %lu\n", num_comparisons);
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
}

void process_trees_weighted(string infilename,  string qfilename, unsigned int distance_option, string outfilename, float * ab_array){
  unsigned int hashtable_pos = 0;
  ofstream fout; //we will use this if we need to print out the hash table
  unsigned long long M1=0;
  unsigned long long M2=0;
  FILE *fp;
  unsigned long uBID = 0;
  tree_counter = 0; //this determines if the trees are multifurcating or binary
  unsigned int threshold = NUM_TREES - (NUM_TREES*P)/100;

  struct timeval hash_start;
  struct timeval hash_end;
  struct timeval ragged_start;
  struct timeval ragged_end;
  struct timeval index_start;
  struct timeval index_end;
  unsigned int hashtablesize = 0;
  //counter information for amazing 
  gettimeofday(&index_start, NULL);

  //allocate helper data structures
  if (QOPT)
    hashtablesize = TOTAL_TREES;
  else
    hashtablesize = NUM_TREES;

  //allocate helper arrays (inverted index)
  windex = (struct weighted_element **)malloc(hashtablesize*sizeof(struct weighted_element *));
  helper_sizes = (unsigned int *)malloc(hashtablesize*sizeof(unsigned int));
  for (unsigned int i =0; i < hashtablesize; i++) { 
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
    mc_procedure(M1, M2, 0, NUM_TREES, fp, lm, ROOTED);
  }
  else {
    lv_procedure(M1, M2, 0, NUM_TREES, fp, lm, ROOTED);
  }

  fclose(fp);

  if (QOPT){
    fp = fopen(qfilename.c_str(), "r");
    if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}
    if (MC){ //do MC procedure
      mc_procedure(M1, M2, NUM_TREES, TOTAL_TREES, fp, lm, ROOTED);
    }
    else{
      lv_procedure(M1, M2, NUM_TREES, TOTAL_TREES, fp, lm, ROOTED);
    }
  }
  gettimeofday(&hash_end, NULL);
  double hash_time = hash_end.tv_sec - hash_start.tv_sec + (hash_end.tv_usec - hash_start.tv_usec) / 1.e6;
  fprintf(stderr, "Hash table Creation Time: %g s\n", hash_time);

  if (QOPT)
    NUM_TREES = TOTAL_TREES;

  //count up the number of unique bipartitions -- note, this function also tells us what dimensions we need for our inverted index
  unsigned int unique_bipartitions = 0;
  calculate_unique_bipartitions(threshold, unique_bipartitions);
  //allocate helper data structures
  gettimeofday(&index_start, NULL);
  //now, allocate helper
  unsigned int newval = 0;
  for (unsigned int i = 0; i < NUM_TREES; ++i) { 
    newval = helper_sizes[i];
    windex[i] = (struct weighted_element*)malloc(newval*sizeof(struct weighted_element));
    helper_sizes[i]=0;
  } 
  gettimeofday(&index_end, NULL);
  temp_time = index_end.tv_sec - index_start.tv_sec + (index_end.tv_usec - index_start.tv_usec) / 1.e6;
  index_time+=temp_time;
  
  weighted_hashtable = (struct weighted_element **)malloc(unique_bipartitions*sizeof(struct weighted_element *));
  hash_lengths = (unsigned int*)malloc(unique_bipartitions*sizeof(unsigned int)); 
  for( unsigned int i = 0; i < unique_bipartitions; ++i) {
    hash_lengths[i] = 0;
  }

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
	weighted_hashtable[hashtable_pos] = (struct weighted_element *)malloc(sizeTreeIdx*sizeof(struct weighted_element));
	hash_lengths[hashtable_pos] = sizeTreeIdx;
	for (unsigned int j=0; j<sizeTreeIdx; ++j) {	
	  unsigned int temp  = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
	  float dist  = vvec_hashrf._hashtab2[hti][i]._vec_dist[j];
	  float val = dist*dist;
	  if (distance_option == 23)
	    val = dist;
	  ab_array[temp] += val;
	  weighted_hashtable[hashtable_pos][j].t_id = temp;
	  weighted_hashtable[hashtable_pos][j].weight = dist;
	  gettimeofday(&index_start, NULL);
	  amazing_loc = helper_sizes[temp];
	  windex[temp][amazing_loc].t_id = hashtable_pos;
	  windex[temp][amazing_loc].weight = dist;
	  helper_sizes[temp]++;
	  gettimeofday(&index_end, NULL);
	  temp_time = index_end.tv_sec - index_start.tv_sec + (index_end.tv_usec - index_start.tv_usec) / 1.e6;
	  index_time+=temp_time;
	} //for everything in sizeTreeIdx
	hashtable_pos++;
      } //for everything in sizeVec
    } //if sizevec
  } //end hashtable creation
 
  assert(uBID == unique_bipartitions);
  gettimeofday(&ragged_end, NULL);
  double ragged_time = ragged_end.tv_sec - ragged_start.tv_sec + (ragged_end.tv_usec - ragged_start.tv_usec) / 1.e6;

  fprintf(stderr, "\nRagged Array Creation Time: %g s\n", ragged_time);

  fprintf(stderr, "UNIQUE BIPARTITIONS: %d\n", unique_bipartitions);
  NUMBIPART = unique_bipartitions;
  
  //determine if trees are multifurcating
  if (tree_counter < ((NUM_TAXA-3)*NUM_TREES)) { 
    fprintf(stderr, "Trees are Multifurcating\n");
    if (!MULTIFURCATING){
      MULTIFURCATING = true;
      fprintf(stderr, "\nWARNING: Detected multifurcated trees!\n Trees will be treated as such through remainder of computation\n");
    }
  }
  else { 
    fprintf(stderr, "Trees are Binary\n");
  } 
  fprintf(stderr,"\nInverted Index Creation Time: %g s\n", index_time);

  //print out data structures
  /*
  cout << "printing out hashtable" << endl;
  for (unsigned int i  = 0; i < unique_bipartitions; i++){
    unsigned int lengths = hash_lengths[i];
    for (unsigned int j = 0; j < lengths; j++){
      cout << weighted_hashtable[i][j].t_id << " ";
    } 
    cout << endl;
  }

  cout << "associated weights: " << endl;
  for (unsigned int i  = 0; i < unique_bipartitions; i++){
    unsigned int lengths = hash_lengths[i];
    for (unsigned int j = 0; j < lengths; j++){
      cout << weighted_hashtable[i][j].weight << " ";
    } 
    cout << endl;
  }

  cout << "printing out inverted index: " << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++){
    unsigned int lengths = helper_sizes[i];
    for (unsigned int j = 0; j < lengths; j++){
      cout << index[i][j].t_id << " ";
    }
    cout << endl;
  }

  cout << "printing out associated weights: " << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++){
    unsigned int lengths = helper_sizes[i];
    for (unsigned int j = 0; j < lengths; j++){
      cout << index[i][j].weight << " ";
    }
    cout << endl;
  }

  cout << "printing out ab_array: " << endl;
  for (unsigned int i = 0; i < NUM_TREES; i++)
    cout << ab_array[i] << " ";
    cout << endl;*/
  //weighted_row_distance(weighted_hashtable, index, ab_array, distance_option, outfilename);
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
  
  if (MULTIFURCATING)
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
	  if (TRIVIAL){
	    a += NUM_TAXA; 
	    b = n_biparts[pos] + NUM_TAXA - a; //b 
	    c = n_biparts[y] + NUM_TAXA - a; //c
	  }
	  else{
	    b = n_biparts[pos] - a; //b 
	    c = n_biparts[y] - a; //c
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
  fprintf(stderr, "Number of comparisons: %lu\n", num_comparisons);
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

void handle_newick_error(int err){
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

void calculate_unique_bipartitions(unsigned int threshold, unsigned int & unique_bipartitions){
  struct timeval index_start;
  struct timeval index_end;
  double temp_time = 0;

  vector<bool> check(NUM_TREES, false);
  unique_bipartitions = 0;
  for (unsigned int hti=0; hti<vvec_hashrf._hashtab2.size(); ++hti) {
    unsigned int sizeVec = vvec_hashrf._hashtab2[hti].size(); 
    if (sizeVec) {
      unique_bipartitions += sizeVec;      
      gettimeofday(&index_start, NULL);
      for (unsigned int i=0; i<sizeVec; ++i) {
	unsigned int sizeTreeIdx = vvec_hashrf._hashtab2[hti][i]._vec_treeidx.size();	  
	if (!WEIGHTED){
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
	}
	else{
	  for (unsigned int j=0; j<sizeTreeIdx; j++) {	
	    unsigned int temp  = vvec_hashrf._hashtab2[hti][i]._vec_treeidx[j];
	    helper_sizes[temp]++;
	  }
	}
      } //end for everything in sizeVec
      gettimeofday(&index_end, NULL);
      temp_time = index_end.tv_sec - index_start.tv_sec + (index_end.tv_usec - index_start.tv_usec) / 1.e6;
      index_time+=temp_time;  
    } //end if sizevec
  } //end index creation
  
  NUMBIPART = unique_bipartitions;
}


void mc_procedure(unsigned long long M1, unsigned long long M2, unsigned int start, unsigned int end, FILE * fp, LabelMap &lm, bool rooted){
  NEWICKTREE *newickTree;
  int err;
  if (rooted){
    for (unsigned int treeIdx=start; treeIdx<end; ++treeIdx) {
      newickTree = loadnewicktree2(fp, &err);
      if (!newickTree) {
	handle_newick_error(err);
      }
      else {
	unsigned int numBitstr=0;
	dfs_compute_hash_mc(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, NULL);
	if (MAXVAL == 0){
	    MAXVAL = numBitstr -1;
	}
	killnewicktree(newickTree);
      }
    } //end for all trees
  }
  else{
   for (unsigned int treeIdx=start; treeIdx<end; ++treeIdx) {
      newickTree = loadnewicktree2(fp, &err);
      if (!newickTree) {
      handle_newick_error(err);
      }
      else {
	unsigned int numBitstr=0;
	bool processed = false;
	dfs_compute_hash_unrooted_mc(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, processed, NULL);
	if (MAXVAL == 0){
	  MAXVAL = NUM_TAXA - 3;
	}
	killnewicktree(newickTree);
      }
   }
  }
}
 
void lv_procedure(unsigned long long M1, unsigned long long M2, unsigned int start, unsigned int end, FILE * fp, LabelMap & lm, bool rooted){
  NEWICKTREE * newickTree;
  int err;
  if (rooted){
    for (unsigned int treeIdx=start; treeIdx<end; ++treeIdx) {
      newickTree = loadnewicktree2(fp, &err);
      if (!newickTree) {
	handle_newick_error(err);
      }
      else {
	unsigned int numBitstr=0;
	dfs_compute_hash(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, NULL);
	if (MAXVAL == 0)
	  MAXVAL = numBitstr -1;
	killnewicktree(newickTree);
      }
    }
  }
  else{
   for (unsigned int treeIdx=start; treeIdx<end; ++treeIdx) {
      newickTree = loadnewicktree2(fp, &err);
      if (!newickTree) {
	handle_newick_error(err);
      }
      else {
	unsigned int numBitstr=0;
	bool processed = false;
	dfs_compute_hash_unrooted2(newickTree->root, lm, vvec_hashrf, treeIdx, numBitstr, M1, M2, processed, NULL);
	if (MAXVAL == 0)
	  MAXVAL = NUM_TAXA - 3;
	killnewicktree(newickTree);
      }
    }
  }
}

void process_trees(string infilename, 
   string qfilename,
   unsigned int & hashtable_length, 
   unsigned int * row){
  unsigned int hashtable_pos = 0;
   unsigned int hashtable_length_correction = 0;
  ofstream fout; //we will use this if we need to print out the hash table
  unsigned long long M1=0;
  unsigned long long M2=0;
  FILE *fp;
  unsigned long uBID = 0;
  tree_counter = 0; //this determines if the trees are multifurcating or binary

  struct timeval hash_start;
  struct timeval hash_end;
  struct timeval ragged_start;
  struct timeval ragged_end;
  struct timeval index_start;
  struct timeval index_end;

  //counter information for amazing 
  gettimeofday(&index_start, NULL);
  //allocate helper data structures
  unsigned int hashtablesize = 0;
  if (QOPT)
    hashtablesize = TOTAL_TREES;
  else
    hashtablesize = NUM_TREES;

  //allocate helper arrays (inverted index)
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
    mc_procedure(M1, M2, 0, NUM_TREES, fp, lm, ROOTED);
  }
  else {
    lv_procedure(M1, M2, 0, NUM_TREES, fp, lm, ROOTED);
  }

  fclose(fp);

  if (QOPT){
    fp = fopen(qfilename.c_str(), "r");
    if (!fp) {cout << "Fatal error: file open error\n";  exit(0);}
    if (MC){ //do MC procedure
      mc_procedure(M1, M2, NUM_TREES, TOTAL_TREES, fp, lm, ROOTED);
    }
    else{
      lv_procedure(M1, M2, NUM_TREES, TOTAL_TREES, fp, lm, ROOTED);
    }
  }
  gettimeofday(&hash_end, NULL);
  double hash_time = hash_end.tv_sec - hash_start.tv_sec + (hash_end.tv_usec - hash_start.tv_usec) / 1.e6;
  fprintf(stderr, "Hash table Creation Time: %g s\n", hash_time);

  if (QOPT)
    NUM_TREES = TOTAL_TREES;

  vector<bool> check(NUM_TREES, false);
  unsigned int threshold = NUM_TREES - (NUM_TREES*P)/100;
  n_biparts = (unsigned int *)malloc(NUM_TREES * sizeof(unsigned int));
  for (unsigned int i = 0; i < NUM_TREES; i++)
    n_biparts[i] = 0;

  //count up the number of unique bipartitions
  unsigned int unique_bipartitions = 0;
  calculate_unique_bipartitions(threshold, unique_bipartitions);
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

  fprintf(stderr, "UNIQUE BIPARTITIONS: %d\n", unique_bipartitions);
  NUMBIPART = unique_bipartitions;
  
  //determine if trees are multifurcating
  if (tree_counter < ((NUM_TAXA-3)*NUM_TREES)) { 
    fprintf(stderr, "Trees are Multifurcating\n");
    if (!MULTIFURCATING){
      MULTIFURCATING = true;
      fprintf(stderr, "\nWARNING: Detected multifurcated trees!\n Trees will be treated as such through remainder of computation\n");
    }
  }
  else { 
    fprintf(stderr, "Trees are Binary\n");
    free(n_biparts);
  }
 
  fprintf(stderr,"\nInverted Index Creation Time: %g s\n", index_time);
}
