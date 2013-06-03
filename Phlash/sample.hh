#ifndef SAMPLE_HH_
#define SAMPLE_HH_

#include <string>
#include <iostream>
using namespace std;

void setup_sampling(string infile, unsigned int dist, string outfile, unsigned int nsample, bool sort_sampling);
void weighted_sampling(string infile, unsigned int dist, string outfile, unsigned int nsamples, unsigned int sample_size, bool sort_sampling);
void unweighted_sampling(string infile, unsigned int dist, string outfile, unsigned int nsamples, unsigned int sample_size, bool sort_sampling);
#endif
