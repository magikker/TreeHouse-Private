#ifndef _QUARTET_H
#define _QUARTET_H

#include <vector>
#include "THGlobals.h"
#include "pql.h"
using namespace std;

typedef pair<int, int> iPair;

class quartet{ //a class to store quartets. I'm open to ideas about how to make this better...
	private:	
		iPair A; //one side of the quartet
		iPair B; //the other side
	public:
		quartet(int, int, int, int);
		quartet(iPair, iPair);
		quartet();
		iPair getA();
		iPair getB();
		void setA(int, int);
		void setB(int, int);
		void print();

};

void printQuartets(vector<quartet>);
vector<quartet> generateQuartetsFromBipart(int);
vector<quartet> generateQuartetsFromBipart(vector<bool>);

vector<quartet> generateQuartetsFromOnesZeros(vector<int> ones, vector<int> zeros);

vector<iPair> nChooseTwo(vector<int> in);
vector<iPair> matchTaxa(vector<int> a, vector<int> b);

void fillDifferenceVectors(vector<char> v1, vector<char> v2, vector<int> &sameOnes, vector<int> &sameZeros, vector<int> &uniqueOnesA, vector<int> &uniqueOnesB, vector<int> &uniqueZerosA, vector<int> &uniqueZerosB, vector<int> &OnesA, vector<int> &OnesB, vector<int> &ZerosA, vector<int> &ZerosB);
void addMatchedPairs(vector<iPair> p1, vector<iPair> p2, vector<quartet> &retVec);
vector<quartet> generateDifferentQuartets(int bipartA, int bipartB);
vector<quartet> generateSameQuartets(int bipartA, int bipartB);

unsigned int getNumQuartets(int);
unsigned int getNumQuartets(vector<bool> bipart);
unsigned int getNumQuartets(vector<char> bipart);
unsigned int getNumQuartets(int nOnes, int nZeros);

unsigned int getNumDifferentQuartets(int a, int b);
unsigned int getNumDifferentQuartets(vector<char> a, vector<char> b);

void TESTSTUFF();

#endif
