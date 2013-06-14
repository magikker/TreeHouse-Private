#ifndef _QUARTET_H
#define _QUARTET_H

#include <vector>
#include "THGlobals.h"
#include "pql.h"
#include <unordered_set>
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

class bPair{ //holds all relevant info about quartets on a pair of bipartitions
  private:
 	boost::dynamic_bitset<> A;
	boost::dynamic_bitset<> B;
  public:
	//meta data about vectors A and B
	vector<int> sameOnes;
	vector<int> sameZeros;
	vector<int> uniqueOnesA;
	vector<int> uniqueZerosA;
	vector<int> uniqueOnesB;
	vector<int> uniqueZerosB;
	vector<int> OnesA;
	vector<int> ZerosA;
	vector<int> OnesB;
	vector<int> ZerosB;
	void printVectors();
	void printMetaData();
	
	void invertB(); //inverts vector B, in case we need to do that in calculate
	void calculate(); //calculates and sets meta data vectors
				  //also inverts the second vector if necessary  	
	boost::dynamic_bitset<> getA(); //returns first bipartiton
	boost::dynamic_bitset<> getB(); //returns second bipartition
        bPair(boost::dynamic_bitset<> a, boost::dynamic_bitset<> b);
	
};

void printQuartets(vector<quartet>);
vector<quartet> generateQuartetsFromBipart(int);
vector<quartet> generateQuartetsFromBipart(vector<bool>);

vector<quartet> generateQuartetsFromOnesZeros(vector<int> ones, vector<int> zeros);

vector<iPair> nChooseTwo(vector<int> in);
vector<iPair> matchTaxa(vector<int> a, vector<int> b);

void fillDifferenceVectors(boost::dynamic_bitset<> v1, boost::dynamic_bitset<> v2, vector<int> &sameOnes, vector<int> &sameZeros, vector<int> &uniqueOnesA, vector<int> &uniqueOnesB, vector<int> &uniqueZerosA, vector<int> &uniqueZerosB, vector<int> &OnesA, vector<int> &OnesB, vector<int> &ZerosA, vector<int> &ZerosB);
void addMatchedPairs(vector<iPair> p1, vector<iPair> p2, vector<quartet> &retVec);
vector<quartet> generateDifferentQuartets(int bipartA, int bipartB);
vector<quartet> generateSameQuartets(int bipartA, int bipartB);

unsigned int getNumQuartets(int);
unsigned int getNumQuartets(vector<bool> bipart);
unsigned int getNumQuartets(boost::dynamic_bitset<> bipart);
unsigned int getNumQuartets(int nOnes, int nZeros);

unsigned int getNumDifferentQuartets(int a, int b);
unsigned int getNumSameQuartets(int a, int b);
unsigned int getNumSameQuartets(boost::dynamic_bitset<> a, boost::dynamic_bitset<> b);
unsigned int getNumDifferentQuartets(boost::dynamic_bitset<> a, boost::dynamic_bitset<> b);



void TESTSTUFF();

#endif
