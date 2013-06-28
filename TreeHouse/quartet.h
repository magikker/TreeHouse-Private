#ifndef _QUARTET_H
#define _QUARTET_H

#include <vector>
#include "THGlobals.h"
#include "pql.h"
#include "AnalysisFunctions.h"
#include <unordered_set>
#include <set>
#include <algorithm>
#include <iterator>
#include <iomanip>
#include <string>
#include <sstream>
using namespace std;


void calculateTrivialBipartitions();

typedef pair<int, int> iPair;
typedef pair< set<iPair>, set<iPair> > qPair;

void printQPair(qPair a);
qPair mergeQPair(qPair a, qPair b);
bool isQPairEmpty(qPair a);

class quartet{ //a class to store quartets. I'm open to ideas about how to make this better...
	private:	
		iPair A; //one side of the quartet
		iPair B; //the other side
	public:
		quartet(int, int, int, int);
		quartet(iPair, iPair);
		quartet();
		iPair getA() const;
		iPair getB() const;
		void setA(int, int);
		void setB(int, int);
		void print() const;
		//bool operator<(const quartet& other) const;
		//bool operator==(const quartet& other);

};

  bool operator<(const quartet& one, const quartet& two);
  bool operator<=(const quartet& one, const quartet& two);
  bool operator>(const quartet& one, const quartet& two);
  bool operator>=(const quartet& one, const quartet& two);
  bool operator==(const quartet& one, const quartet& two);
  //void operator=(quartet& one, const quartet& two);


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
	void printMetaData();
	
	void invertB(); //inverts vector B, in case we need to do that in calculate
	void calculate(); //calculates and sets meta data vectors
				  //also inverts the second vector if necessary  	
	boost::dynamic_bitset<> getA(); //returns first bipartiton
	boost::dynamic_bitset<> getB(); //returns second bipartition
        bPair(boost::dynamic_bitset<> a, boost::dynamic_bitset<> b);
	
};


char isQuartetImplied(quartet x, unsigned int bipart);
char isQuartetImplied(quartet x, set<unsigned int> biparts);

void printQuartets(vector<quartet>);
vector<quartet> generateQuartetsFromBipart(int);
set<quartet> generateQuartetsFromBipartSet(int);
vector<quartet> generateQuartetsFromBipart(vector<bool>);
void insertQuartetsFromBipart(int b, set<quartet> &setty);

vector<quartet> generateQuartetsFromOnesZeros(vector<int> ones, vector<int> zeros);
set<quartet> generateQuartetsFromOnesZerosSet(vector<int> ones, vector<int> zeros);


vector<iPair> nChooseTwo(vector<int> in);
vector<iPair> matchTaxa(vector<int> a, vector<int> b);


void fillDifferenceVectors(boost::dynamic_bitset<> v1, boost::dynamic_bitset<> v2, vector<int> &sameOnes, vector<int> &sameZeros, vector<int> &uniqueOnesA, vector<int> &uniqueOnesB, vector<int> &uniqueZerosA, vector<int> &uniqueZerosB, vector<int> &OnesA, vector<int> &OnesB, vector<int> &ZerosA, vector<int> &ZerosB);

void addMatchedPairs(vector<iPair> p1, vector<iPair> p2, vector<quartet> &retVec);
vector<quartet> generateDifferentQuartets(int bipartA, int bipartB);
vector<quartet> generateSameQuartets(int bipartA, int bipartB);


set<quartet> generateConflictingQuartets(int bipart1, int bipart2);
set<quartet> generateConflictingQuartets2(int bipart1, int bipart2);
qPair generateConflictingQuartets3(int bipart1, int bipart2);
unsigned int numConflictingQuartets(int bipart1, int bipart2);
set<quartet> generateConflictingQuartetsGroup(int bipart1, set<unsigned int> bipart2);

unsigned int getNumQuartets(int);
unsigned int getNumQuartets(vector<bool> bipart);
unsigned int getNumQuartets(boost::dynamic_bitset<> bipart);
unsigned int getNumQuartets(int nOnes, int nZeros);

unsigned int getNumDifferentQuartets(int a, int b);
unsigned int getNumSameQuartets(int a, int b);
unsigned int getNumSameQuartets(boost::dynamic_bitset<> a, boost::dynamic_bitset<> b);
unsigned int getNumDifferentQuartets(boost::dynamic_bitset<> a, boost::dynamic_bitset<> b);

set<quartet> generateDifferentQuartetsFromTrees(int a, int b);
set<quartet> generateDifferentQuartetsFromTrees2(int a, int b);
set<quartet> generateDifferentQuartetsFromTrees3(int a, int b);
set<quartet> generateDifferentQuartetsFromTrees4(int a, int b);
set<quartet> generateSameQuartetsFromTrees(int a, int b);
set<quartet> generateQuartetsFromTree(int t);
set<quartet> generateQuartetsFromBiparts(set<unsigned int> s);

unsigned int quartet_distance(int tree1, int tree2);

void shared_quartets_strict(set<unsigned int> trees);
void shared_quartets_majority(set<unsigned int> trees);

void ktetAnalysis(int t);
void quartetAnalysis(int, int); //not fully implemented
void printSet(set<quartet> s);
void TESTSTUFF();
void testGenerateConflictingQuartets();
void testNumConflictingQuartets();
void testGenerateDifferentQuartets();
void testOperatorsForQuartets();
void testGeenerateDifferentQuartetsFromTrees();

#endif
