#include "quartet.h"

using namespace std;

set<int> trivial_bipartitions;

quartet::quartet(int a0, int a1, int b0, int b1){

A = make_pair(a0, a1);
B = make_pair(b0, b1);
//  a.first = a0;
//  a.second = a1;
//  b.first = b0;
//  b.second = b1;
}

quartet::quartet(){
 A = make_pair(-1, -1); //assign values to -1 to indicate there is not yet data
 B = make_pair(-1, -1);
}

quartet::quartet(iPair a, iPair b){
	//sort both quartets so they can be more easily compared
	int temp;
	if(a.second < a.first){
		temp = a.first;
		a.first = a.second;
		a.second = a.first;
	}
	if(b.second < b.first){
		temp = b.first;
		b.first = b.second;
		b.second = b.first;
	}
	if(b.first < a.first || (b.first==a.first && b.second<a.second)){
		A = b;
		B = a;
	}
	else{
		A = a;
		B = b;
	}

}

bool operator<(const quartet& one, const quartet& two){
  if(one.getA() < two.getA()){
	return true;
	}
  else{
	return(one.getA()==two.getA() && one.getB()<two.getB());
	}

}

bool operator<=(const quartet& one, const quartet& two){
  if(one.getA() < two.getA()){
	return true;
	}
  else{
	return(one.getA()==two.getA() && one.getB()<=two.getB());
	}

}

bool operator>(const quartet& one, const quartet& two){
  return !(one<=two);
}

bool operator>=(const quartet& one, const quartet& two){
  return !(one<two);
}

bool operator==(const quartet& one, const quartet& two){
  return(one.getA()==two.getA() && one.getB()==two.getB());
}

/*
void operator=(quartet& one, const quartet& two) const{
 one.setA(two.getA());
 one.setB(two.getB());
}*/

/*
bool quartet::operator<(const quartet& other) const{
  if(A < other.getA()){
	return true;
	}
  else{
	return(A==other.getA() && B<other.getB());
	}
}

bool quartet::operator==(const quartet& other){
  return(A==other.getA() && B==other.getB());
}
*/


iPair quartet::getA() const{
 return A;
}

iPair quartet::getB() const{
 return B;
}

void quartet::setA(int x, int y){
  A = make_pair(x, y);
} 

void quartet::setB(int x, int y){
  B = make_pair(x, y);
}

void quartet::print() const{
 int vals[4];
 vals[0] = A.first; vals[1] = A.second;
 vals[2] = B.first; vals[3] = B.second;
 for(int i = 0; i < 4; i++){
	if(vals[i] == -1){
		cout << "? ";
		}
	else{
		cout << vals[i];
		if(i!=3) {cout << " ";}
		}
	if(i==1){cout << "|";}
	}
 cout << endl;
}

void printQuartets(vector<quartet> input){
 for(unsigned int i = 0; i < input.size(); i++){
	cout << i+1 << ": ";
	input.at(i).print();
	}
}

bPair::bPair(boost::dynamic_bitset<> a, boost::dynamic_bitset<> b){
 if(a.size()!=b.size()){
	cerr << "bPair constructor error! You gave vectors of different size! Unpredictable behavior may ensue!";
	//TODO- pad smaller vector with 0s
	return;
	}
A = a;
B = b;
calculate();
}

void bPair::invertB(){
  B.flip();
  cout << "bPair: B inverted\n";


}

boost::dynamic_bitset<> bPair::getA(){
return A;
}

boost::dynamic_bitset<> bPair::getB(){
return B;
}

void bPair::calculate(){
//this finds out if we need to invert the second bipartition. Also, it fills all of the meta data
//first, call getNumQuartets on sameOnes and sameZeros. 

//Then, see if inverting one of the vectors leads to a greater number of similar quartets
//if it is the case, we will invert the vector and recalculate all of the meta vectors
  for(unsigned int i = 0; i < A.size(); i++){
	if(A[i]==1){
		OnesA.push_back(i);
		if(B[i]==1) {sameOnes.push_back(i); OnesB.push_back(i);}
		else{uniqueOnesA.push_back(i); uniqueZerosB.push_back(i); ZerosB.push_back(i);}
		}
	else{//A[i]==0
		ZerosA.push_back(i);
		if(B[i]==0){sameZeros.push_back(i); ZerosB.push_back(i);}
		else{uniqueZerosA.push_back(i); uniqueOnesB.push_back(i); OnesB.push_back(i);}
		}
	}
  unsigned int similarQuartetsWithoutInverting = getNumQuartets(sameOnes.size(), sameZeros.size());
  unsigned int similarQuartetsInverting = getNumQuartets(OnesA.size()-sameOnes.size(), ZerosA.size()- sameZeros.size());
  //cout << "bPair calculate: without inverting, number of similar quartets is " << similarQuartetsWithoutInverting << endl;
 // cout << "IF we invert, number of similar quartets is " << similarQuartetsInverting << endl;

  if(similarQuartetsInverting > similarQuartetsWithoutInverting){
   //if we get more similar vectors, invertB and recalculate everything
   //TODO- possible optimization- don't iterate through loops again, just calculate using existing metadata
	  invertB();
	  vector<int> temp, temp2;  
	//onesA and zerosA unchanged
          //unique and same switch around

	/* Table showing vectors before and after
	BEFORE		AFTER	
	sameOnes  	uniqueOnesA		
	sameZeros	uniqueZerosA
	uniqueOnesA	sameOnes
	uniqueZerosA	sameZeros
	onesB		ZerosB	
	zerosB		OnesB
	uniqueOnesB	sameZeros
	uniqueZerosB	sameOnes
	OnesA, ZerosA	unchanged
*/
	  temp = uniqueOnesA;
	  temp2 = uniqueZerosA;
 	  uniqueOnesA = sameOnes;
	  uniqueZerosB = sameOnes;
	  uniqueOnesB = sameZeros;
	  uniqueZerosA = sameZeros;
	  sameOnes = temp;
	  sameZeros = temp2;
	  temp = OnesB;
	  OnesB = ZerosB;
	  ZerosB = temp;
	}	
}

void bPair::printMetaData(){
	cout << "printing vectors: Same ones ";
	printVectorCompact(sameOnes);
	cout << endl << "printing vectors: unique ones A";
	printVectorCompact(uniqueOnesA);
	cout << endl << "printing vectors: unique ones B ";
	printVectorCompact(uniqueOnesB);
	cout << endl << "printing vectors: Same zeros ";
	printVectorCompact(sameZeros);
	cout << endl << "printing vectors: uniqueZerosA ";
	printVectorCompact(uniqueZerosA);
	cout << endl << "printing vectors: uniqueZerosB ";
	printVectorCompact(uniqueZerosB);
	cout << endl << "printing vectors: OnesA ";
	printVectorCompact(OnesA);
	cout << endl << "printing vectors: OnesB ";
	printVectorCompact(OnesB);
	cout << endl << "printing vectors: ZerosA ";
	printVectorCompact(ZerosA);
	cout << endl << "printing vectors: ZerosB ";
	printVectorCompact(ZerosB);
}

vector<quartet> generateQuartetsFromBipart(vector<bool> b){
  vector<int> ones, zeros; 
  for(unsigned int i = 0; i < b.size(); i++){
	if(b.at(i)==1){
		ones.push_back(i);
		}
	else{
		zeros.push_back(i);
		}
	}
  return generateQuartetsFromOnesZeros(ones, zeros);
}


//fillDifferenceVectors was used for generateDifferentQuartets, but has been
//	replaced by the bPair class
void fillDifferenceVectors(boost::dynamic_bitset<> v1, boost::dynamic_bitset<> v2, vector<int> &sameOnes, vector<int> &sameZeros, vector<int> &uniqueOnesA, vector<int> &uniqueOnesB, vector<int> &uniqueZerosA, vector<int> &uniqueZerosB, vector<int> &OnesA, vector<int> &OnesB, vector<int> &ZerosA, vector<int> &ZerosB){
 
  //iterate through vectors, fill appropriate vectors
    if(v1.size()!=v2.size()){
	cerr << "fillDifferenceVectors Error!: the vectors are of different sizes!";
	return;
	}
  for(unsigned int i = 0; i < v1.size(); i++){
	if(v1[i]==1){
		OnesA.push_back(i);
		if(v2[i]==1) {sameOnes.push_back(i); OnesB.push_back(i);}
		else{uniqueOnesA.push_back(i); uniqueZerosB.push_back(i); ZerosB.push_back(i);}
		}
	else{//v1[i]==0
		ZerosA.push_back(i);
		if(v2[i]==0){sameZeros.push_back(i); ZerosB.push_back(i);}
		else{uniqueZerosA.push_back(i); uniqueOnesB.push_back(i); OnesB.push_back(i);}
		}
	}
}

void addMatchedPairs(vector<iPair> p1, vector<iPair> p2, vector<quartet> &retVec){
  for(vector<iPair>::iterator i = p1.begin(); i!=p1.end(); i++){
	for(vector<iPair>::iterator j = p2.begin(); j!=p2.end(); j++){
			retVec.push_back(quartet(*i, *j));
			}
		}
}

vector<quartet> generateDifferentQuartets(int bipartA, int bipartB){
  vector<quartet> retVec;

  //for algorithm description, see comment at the bottom of the file


 boost::dynamic_bitset<> A = biparttable.non_trunc_bitstring(bipartA);
 boost::dynamic_bitset<> B = biparttable.non_trunc_bitstring(bipartB);
 //optimization note- instead of calling sharedOnes then sharedZeros, we could make a function
	//that passes once through the vectors and puts elements in sharedOnes or sharedZero (or neither)
 
 bPair b = bPair(A, B);


  //ALL QUARTETS IN VECTOR A BUT NOT IN VECTOR B
  //generate all iPairs of nChoose2 from uniqueOnesA, then match each uniqueOnesA with each sameOnes 
  //pairsA is all of the iPairs that can be generated from 1s in A that AREN'T shared
 
 vector<iPair> uniquePairs1A = nChooseTwo(b.uniqueOnesA);
  vector<iPair> matchedOnesA = matchTaxa(b.uniqueOnesA, b.sameOnes);
  uniquePairs1A.insert(uniquePairs1A.end(), matchedOnesA.begin(), matchedOnesA.end());
  //now, combine pairsOneA with all pairs of zeros from A
  vector<iPair> pairs0A = nChooseTwo(b.ZerosA);
  addMatchedPairs(uniquePairs1A, pairs0A, retVec);
  //now we need to match all permutations of shared 1s with (unique 0s + matched(unique0s, shared0s)
  vector<iPair> uniquePairs0A = nChooseTwo(b.uniqueZerosA);
  vector<iPair> matchedZerosA = matchTaxa(b.uniqueZerosA, b.sameZeros); 
  uniquePairs0A.insert(uniquePairs0A.end(), matchedZerosA.begin(), matchedZerosA.end());
  addMatchedPairs(nChooseTwo(b.sameOnes), uniquePairs0A, retVec);

  //ALL QUARTETS IN VECTOR B BUT NOT IN VECTOR A

 vector<iPair> uniquePairs1B = nChooseTwo(b.uniqueOnesB);
  vector<iPair> matchedOnesB = matchTaxa(b.uniqueOnesB, b.sameOnes);
  uniquePairs1B.insert(uniquePairs1B.end(), matchedOnesB.begin(), matchedOnesB.end());
  //now, combine pairsOneB with all pairs of zeros from B
  vector<iPair> pairs0B = nChooseTwo(b.ZerosB);
  addMatchedPairs(uniquePairs1B, pairs0B, retVec);
  //now we need to match all permutations of shared 1s with (unique 0s + matched(unique0s, shared0s)
  vector<iPair> uniquePairs0B = nChooseTwo(b.uniqueZerosB);
  vector<iPair> matchedZerosB = matchTaxa(b.uniqueZerosB, b.sameZeros); 
  uniquePairs0B.insert(uniquePairs0B.end(), matchedZerosB.begin(), matchedZerosB.end());
  addMatchedPairs(nChooseTwo(b.sameOnes), uniquePairs0B, retVec);


  return retVec;
}


//Doesn't always return a value
vector<iPair> nChooseTwo(vector<int> in){
//generates all pairs of nChoose2 from input vector
//size of return vector is   in.size() choose 2
 vector<iPair> ret;
 if(in.size()<2){ //make sure our input is big enough
	return ret;
	}
 for(vector<int>::iterator i = in.begin(); i!=in.end()-1; i++){
	for(vector<int>::iterator j = i+1; j!=in.end(); j++){
		ret.push_back(make_pair(*i, *j));	
		}
	}
  return ret;
}

vector<iPair> matchTaxa(vector<int> A, vector<int> B){
//makes all iPairs by matching one taxon from A and one taxon from B
//size of return vector is A.size()*B.size()
 vector<iPair> ret;
  for(vector<int>::iterator i = A.begin(); i!=A.end(); i++){
	for(vector<int>::iterator j=B.begin(); j!=B.end(); j++){
		ret.push_back(make_pair(*i, *j));	
		}
	}
  return ret;
}

vector<quartet> generateSameQuartets(int bipartA, int bipartB){
   boost::dynamic_bitset<> A = biparttable.non_trunc_bitstring(bipartA);
   boost::dynamic_bitset<> B = biparttable.non_trunc_bitstring(bipartB);
   boost::dynamic_bitset<> C;
   boost::dynamic_bitset<> D;
   
  C = A & B;
  D = A ^ B;
  
  vector<int> ones;
  vector<int> zeros;
 
 for(unsigned int i = 0; i < A.size(); i++){
	if (C[i] == true){
		ones.push_back(i);
	}
	if (D[i] == true){
		zeros.push_back(i);
	}
 }
 
 return generateQuartetsFromOnesZeros(ones, zeros);
}


set<quartet> generateQuartetsFromBipartSet(int b){
  vector<int> ones;
  vector<int> zeros;
  for(int i = 0; i < biparttable.BipartitionTable.at(b).bitstring_size(); i++){
	if(biparttable.BipartitionTable.at(b).get_bit(i)==false){
		zeros.push_back(i);
		}
	else{
		ones.push_back(i);
		}
	}
  if(biparttable.BipartitionTable.at(b).bitstring_size() > ::NUM_TAXA){ //make sure the next for loop terminates
	cerr << "generateQuartetsFromBipart error: Length of bitstring is greater than NUM_TAXA! This should never happen\n";
   	}
  else{
  	for(int i = biparttable.BipartitionTable.at(b).bitstring_size(); i < ::NUM_TAXA; i++){
		zeros.push_back(i);
		}
	}
//now we have ones and zeros, pass it to generateQuartetsFromOnesZeros
  return generateQuartetsFromOnesZerosSet(ones, zeros);
}

vector<quartet> generateQuartetsFromBipart(int b){
  //b is the index of the bipartition in biparttable.bipartitions
  //first, get the bitstring and map all of the 1s to one vector and all of the 
	//0s to another. NOTE with hetero trees this can be tricky, since taxa might
	//not exist even if the bitstring has a '0' for them
  vector<int> ones;
  vector<int> zeros;
  for(unsigned int i = 0; i < biparttable.bitstring_size(b); i++){
	if(biparttable.get_bit(b,i)==false){
		zeros.push_back(i);
		}
	else{
		ones.push_back(i);
		}
	}
  if(biparttable.bitstring_size(b) > ::NUM_TAXA){ //make sure the next for loop terminates
	cerr << "generateQuartetsFromBipart error: Length of bitstring is greater than NUM_TAXA! This should never happen\n";
   	}
  else{
  	for(unsigned int i = biparttable.bitstring_size(b); i < ::NUM_TAXA; i++){
		zeros.push_back(i);
		}
	}
//now we have ones and zeros, pass it to generateQuartetsFromOnesZeros
  return generateQuartetsFromOnesZeros(ones, zeros);
}


set<quartet> generateQuartetsFromOnesZerosSet(vector<int> ones, vector<int> zeros){
        set<quartet> ret;
	if(ones.size() < 2 || zeros.size() < 2)
	{
	return ret;
	}
	//first, store all of the pairs from zeros in a vector so we don't have too many
	//	nested loops
	  vector<iPair> zerosPairs;
	  for(int i = 0; i < zeros.size()-1; i++){
		for(int j = i+1; j < zeros.size(); j++){	
			zerosPairs.push_back(make_pair(zeros.at(i), zeros.at(j)));
			}
	 	}

	//now we have to choose 2 from each of these vectors. We know both are of at least size 2
	  for(int i = 0; i < ones.size()-1; i++){
		for(int j = i+1; j < ones.size(); j++){
			iPair p = make_pair(ones.at(i), ones.at(j));
			for(vector<iPair>::iterator k = zerosPairs.begin(); k!= zerosPairs.end(); k++){
				//add a quartet of both pairs
				ret.insert(quartet(p, *k));
				}
			}
		}
	return ret;
}

void insertQuartetsFromBipart(int b, set<quartet> &setty){
  vector<int> ones;
  vector<int> zeros;
  for(int i = 0; i < biparttable.BipartitionTable.at(b).bitstring_size(); i++){
	if(biparttable.BipartitionTable.at(b).get_bit(i)==false){
		zeros.push_back(i);
		}
	else{
		ones.push_back(i);
		}
	}
  if(biparttable.BipartitionTable.at(b).bitstring_size() > ::NUM_TAXA){ //make sure the next for loop terminates
	cerr << "generateQuartetsFromBipart error: Length of bitstring is greater than NUM_TAXA! This should never happen\n";
   	}
  else{
  	for(int i = biparttable.BipartitionTable.at(b).bitstring_size(); i < ::NUM_TAXA; i++){
		zeros.push_back(i);
		}
	}
	if(ones.size() < 2 || zeros.size() < 2)
	{
	return;
	}

	//first, store all of the pairs from zeros in a vector so we don't have too many
	//	nested loops
	  vector<iPair> zerosPairs;
	  for(int i = 0; i < zeros.size()-1; i++){
		for(int j = i+1; j < zeros.size(); j++){	
			zerosPairs.push_back(make_pair(zeros.at(i), zeros.at(j)));
			}
	 	}

	//now we have to choose 2 from each of these vectors. We know both are of at least size 2
	  for(int i = 0; i < ones.size()-1; i++){
		for(int j = i+1; j < ones.size(); j++){
			iPair p = make_pair(ones.at(i), ones.at(j));
			for(vector<iPair>::iterator k = zerosPairs.begin(); k!= zerosPairs.end(); k++){
				setty.insert(quartet(p, *k));
				}
			}
		}
}






vector<quartet> generateQuartetsFromOnesZeros(vector<int> ones, vector<int> zeros){
  vector<quartet> ret;
	if(ones.size() < 2 || zeros.size() < 2)
	{
	return ret;
	}

	//first, store all of the pairs from zeros in a vector so we don't have too many
	//	nested loops
	  vector<iPair> zerosPairs;
	  for(unsigned int i = 0; i < zeros.size()-1; i++){
		for(unsigned int j = i+1; j < zeros.size(); j++){	
			zerosPairs.push_back(make_pair(zeros.at(i), zeros.at(j)));
			}
	 	}

	//now we have to choose 2 from each of these vectors. We know both are of at least size 2
	  for(int i = 0; i < ones.size()-1; i++){
		for(int j = i+1; j < ones.size(); j++){
			iPair p = make_pair(ones.at(i), ones.at(j));
			for(vector<iPair>::iterator k = zerosPairs.begin(); k!= zerosPairs.end(); k++){
				//add a quartet of both pairs
				ret.push_back(quartet(p, *k));
				}
			}
		}
	return ret;
}

unsigned int getNumDifferentQuartets(int a, int b){
//quick optimization- if both ints are the same, we're comparing the same vector
 if(a==b){
  return 0;
	}
 boost::dynamic_bitset<> A = biparttable.non_trunc_bitstring(a);
 boost::dynamic_bitset<> B = biparttable.non_trunc_bitstring(b);
  return getNumDifferentQuartets(A, B);
}

unsigned int getNumSameQuartets(int a, int b){
//takes two bipartitions as input, returns the number of quartets that they do NOT share
 boost::dynamic_bitset<> A = ::biparttable.non_trunc_bitstring(a);
 boost::dynamic_bitset<> B = ::biparttable.non_trunc_bitstring(b);
 return getNumSameQuartets(A, B);
}


unsigned int getNumSameQuartets(boost::dynamic_bitset<> a, boost::dynamic_bitset<> b){
//takes two bipartitions as input, returns the number of quartets that they share
  //first, take the AND and the OR of the vectors. AND tells us which taxa are 1s in both
  //OR tells us which taxa are 0s in both

  //then, count the number of quartets that as if the bipartition was only the common 1s and 0s
  //this is the number of quartets both biparts have in common

  bPair biparts = bPair(a, b); 
  return getNumQuartets(biparts.sameOnes.size(), biparts.sameZeros.size());


//TODO- check if inverting the 1s and 0s produces a different result
}

unsigned int getNumDifferentQuartets(boost::dynamic_bitset<> a, boost::dynamic_bitset<> b){
//takes two bipartitions as input, returns the number of quartets that they do NOT share
  //first, take the AND and the OR of the vectors. AND tells us which taxa are 1s in both
  //OR tells us which taxa are 0s in both

  //then, count the number of quartets that as if the bipartition was only the common 1s and 0s
  //this is the number of quartets both biparts have in common
  bPair biparts = bPair(a, b);
  int commonOnes = biparts.sameOnes.size();
//  int commonZeros = chNumberOfZeros(chOR(a, b));  
  int commonZeros = biparts.sameZeros.size(); 
  unsigned int commonQuartets = getNumQuartets(commonOnes, commonZeros);
  unsigned int quartetsA = getNumQuartets(a);
  unsigned int quartetsB = getNumQuartets(b);

 return (quartetsA + quartetsB - (2*commonQuartets));

//TODO- check if inverting the 1s and 0s produces a different result
}

unsigned int getNumQuartets(int b){
  //takes the index of a bipartition and counts the number of quartets implied by it
  int nOnes = biparttable.number_of_ones(b);
  int nZeros = biparttable.number_of_zeros(b) + ::NUM_TAXA - biparttable.bitstring_size(b);
  //quartets = (nOnes choose 2) * (nZeros choose 2) 
//cout << "num ones is: " << nOnes << ", number of 0s is " << nZeros << endl;
  return(getNumQuartets(nOnes,nZeros));

}

unsigned int getNumQuartets(vector<bool> bipart){
  int nOnes, nZeros;
  nOnes = nZeros = 0;
  for(vector<bool>::iterator it = bipart.begin(); it!=bipart.end(); it++){
	if(*it){
		nOnes++;
		}
	else{
		nZeros++;
		}
	}
  return(getNumQuartets(nOnes,nZeros));
}

unsigned int getNumQuartets(boost::dynamic_bitset<> bipart){
  int nOnes, nZeros;
  nOnes = nZeros = 0;
  
  for(unsigned int i = 0; i< bipart.size(); i++){
	if(bipart[i]==1){
		nOnes++;
		}
	else{
		nZeros++;
		}
	}
 // cout << "get num quartets: num 0s is " << nZeros << " num 1s is: " << nOnes << endl;
  return(getNumQuartets(nOnes,nZeros));
}



unsigned int getNumQuartets(int nOnes, int nZeros){
  //NOTE- nOnes choose 2 simplifies to (nOnes * nOnes-1) / 2
  return (((nOnes * (nOnes-1)) / 2) * ((nZeros * (nZeros-1)) / 2));
  //return (  (factorial(nOnes) / (2 * factorial(nOnes-2))) * (factorial(nZeros) / (2 * factorial(nZeros-2))) );
}

set<quartet> generateQuartetsFromTree(int t){ //takes an int index of the tree from biparttable
  vector<unsigned int> bipartsInTree = biparts_in_tree(t); 
  set<quartet> retSet;
  for(vector<unsigned int>::iterator it = bipartsInTree.begin(); it!=bipartsInTree.end(); it++){
   	insertQuartetsFromBipart(*it, retSet);
	}
return retSet;
}

//NOTE- for one directional set different, use set_difference
set<quartet> generateDifferentQuartetsFromTrees(int a, int b){
 set<quartet> qa = generateQuartetsFromTree(a);
 set<quartet> qb = generateQuartetsFromTree(b);
 set<quartet> result;
 set_symmetric_difference(qa.begin(), qa.end(), qb.begin(), qb.end(), inserter(result, result.end()));
 return result;
}


set<quartet> generateSameQuartetsFromTrees(int a, int b){
 set<quartet> qa = generateQuartetsFromTree(a);
 set<quartet> qb = generateQuartetsFromTree(b);
 set<quartet> result;
 set_intersection(qa.begin(), qa.end(), qb.begin(), qb.end(), inserter(result, result.end()));
 return result;
}


void printSet(set<quartet> s){
	for(set<quartet>::iterator it = s.begin(); it!=s.end(); it++){
		cout << distance(s.begin(), it)+1 << ". ";		
		(*it).print();
		}
	}

void testOperatorsForQuartets(){
/*
  iPair a(4,8);
  iPair b(4,7);
  iPair d(3,20);
  iPair e(4,9);

  quartet C(b, a);
  quartet D(b, e);
 
  quartet A(a,b);
  quartet B(a,d);
  cout << "is A less than B? " << (A < B) << endl;
  cout << "is C less than D? " << (C < D) << endl;
*/
 /* set<quartet> Qset;
  set<quartet> result;
  Qset.insert(A);
  Qset.insert(B);
  set<quartet> Qset2; 
  Qset2.insert(C);
  Qset2.insert(D);
  Qset.insert(Qset2.begin(), Qset2.end());
 // merge(Qset.begin(), Qset.end(), Qset2.begin(), Qset2.end(), result.begin());
*/ 
// set<quartet> Qset = generateQuartetsFromTree(0);
 //printSet(Qset);

}

void bipartAnalysis(){ //dumps info about quartets in a bipartition
//the goal is to generate the number of quartets implied by each bipartition, and the similarity/difference matrix with each other
cout << "DIFFERENT/SIMILAR\n";
cout << "# qts:";
for(int i = 0; i < biparttable.BipartitionTable.size(); i++)
{
if(trivial_bipartitions.find(i)==trivial_bipartitions.end()) {cout << setw(8) << i;}

}

cout << endl;

for(int i = -1; ++i < biparttable.BipartitionTable.size() && trivial_bipartitions.find(i)==trivial_bipartitions.end();){
	cout << setw(6) << getNumQuartets(i);
	for(int j = 0; j < biparttable.BipartitionTable.size(); j++){	  
		string outty = "";
		stringstream s;
		s << getNumDifferentQuartets(i, j) << "/" << getNumSameQuartets(i,j);
		outty = s.str();
		cout << setw(8) << outty;

		}
cout << endl;
	}


}


void quartetAnalysis(int tree1, int tree2){ //info dump of quartet stuff from both trees
 set<quartet> q;
 q = generateQuartetsFromTree(tree1);
 cout << "Quartets from Tree : " << tree1 << endl;
 printSet(q);

 q = generateQuartetsFromTree(tree2);
 cout << "Quartets from Tree : " << tree2 << endl;
 printSet(q);
 q = generateSameQuartetsFromTrees(tree1, tree2);
 cout << "Same quartets across trees 0 and 1: \n";
 printSet(q);
  q = generateDifferentQuartetsFromTrees(tree1, tree2);
 cout << "different quartets across trees 0 and 1: \n";
 printSet(q);
 
 cout << "biparts in tree " << tree1 << ":\n";
 printVector(biparts_in_tree(tree1));

 cout << "biparts in tree " << tree2 << ":\n";
 printVector(biparts_in_tree(tree2));

}

//calculate trivial bipartitions

void calculateTrivialBipartitions(){
	for(int i = 0; i < biparttable.BipartitionTable.size(); i++){
	  if(biparttable.number_of_ones(i)<2 || ((biparttable.number_of_zeros(i)+::NUM_TAXA-biparttable.BipartitionTable.at(i).bitstring_size()) < 2)){
		trivial_bipartitions.insert(i);
		}
	}
}

void TESTSTUFF(){

/*
 calculateTrivialBipartitions();
// bipartAnalysis();
  quartetAnalysis(2, 4);

  testOperatorsForQuartets();  

  boost::dynamic_bitset<> a(4);
  boost::dynamic_bitset<> b(4);
  a.push_back(0);   a.push_back(0);   a.push_back(1);   a.push_back(1); 
  b.push_back(1);   b.push_back(1);   b.push_back(0);   b.push_back(0); 

  bPair tester = bPair(a,b);
  tester.printMetaData();
  


cout << "QUARTETS FROM BIPARTITION 7\n";
printQuartets(generateQuartetsFromBipart(7));
cout << "QUARTETS FROM BIPARTITION 12\n";
printQuartets(generateQuartetsFromBipart(12));
cout << "Number of quartets in bipartition 7: " << getNumQuartets(7) << endl;
cout << "Number of quartets in bipartition 12: " << getNumQuartets(12) << endl;


boost::dynamic_bitset<> bipart7(43ul);

//char bipart7[] = {0,0,1,1,0,1,0,1,1,1};
//boost::dynamic_bitset<> bipartFive;
//bipartFive.assign(bipart7, bipart7+10);
cout << "Number of quartets in bipartition 7: " << getNumQuartets(bipart7) << endl;

int diffQuartets = getNumDifferentQuartets(7,7);
cout << "different quartets between bipart 7 and 7 is: " << diffQuartets << endl;

diffQuartets = getNumDifferentQuartets(12,12);
cout << "different quartets between bipart 12 and 12 is: " << diffQuartets << endl;

diffQuartets = getNumDifferentQuartets(7,12);
cout << "different quartets between bipart 7 and 12 is: " << diffQuartets << endl;

cout << "about to get same quartets between 7 and 12" << endl;
printQuartets(generateSameQuartets(7, 12));
cout << "about to get different quartets between 7 and 12" << endl;
printQuartets(generateDifferentQuartets(7, 12));
*/

}





/* ALGORITHM DESCRIPTION FOR generateDifferentQuartets
  first, we need to know a few vectors:
	Same1- the taxa which are one across both biparts
	Same0- the taxa which are zero across both biparts	
	unique1A- taxa which are 1 only in bipart A
	unique1B- taxa which are 1 only in bipart B
	unique0A and unique0B

Algorithm description:

we want to generate all quartets from one bipartition but not the other

This means we can match any pair except if the both ones are generated by same1 and the 0s are generated by same0

Step 1- match all 0s with 1s not generated by same1s
first, take all permutations of taxa in unique1A, and also match each taxa in unique1A with a single taxon in same1. 
	This ensures that all of these pairs of 1s are those that can be generated by bipart1 but NOT bipart2
	Now, match all those pairs with the pairs that can be generated by 0s in bipart 1

For instance, suppose bipart A has 1s for taxa {0, 2, 3, 5, and 7} and both biparts had 1 at {5 7}. This means the
	unique pairs of 1s are {0/2, 0/3, 2/3,   |||  0/5, 2/5, 3/5, 0/7, 2/7, 3/7}. The ||| represents
	where the algorithm goes from making pairs from uniqueA to matching 1 taxon from uniqueA with each taxon 	 in shared1. 

Step 2- match all 1s generated with same1s with 0s unique to bipartition (i.e. not in same 0)
	Take same1s and create pairs out of it. Then generate pairs of 0s not possible by same 0 by 
	taking all permutations of unique0, and also matching each taxon in same0 with each taxon in unique0 
Step 3- repeat steps 1 and 2 with second bipartition


*/

