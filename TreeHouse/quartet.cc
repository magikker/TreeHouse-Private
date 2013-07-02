#include "quartet.h"

using namespace std;



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

void getOnesZerosFromBipart(int bipart, vector<unsigned int> &ones, vector<unsigned int> &zeros){
  boost::dynamic_bitset<> b = biparttable.BipartTable.at(bipart).get_bitstring();



}

char isQuartetImplied(quartet x, unsigned int bipart){
//note- this can return three values
//0- quartet is DISPROVEN
//1- quartet is PROVEN
//2- we don't know (e.g. we're looking for a quartet and we find three 0s and a 1
  boost::dynamic_bitset<> b = biparttable.BipartTable.at(bipart).get_bitstring();
  iPair A = x.getA();
  iPair B = x.getB();
  int numOnes = (int)b[A.first] + (int)b[A.second] + (int)b[B.first] + (int)b[B.second];
  if(numOnes!=2){
	//in the significant places we're looking at, there is not a 2 and 2 pairing
	//for instance, in the bipartition 11000, we can't check a quartet over BCDE because there are 3 0s and a 1
	return 2;
	}
  else{
	bool left = b[A.first];
	return (b[A.second]==left && b[B.first]!=left && b[B.second]!=left) ? 1 : 0;
		
	}
}

char isQuartetImplied(quartet x, set<unsigned int> biparts){
  for(set<unsigned int>::iterator it = biparts.begin(); it!=biparts.end(); it++){ //for each bipartition
	//cout << "bipart is: "; biparttable.BipartitionTable.at(*it).print_bitstring(true);
	//cout << "quartet is: "; x.print();
	//cout << "isQuartetImplied? " << (int)isQuartetImplied(x,*it) << endl;
	if(isQuartetImplied(x, *it)==1){
		return 1;
		}
	else if(isQuartetImplied(x, *it)==0){
		return 0;
		}
	}
  
  return 2;
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

set<quartet> generateConflictingQuartets(int bipart1, int bipart2){
  //right now, we will brute force this because we need to results for testing
  set<quartet> retSet;
  set<unsigned int> ones, zeros;
  //ones = biparttable.BipartitionTable.at(bipart1).getOnes();
  //zeros = biparttable.BipartitionTable.at(bipart1).getZeros();
  set<quartet> bipartsA = generateQuartetsFromBipartSet(bipart1);
  for(set<quartet>::iterator it = bipartsA.begin(); it!=bipartsA.end(); it++){
	if(isQuartetImplied(*it, bipart2)==0){
		retSet.insert(*it);
		}
	}
  return retSet;
}

set<quartet> generateConflictingQuartetsGroup(int bipart1, set<unsigned int> biparts){
	//right now, we will brute force this because we need to results for testing
	set<quartet> retSet;
	set<unsigned int> ones, zeros;
	//ones = biparttable.BipartitionTable.at(bipart1).getOnes();
	//zeros = biparttable.BipartitionTable.at(bipart1).getZeros();
	set<quartet> bipartsA = generateQuartetsFromBipartSet(bipart1);
	for(set<quartet>::iterator it = bipartsA.begin(); it!=bipartsA.end(); it++){
	if(isQuartetImplied(*it, biparts)==0){
		retSet.insert(*it);
		}
	}
	return retSet;
}

set<quartet> generateConflictingQuartets2(int bipart1, int bipart2){
	set<quartet> retSet;
	boost::dynamic_bitset<> a, b;
	a = biparttable.BipartTable.at(bipart1).get_bitstring();
	b = biparttable.BipartTable.at(bipart2).get_bitstring();
	a.resize(::biparttable.lm.size()); //restore trailing 0s in bipartitions
	b.resize(::biparttable.lm.size());
	vector<unsigned int> sharedOnes, sharedZeros, unsharedOnes, unsharedZeros;
	for(int i = 0; i < a.size(); i++){
	if(a[i]==1){
		if(b[i]==1){
			sharedOnes.push_back(i);
			}
		else{
			unsharedOnes.push_back(i);
			}
		}
	else{
		if(b[i]==0){
			sharedZeros.push_back(i);
			}
		else{
			unsharedZeros.push_back(i);
			}
		}
	}
	//now we have all the sets of ones and zeros, create quartets by taking
	//  shared1s/unshared1s |  shared0s/unshared0s
	for(vector<unsigned int>::iterator a = sharedOnes.begin(); a!=sharedOnes.end(); a++){
	for(vector<unsigned int>::iterator b = unsharedOnes.begin(); b!=unsharedOnes.end(); b++){
		for(vector<unsigned int>::iterator c = sharedZeros.begin(); c!=sharedZeros.end(); c++){
			for(vector<unsigned int>::iterator d = unsharedZeros.begin(); d!=unsharedZeros.end(); d++){
				retSet.insert(quartet(*a,*b,*c,*d));
				}
			}
		}
	} 
	return retSet;
}

void printQPair(qPair a){
	  set<iPair> one = a.first;
	  set<iPair> two = a.second;
	  if(one.size()==0 || two.size()==0){
		cout << "empty!\n";
		}
	  else{
	  cout << "{ ";
	  for(set<iPair>::iterator i = one.begin(); i!=one.end(); i++){ //for each iPair in first
		cout << (*i).first << "," << (*i).second << "  ";
		}
	  cout << "}\n{ ";
	  for(set<iPair>::iterator j = two.begin(); j!=two.end(); j++){ //for each iPair in first
		cout << (*j).first << "," << (*j).second << "  ";
		}
	 cout << "}\n";
	}
}

qPair mergeQPair(qPair a, qPair b){

         set<iPair> one, two;
	 merge(a.first.begin(), a.first.end(), b.first.begin(), b.first.end(), inserter(one, one.begin()));
	 merge(a.second.begin(), a.second.end(), b.second.begin(), b.second.end(), inserter(two, two.begin()));
	 return make_pair(one, two);
}

bool isQPairEmpty(qPair a){
  return (a.first.size()==0 || a.second.size()==0);
}


qPair generateConflictingQuartets3(int bipart1, int bipart2){
	  qPair retPair;
	  boost::dynamic_bitset<> a, b;
	  a = biparttable.BipartTable.at(bipart1).get_bitstring();
	  b = biparttable.BipartTable.at(bipart2).get_bitstring();
	  a.resize(::biparttable.lm.size()); //restore trailing 0s in bipartitions
	  b.resize(::biparttable.lm.size());
	  vector<unsigned int> sharedOnes, sharedZeros, unsharedOnes, unsharedZeros;
	  for(int i = 0; i < a.size(); i++){
		if(a[i]==1){
			if(b[i]==1){
				sharedOnes.push_back(i);
				}
			else{
				unsharedOnes.push_back(i);
				}
			}
		else{
			if(b[i]==0){
				sharedZeros.push_back(i);
				}
			else{
				unsharedZeros.push_back(i);
				}
			}
		}
	  for(vector<unsigned int>::iterator a = sharedOnes.begin(); a!=sharedOnes.end(); a++){
		for(vector<unsigned int>::iterator b = unsharedOnes.begin(); b!=unsharedOnes.end(); b++){
				retPair.first.insert(make_pair(*a, *b));
			}
		} 

	  for(vector<unsigned int>::iterator c = sharedZeros.begin(); c!=sharedZeros.end(); c++){
		for(vector<unsigned int>::iterator d = unsharedZeros.begin(); d!=unsharedZeros.end(); d++){
				retPair.second.insert(make_pair(*c, *d));
				}
		}

	  return retPair;
}

qPair generateConflictingQuartetsGroup2(int bipart1, set<unsigned int> biparts){
  qPair retPair;
  for(set<unsigned int>::iterator it = biparts.begin(); it!=biparts.end(); it++){
	retPair = mergeQPair(retPair, generateConflictingQuartets3(bipart1, *it));
	}
  return retPair;
}

unsigned int numConflictingQuartets(int bipart1, int bipart2){
  boost::dynamic_bitset<> a, b;
  a = biparttable.BipartTable.at(bipart1).get_bitstring();
  b = biparttable.BipartTable.at(bipart2).get_bitstring();
  a.resize(::biparttable.lm.size()); //restore trailing 0s in bipartitions
  b.resize(::biparttable.lm.size());
  unsigned int sharedOnes, sharedZeros, unsharedOnes, unsharedZeros;
  sharedOnes = sharedZeros = unsharedOnes = unsharedZeros = 0;
  for(int i = 0; i < a.size(); i++){
	if(a[i]==1){
		if(b[i]==1){
			sharedOnes++;
			}
		else{
			unsharedOnes++;
			}
		}
	else{
		if(b[i]==0){
			sharedZeros++;
			}
		else{
			unsharedZeros++;
			}
		}
	}

  return (sharedOnes*unsharedOnes) * (sharedZeros*unsharedZeros);
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
  for(int i = 0; i < biparttable.BipartTable.at(b).bitstring_size(); i++){
	if(biparttable.BipartTable.at(b).get_bit(i)==false){
		zeros.push_back(i);
		}
	else{
		ones.push_back(i);
		}
	}
  if(biparttable.BipartTable.at(b).bitstring_size() > ::biparttable.lm.size()){ //make sure the next for loop terminates
	cerr << "generateQuartetsFromBipart error: Length of bitstring is greater than NUM_TAXA! This should never happen\n";
   	}
  else{
  	for(int i = biparttable.BipartTable.at(b).bitstring_size(); i < ::biparttable.lm.size(); i++){
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
  if(biparttable.bitstring_size(b) > ::biparttable.lm.size()){ //make sure the next for loop terminates
	cerr << "generateQuartetsFromBipart error: Length of bitstring is greater than NUM_TAXA! This should never happen\n";
   	}
  else{
  	for(unsigned int i = biparttable.bitstring_size(b); i < ::biparttable.lm.size(); i++){
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
	  for(unsigned int i = 0; i < zeros.size()-1; i++){
		for(unsigned int j = i+1; j < zeros.size(); j++){	
			zerosPairs.push_back(make_pair(zeros.at(i), zeros.at(j)));
			}
	 	}

	//now we have to choose 2 from each of these vectors. We know both are of at least size 2
	  for(unsigned int i = 0; i < ones.size()-1; i++){
		for(unsigned int j = i+1; j < ones.size(); j++){
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
  for(unsigned int i = 0; i < biparttable.BipartTable.at(b).bitstring_size(); i++){
	if(biparttable.BipartTable.at(b).get_bit(i)==false){
		zeros.push_back(i);
		}
	else{
		ones.push_back(i);
		}
	}
  if(biparttable.BipartTable.at(b).bitstring_size() > ::biparttable.lm.size()){ //make sure the next for loop terminates
	cerr << "generateQuartetsFromBipart error: Length of bitstring is greater than NUM_TAXA! This should never happen\n";
   	}
  else{
  	for(unsigned int i = biparttable.BipartTable.at(b).bitstring_size(); i < ::biparttable.lm.size(); i++){
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
	  for(unsigned int i = 0; i < zeros.size()-1; i++){
		for(unsigned int j = i+1; j < zeros.size(); j++){	
			zerosPairs.push_back(make_pair(zeros.at(i), zeros.at(j)));
			}
	 	}

	//now we have to choose 2 from each of these vectors. We know both are of at least size 2
	  for(unsigned int i = 0; i < ones.size()-1; i++){
		for(unsigned int j = i+1; j < ones.size(); j++){
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
	  for(unsigned int i = 0; i < ones.size()-1; i++){
		for(unsigned int j = i+1; j < ones.size(); j++){
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
  int nZeros = biparttable.number_of_zeros(b) + ::biparttable.lm.size() - biparttable.bitstring_size(b);
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

set<quartet> generateQuartetsFromBiparts(set<unsigned int> s){
  set<quartet> retSet;
  for(set<unsigned int>::iterator it = s.begin(); it!=s.end(); it++){
   	insertQuartetsFromBipart(*it, retSet);
	}
  return retSet;

}


set<quartet> generateDifferentQuartetsFromTrees(int a, int b){
 //COMPARISON OF ALL-TREE1 TO ALL-TREE2
 set<quartet> qa = generateQuartetsFromTree(a);
 set<quartet> qb = generateQuartetsFromTree(b);
 set<quartet> result;
 set_difference(qa.begin(), qa.end(), qb.begin(), qb.end(), inserter(result, result.end()));
 return result;
}

set<quartet> generateDifferentQuartetsFromTrees2(int a, int b){
//first, take all bipartitions from each tree and add them into sets
 set<unsigned int> t1Biparts(::biparttable.inverted_index.at(a).begin(), ::biparttable.inverted_index.at(a).end());
 set<unsigned int> t2Biparts(::biparttable.inverted_index.at(b).begin(), ::biparttable.inverted_index.at(b).end());
 set<unsigned int> unique;
 set_difference(t1Biparts.begin(), t1Biparts.end(), t2Biparts.begin(), t2Biparts.end(), inserter(unique, unique.end()));
//now we have unique, a set of all bipartitions that are unique to tree 1 (note that this set_difference ISNT symmetric)
 set<quartet> retSet;
 for(set<unsigned int>::iterator it = unique.begin(); it!=unique.end(); it++){ //for each bipartition unique to tree1
	set<quartet> qs = generateQuartetsFromBipartSet(*it);
	for(set<quartet>::iterator j = qs.begin(); j != qs.end(); j++){
		if(!isQuartetImplied(*j, t2Biparts)){
			retSet.insert(*j);
			}
		}
	} 
  return retSet;
}

set<quartet> generateDifferentQuartetsFromTrees3(int a, int b){
 //COMPARISON OF ALL UNIQUE-TREE1 TO ALL OF TREE2
 set<unsigned int> t1Biparts(::biparttable.inverted_index.at(a).begin(), ::biparttable.inverted_index.at(a).end());
 set<unsigned int> t2Biparts(::biparttable.inverted_index.at(b).begin(), ::biparttable.inverted_index.at(b).end());
 set<unsigned int> unique;
 set_difference(t1Biparts.begin(), t1Biparts.end(), t2Biparts.begin(), t2Biparts.end(), inserter(unique, unique.end()));
 set<quartet> retSet, qUnique, qb;
 qUnique = generateQuartetsFromBiparts(unique);
 qb = generateQuartetsFromTree(b);
 set_difference(qUnique.begin(), qUnique.end(), qb.begin(), qb.end(), inserter(retSet, retSet.end()));
  
  return retSet;
}

set<quartet> generateDifferentQuartetsFromTrees4(int a, int b){
 //COMPARISON OF UNIQUE-TREE1 TO UNIQUE-TREE2
 set<unsigned int> t1Biparts(::biparttable.inverted_index.at(a).begin(), ::biparttable.inverted_index.at(a).end());
 set<unsigned int> t2Biparts(::biparttable.inverted_index.at(b).begin(), ::biparttable.inverted_index.at(b).end());
 set<unsigned int> unique, unique2;
 set_difference(t1Biparts.begin(), t1Biparts.end(), t2Biparts.begin(), t2Biparts.end(), inserter(unique, unique.end()));
 set_difference(t2Biparts.begin(), t2Biparts.end(), t1Biparts.begin(), t1Biparts.end(), inserter(unique2, unique2.end()));
 set<quartet> retSet, qUnique, qUnique2;
 qUnique = generateQuartetsFromBiparts(unique);
 qUnique2 = generateQuartetsFromBiparts(unique2);
 set_difference(qUnique.begin(), qUnique.end(), qUnique2.begin(), qUnique2.end(), inserter(retSet, retSet.end()));
  
  return retSet;
}


set<quartet> generateSameQuartetsFromTrees(int a, int b){
 set<quartet> qa = generateQuartetsFromTree(a);
 set<quartet> qb = generateQuartetsFromTree(b);
 set<quartet> result;
 set_intersection(qa.begin(), qa.end(), qb.begin(), qb.end(), inserter(result, result.end()));
 return result;
}

unsigned int quartet_distance(int tree1, int tree2){
//first, get set of non-trivial bipartitions which are shared by the trees

  vector<unsigned int> biparts1 = ::biparttable.inverted_index.at(tree1);
  vector<unsigned int> biparts2 = ::biparttable.inverted_index.at(tree2);
  
 
 set<unsigned int> b1(biparts1.begin(), biparts1.end());
 set<unsigned int> b2(biparts2.begin(), biparts2.end()); 
 //copy(biparts1.begin(), biparts1.end(), inserter(b1, b1.end()));
  //copy(biparts2.begin(), biparts2.end(), inserter(b2, b2.end()));

  cout << "printing biparts 1";
  for(set<unsigned int>::iterator it = b1.begin(); it!=b1.end(); it++){
	cout << *it << "   ";
	}cout << endl;
  
  set<unsigned int> b1nonTrivial, b2nonTrivial, b1Unique, sharedBiparts; 
 
 set_difference(b1.begin(), b1.end(), biparttable.trivial_bipartitions.begin(), biparttable.trivial_bipartitions.end(), inserter(b1nonTrivial, b1nonTrivial.begin()));

  set_difference(b2.begin(), b2.end(), biparttable.trivial_bipartitions.begin(), biparttable.trivial_bipartitions.end(), inserter(b2nonTrivial, b2nonTrivial.begin()));

  set_difference(b1.begin(), b1.end(), b2.begin(), b2.end(), inserter(b1Unique, b1Unique.begin()));
  
  set_intersection(b1.begin(), b1.end(), b2.begin(), b2.end(), inserter(sharedBiparts, sharedBiparts.begin()));
  
  cout << "printing biparts 1 non trivial  ";
  for(set<unsigned int>::iterator it = b1nonTrivial.begin(); it != b1nonTrivial.end(); it++){
	cout << *it << "   ";
	}cout << endl;

  //now, compare all bipartitions in b1Unique to those in b1Unique and in set_difference
  unsigned int total = 0;
  for(unsigned int i = 0; i < b1Unique.size(); i++){
        
  	for(int k = 0; k < sharedBiparts.size(); k++){
		total=getNumDifferentQuartets(i,k);
		}  
	for(unsigned int j = i+1; j<b1Unique.size(); j++){
		total+=getNumDifferentQuartets(i,j);
		}
	}


  return total;
}


void shared_quartets_strict(set<unsigned int> trees){
//takes a set of trees, returns via print the quartets which are present in all of the trees
  set<quartet> strictly_shared;
  
  //first, generate all quartets from the first tree
  strictly_shared = generateQuartetsFromTree(*trees.begin());
  //now, go through all remaining trees and yank out values that aren't found
  for(set<unsigned int>::iterator it = trees.begin(); it!=trees.end(); it++){
	if(it!=trees.begin()){
	set<quartet> treeSet = generateQuartetsFromTree(*it);
	set<quartet> temp;
	set_intersection(strictly_shared.begin(), strictly_shared.end(), treeSet.begin(), treeSet.end(), inserter(temp, temp.end()));
	strictly_shared = temp;}
	}
  cout << "Now printing quartets which are common to all trees in the treeset:\n";
  printSet(strictly_shared);
}

void shared_quartets_majority(set<unsigned int> trees){
//takes a set of trees, returns via print the quartets which are present in all of the trees
  set<quartet> quartets;
  map<quartet, unsigned int> counter;
  
  for(set<unsigned int>::iterator it = trees.begin(); it!=trees.end(); it++){
	set<quartet> treeSet = generateQuartetsFromTree(*it);
	for(set<quartet>::iterator j = treeSet.begin(); j!= treeSet.end(); j++){
		quartet q = *j;
		if(counter.count(q)==0){ //we have a new quartet we haven't seen
			counter[q] = 1;
			}
		else{ //add to the count
			unsigned int count = counter[q];
			counter.erase(q);
			counter[q] = count+1;
			}
		}
	}
  int threshhold = trees.size()/2 + trees.size()%2;
  //copy quartets in counter to set quartets if they meet the threshhold
  for(map<quartet, unsigned int>::iterator it = counter.begin(); it!=counter.end(); it++){
	if(it->second >= threshhold){
		quartets.insert(it->first);
		}
	}
  cout << "Now printing quartets which are in at least 50% of all trees in the treeset:\n";
  printSet(quartets);


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

void ktetAnalysis(int t){
//gets all the bipartitions of a tree and gets the ones and zero sets
vector<unsigned int> biparts = ::biparttable.inverted_index.at(t);
set<set<unsigned int>> bSets; 
for(int i = 0; i < biparts.size(); i++){
	
	if(!::biparttable.is_trivial(biparts.at(i))){
		set<unsigned int> ones, zeros;
		ones = biparttable.getOnes(biparts.at(i));
		zeros = biparttable.getZeros(biparts.at(i));
		cout << "Bipartition " << biparts.at(i) << ":" << endl;	
		cout << "ones: "; printSetCompactTwo(ones);
		cout << "zeros: "; printSetCompactTwo(zeros);
		bSets.insert((ones.size() < zeros.size()) ? ones : zeros);
		}
	}
cout << "PRINTING SET OF SETS:\n";
for(int i = 2; i < 26; i++){
	for(set<set<unsigned int>>::iterator it = bSets.begin(); it!=bSets.end(); it++)
	{
		if((*it).size()==i) {printSetCompactTwo(*it); cout << endl;}
		}
	}
}

void TESTSTUFF(){
testGenerateConflictingQuartets();
//testNumConflictingQuartets();
//testGenerateDifferentQuartets();
//testGeenerateDifferentQuartetsFromTrees(){
}

void testGenerateConflictingQuartets(){
  
 
unsigned int numNonEmpty = 0;
unsigned int total = 0;

for(int i = 0; i < 100; i++){
	for(int k = i + 1; k < 100; k++){
		total++;
		if(!isQPairEmpty(generateConflictingQuartets3(i, k))){
			numNonEmpty++;
			//cout << i << "," << k << " is nonempty\n";
			}
		}
	}
	cout << "total: " << total << endl;
	cout << "nonEmpty: " << numNonEmpty << endl;
	cout << "ratio: " << (double)numNonEmpty/(double)total << endl;

/*
  vector<unsigned int> rf, rf2;
  //rf is the one sided difference from t1 to t2. rf2 is from t2 to t1
  vector<unsigned int> t1 = ::biparttable.inverted_index.at(258);
  vector<unsigned int> t2 = ::biparttable.inverted_index.at(561);
  set_difference(t1.begin(), t1.end(), t2.begin(), t2.end(), inserter(rf, rf.begin()));
  set_difference(t2.begin(), t2.end(), t1.begin(), t1.end(), inserter(rf2, rf2.begin()));
  set<unsigned int> group;
  for(int r = 0; r < rf2.size(); r++){
	group.insert(rf2.at(r));
	}

  for(int k = 0; k < rf.size(); k++){
	set<quartet> conflict = generateConflictingQuartetsGroup(k, group);
	if(conflict.size() > 0){
		for(int i = 0; i < rf2.size(); i++){
		cout << "conflicting quartets between bipartitions " << rf.at(k) << " and " << rf2.at(i) << " are:\n";
		printQPair(generateConflictingQuartets3(rf.at(k), rf2.at(i)));
		group.insert(rf2.at(i));
		cout << endl;
		}
		cout << "conflicting quartets between " << rf.at(k) << " and group:\n";
		printSet(generateConflictingQuartetsGroup(rf.at(k), group));
		//return;
		}
	else{
		cout << "bipartition " << rf.at(k) << " created no conflict!\n";
		}
	}  
*/
/*
  int a, g1, g2, g3, g4;
  a = 9;
  g1 =3; g2 = 10; g3 = 11; g4 = 12;
  set<unsigned int> group;
  group.insert(g1); group.insert(g2); group.insert(g3); group.insert(g4); 
  cout << "conflicting quartets between " << a << " and " << g1 << ":\n";
  printQPair(generateConflictingQuartets3(a, g1));

  cout << "conflicting quartets between " << a << " and " << g2 << ":\n";
  printQPair(generateConflictingQuartets3(a, g2));

  cout << "conflicting quartets between " << a << " and " << g3 << ":\n";
  printQPair(generateConflictingQuartets3(a, g3));

  cout << "conflicting quartets between " << a << " and " << g4 << ":\n";
  printQPair(generateConflictingQuartets3(a, g4));

  cout << "conflicting quartets between " << a << " and group:\n";
  printSet(generateConflictingQuartetsGroup(a, group));
  //qPair x = generateConflictingQuartets3(2,4);
  //printQPair(x);
  */


  //cout << "quartet distance is:\n";
 // printSet(generateDifferentQuartetsFromTrees(1, 8));

/*  set<quartet> conflictingQuartets;
  set<unsigned int> group;
  group.insert(4); group.insert(13);

  cout << "conflicting quartets between 2 and 4:\n";
  printSet(generateConflictingQuartets(2,4));
  cout << "conflicting quartets between 2 and 13:\n";
  printSet(generateConflictingQuartets(2,13));
  cout << "conflicting quartets between 2 and group:\n";
  printSet(generateConflictingQuartetsGroup(2, group));
  cout << "conflicting quartets between 6 and 4:\n";
  printSet(generateConflictingQuartets(6,4));
  cout << "conflicting quartets between 6 and 13:\n";
  printSet(generateConflictingQuartets(6,13));
  cout << "conflicting quartets between 6 and group:\n";
  printSet(generateConflictingQuartetsGroup(6, group));


  cout << "quartet distance is:\n";
  printSet(generateDifferentQuartetsFromTrees(0, 1));
*/
 /*
 vector<unsigned int> rf;
 for(int i = 0; i < 10; i++){
 	for(int j = 2000; j < 2050; j++){
	 vector<unsigned int> t1 = ::biparttable.inverted_index.at(i);
	 vector<unsigned int> t2 = ::biparttable.inverted_index.at(j);
	 set_difference(t1.begin(), t1.end(), t2.begin(), t2.end(), inserter(rf, rf.begin()));
	 cout << "rf distance between " << i << " and " << j << " is: " << rf.size() << endl;
	 rf.clear();
	}
}
*/



}


void testNumConflictingQuartets(){
	//TO TEST NUMCONFLICTINGQUARTETS ON LARGE DATA SETS
	unsigned long numBipartsLookedAt = 0;
	unsigned long nonZero = 0;
	if(::biparttable.lm.size()>500){
	for(int i = 0; i < 3000; i++){
		for(int j = 0; j < 3000; j++){
			//cout << "conflicting quartets between bipartition " << i << " and " << j << ":";
			numBipartsLookedAt++;		
			unsigned int x =  numConflictingQuartets(i,j);
			if(x!=0) {
				cout << x; 
				nonZero++;
					}
		numConflictingQuartets(i,j);
			}
		cout << "done with " << i << endl << endl;
		}
	}
	cout << "done! Total number of bipartition groups looked at is " << numBipartsLookedAt << endl;
	cout << "Of those, there were " << nonZero << " number of non-zero pairings\n";
	double ratio = (double)nonZero/numBipartsLookedAt;
	cout << "The proportion of non-zero pairings is " << ratio*100 << "%" << endl;

	//TO TEST numConflictingQuartets ON SMALL DATA SETS
	/*for(int i = 0; i < 10; i++){
		for(int j = 0; j < 10; j++){
			cout << "conflicting quartets between bipartition " << i << " and " << j << ":";
			cout << numConflictingQuartets(i,j) << endl;
			}
		cout << "done with " << i << endl << endl;	
		}*/

}

void testGenerateDifferentQuartets(){
	cout << "TREE 0:" << endl;
	ktetAnalysis(0);
	cout << "TREE 1:" << endl;
	ktetAnalysis(1);
	set<quartet> q2, q6, q4, q13, intersect1, intersect2, intersect3, conflict1, conflict2, conflict3, conflict4;
	cout << "Different quartets from 0 and 1:" << endl;
	printSet(generateDifferentQuartetsFromTrees3(0,1));

	cout << "quartets implied by bipartition 2:\n";
	printSet(q2 = generateQuartetsFromBipartSet(2));
	cout << "quartets implied by bipartition 6:\n";
	printSet(q6 = generateQuartetsFromBipartSet(6));
	cout << "quartets implied by bipartition 4:\n";
	printSet(q4 = generateQuartetsFromBipartSet(4));
	cout << "quartets implied by bipartition 13:\n";
	printSet(q13 = generateQuartetsFromBipartSet(13));



	cout << "conflicting quartets between biparts 2 and 4:\n";
	printSet(conflict1 = generateConflictingQuartets(2,4));
	cout << "conflicting quartets between biparts 2 and 13:\n";
	printSet(conflict2 = generateConflictingQuartets(2,13));
	cout << "conflicting quartets between biparts 6 and 4:\n";
	printSet(conflict3 = generateConflictingQuartets(6,4));
	cout << "conflicting quartets between biparts 6 and 13:\n";
	printSet(conflict4 = generateConflictingQuartets(6,13));
	


	


}

void testGeenerateDifferentQuartetsFromTrees(){
	
	cout << "generate different quartets from trees VERSION 0:\n";
	for(int i = 0; i < 9; i++){
		set<quartet> qs = generateDifferentQuartetsFromTrees(0,i);
		cout << "qs size for trees 0 and " << i << " is: " << qs.size() << endl << endl;

	}

	cout << "generate different quartets from trees VERSION 1:\n";
	for(int i = 0; i < 9; i++){
		set<quartet> qs = generateDifferentQuartetsFromTrees3(0,i);
		cout << "qs size for trees 0 and " << i << " is: " << qs.size() << endl << endl;

	}

	cout << "generate different quartets from trees VERSION 2:\n";
	for(int i = 0; i < 9; i++){
		set<quartet> qs = generateDifferentQuartetsFromTrees4(0,i);
		cout << "qs size for trees 0 and " << i << " is: " << qs.size() << endl << endl;

	}
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

