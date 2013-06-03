#include "InterfaceFunctions.h"

int HashCS(vector<int> input_from_int) {
	 vector<string> temp;
	 ofstream newick_file;
	 newick_file.open("temp/newick.txt");
	 temp = to_newick(input_from_int);
	 for(unsigned int i=0; i<temp.size(); i++)    
	      newick_file << temp[i] << endl;
	 newick_file.close();
	   
	 std::stringstream consensus_string;
	 consensus_string << "cd ..; cd HashCS/; ./hashcs ../TreeHouse/temp/newick.txt -o ../TreeHouse/temp/consensus.txt";
	 system(consensus_string.str().c_str());
	 
	 cout << "***Results will be saved to consensus.txt in temp directory***" << endl;
	 //return consensus_string.str();
	 return 0;
 }
 
 int HashCS(int x) {
	 ofstream newick_file;
	 newick_file.open("temp/newick.txt");
	 newick_file << to_newick(x);
	 newick_file.close();
 
	 std::stringstream consensus_string;
	 consensus_string << "cd ..; cd HashCS/; ./hashcs ../TreeHouse/temp/newick.txt -o ../TreeHouse/temp/consensus.txt";
	 system(consensus_string.str().c_str());
	 
	 cout << "***Results will be saved to consensus.txt in temp directory***" << endl;
	 //return consensus_string.str();
	 return 0;
 }
 

int phlash(vector<int> input_from_int) {
	vector<string> temp;
	ofstream newick_file;
	newick_file.open("temp/newick.txt");
	temp = to_newick(input_from_int);
	for(unsigned int i=0; i<temp.size(); i++)    
	     newick_file << temp[i] << endl;
	newick_file.close();
	  
	std::stringstream phlash_string;
	cout << temp.size() << endl;
	phlash_string << "cd ..; cd Phlash/; ./phlash ../TreeHouse/temp/newick.txt " << temp.size() << " -d 0 -m -o ../TreeHouse/temp/phlash.out";
	system(phlash_string.str().c_str());
	
	cout << "***Results will be saved to phlash.out in temp directory***" << endl;
	
	return 0;
}

int phlash(set<unsigned int> input_from_int){
	vector<string> temp;
	ofstream newick_file;
	newick_file.open("temp/newick.txt");
	temp = to_newick(input_from_int);
	for(unsigned int i=0; i<temp.size(); i++)    
	     newick_file << temp[i] << endl;
	newick_file.close();
	  
	std::stringstream phlash_string;
	cout << temp.size() << endl;
	phlash_string << "cd ..; cd Phlash/; ./phlash ../TreeHouse/temp/newick.txt " << temp.size() << " -d 0 -m -o ../TreeHouse/temp/phlash.out";
	system(phlash_string.str().c_str());
	
	cout << "***Results will be saved to phlash.out in temp directory***" << endl;
	
	return 0;
}

int phlash(int x) {
	ofstream newick_file;
	newick_file.open("temp/newick.txt");
	newick_file << to_newick(x);
	newick_file.close();

	std::stringstream phlash_string;
	phlash_string << "cd ..; cd Phlash/; ./phlash ../TreeHouse/temp/newick.txt " << x << " -d 0 -m -o ../TreeHouse/temp/phlash.out";
	system(phlash_string.str().c_str());
		
	cout << "***Results will be saved to phlash.out in temp directory***" << endl;

	return 0;
}


int quick_quartet(vector<int> input_from_int) {
	vector<string> temp;
	ofstream newick_file;
	newick_file.open("temp/newick.txt");
	temp = to_newick(input_from_int);
	for(unsigned int i=0; i<temp.size(); i++)    
	     newick_file << temp[i] << endl;
	newick_file.close();
	  
	std::stringstream qq_string;
	qq_string << "cd ..; cd QuickQuartet/; ./qq ../TreeHouse/temp/newick.txt -b | tee ../TreeHouse/temp/qq.txt";
	system(qq_string.str().c_str());
	
	cout << "***Results will be saved to qq.txt in temp directory***" << endl;
	
	return 0;
}

int quick_quartet(set<unsigned int> input_from_int) {
	vector<string> temp;
	ofstream newick_file;
	newick_file.open("temp/newick.txt");
	temp = to_newick(input_from_int);
	for(unsigned int i=0; i<temp.size(); i++)    
	     newick_file << temp[i] << endl;
	newick_file.close();
	  
	std::stringstream qq_string;
	qq_string << "cd ..; cd QuickQuartet/; ./qq ../TreeHouse/temp/newick.txt -b | tee ../TreeHouse/temp/qq.txt";
	system(qq_string.str().c_str());
	
	cout << "***Results will be saved to qq.txt in temp directory***" << endl;
	
	return 0;
}

int quick_quartet(int x) {
	ofstream newick_file;
	newick_file.open("temp/newick.txt");
	newick_file << to_newick(x);
	newick_file.close();

	std::stringstream qq_string;
	qq_string << "cd ..; cd QuickQuartet/; ./qq ../TreeHouse/temp/newick.txt -b | tee ../TreeHouse/temp/qq.txt";
	system(qq_string.str().c_str());
	
	cout << "***Results will be saved to qq.txt in temp directory***" << endl;

	return 0;
}
