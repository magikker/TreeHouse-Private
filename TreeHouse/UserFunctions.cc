#include "UserFunctions.h"

//Return the type of each arg. 
string get_arg_types(vector<pqlsymbol * > arglist) {
  string argtypes = "(";
  string arg;
  int otype;
  int dtype;
  for (unsigned int x = 0; x < arglist.size(); x++) {
    otype = arglist[x]->get_object_type();
    dtype = arglist[x]->get_data_type();
    
    switch (otype){
			
    case TREESET: 
      arg = "TREESET";
      break;				
    case LIST:
      arg = "LIST";
      break;
    default:
      switch (dtype){			
      case THE_EMPTY_LIST: 
	arg = "THE_EMPTY_LIST";
	break;
      case BOOLEAN:
	arg = "BOOLEAN";
	break;
      case CHAR:
	arg = "CHAR";
	break;
      case INT:
	arg = "INT";
	break;
      case FLOAT:
	arg = "FLOAT";
	break;
      case DOUBLE:
	arg = "DOUBLE";
	break;
      case STRING:
	arg = "STRING";
	break;
      case VECTOR:
	arg = "VECTOR";
	break;
      }
    }

    argtypes += arg;
    if (x < (arglist.size() - 1)) {
      argtypes += ", ";
    }
  }
  return (argtypes + ")");
}


pqlsymbol * u_get_trees_by_subtree(vector< pqlsymbol * > arglist) {  
	//type check and catch errors and handle any method overloading. 
	if (arglist.size() == 1 && arglist[0]->is_string() ) {
		return new pqlsymbol(get_trees_by_subtree(arglist[0]->get_string() ), NUM_TREES );
	}
	else {
		cout << "get_trees_by_subtree expects 1 string argument." << "Found " << get_arg_types(arglist) << endl;
		return  new pqlsymbol(ERROR, "Type Error");
	}
}

pqlsymbol * u_get_trees_by_taxa(vector< pqlsymbol * > arglist){  
	pqlsymbol * result;

	//type check and catch errors and handle any method overloading. 
	if (arglist.size() == 2 && arglist[0]->is_vect() && arglist[0]->is_string() && arglist[1]->is_vect() && arglist[1]->is_string() ) {

		vector<string> missednames1 = ::biparttable.lm.catchDeclaredTaxa(arglist[0]->get_string_vect());
		vector<string> missednames2 = ::biparttable.lm.catchDeclaredTaxa(arglist[1]->get_string_vect());
		
		if(missednames1.size() > 0 || missednames2.size() > 0){
			string err = "Taxa ";
			for (unsigned int i = 0; i < missednames1.size(); i++){
				err += missednames1[i];
				err += ", ";
			}
			for (unsigned int i = 0; i < missednames2.size(); i++){
				err += missednames2[i];
				err += ", ";
			}
			err += "are contained in no trees";
			
			result = new pqlsymbol(ERROR, err);
		}
		else{
			result = new pqlsymbol(get_trees_by_taxa(arglist[0]->get_string_vect(), arglist[1]->get_string_vect()), NUM_TREES );
		}
	}
	else {	
		cout << "get_trees_by_taxa expects 2 string vector arguments." << "Found " << get_arg_types(std::move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}


string missedTaxaErrorMessage(vector<string> missedtaxanames){
	string errorstring = "The following taxa don't appear in any trees: ";
	for (unsigned int i = 0; i < missedtaxanames.size(); i++){
				errorstring += missedtaxanames[i];
				if (i != missedtaxanames.size()-1)
					errorstring += ", ";
	}
	return errorstring;
}


//should be tested for use
string missedTaxaErrorMessage(vector<string> missedtaxanames1, vector<string> missedtaxanames2 ){
	string errorstring = "The following taxa don't appear in any trees: ";
	for (unsigned int i = 0; i < missedtaxanames1.size(); i++){
				errorstring += missedtaxanames1[i];
				errorstring += ", ";
	}

	for (unsigned int i = 0; i < missedtaxanames2.size(); i++){
				errorstring += missedtaxanames2[i];
				if (i != missedtaxanames2.size()-1)
					errorstring += ", ";
	}
	return errorstring;
}
//set<unsigned int> clade_size_search(vector<int> required, int size)
pqlsymbol * u_clade_size_search(vector< pqlsymbol * > arglist) {
	pqlsymbol * result;
	//make sure we have two arguments: first is a vector of ints, second is an int
	if(arglist.size() != 2) {
		cout << "clade_size_search expects 2 arguments: a StringVec/IntVec and an integer. " << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR,"Type Error");
	}
	else if (arglist[0]->is_atom() && arglist[0]->is_string() && arglist[1]->is_int()) {
		vector<string> temp;
		temp.push_back(arglist[0]->get_string());
		vector<string> missednames = ::biparttable.lm.catchDeclaredTaxa(temp);
		if(missednames.size() > 0){
			result = new pqlsymbol(ERROR, missedTaxaErrorMessage(missednames));
		}
		else{
			result = new pqlsymbol(clade_size_search(temp, arglist[1]->get_int() ), ::NUM_TREES );
		}
	}
	else if (arglist[0]->is_vect() && arglist[0]->is_string() && arglist[1]->is_int()) {
		vector<string> missednames = ::biparttable.lm.catchDeclaredTaxa(arglist[0]->get_string_vect());
		if(missednames.size() > 0){
			result = new pqlsymbol(ERROR, missedTaxaErrorMessage(missednames));
		}
		else{
			result = new pqlsymbol(clade_size_search(arglist[0]->get_string_vect(), arglist[1]->get_int() ), ::NUM_TREES );
		}
		
	}
	else if (arglist[0]->is_vect() && arglist[0]->is_int() && arglist[1]->is_int()) {
		result = new pqlsymbol(clade_size_search(arglist[0]->get_int_vect(), arglist[1]->get_int() ), ::NUM_TREES );
	}
	else {
		cout << "clade_size_search expects a StringVec/Intvec and an Integer" << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

pqlsymbol * u_smallest_clade(vector< pqlsymbol * > arglist) {
	pqlsymbol * result;

	//make sure we have two arguments: first is a vector of ints, second is an int
	if(arglist.size() != 1) {
		result = new pqlsymbol(ERROR,"Type Error");
	}
	else if (arglist[0]->is_atom() && arglist[0]->is_string()) {
		vector<string> temp;
		temp.push_back(arglist[0]->get_string());
		vector<string> missednames = ::biparttable.lm.catchDeclaredTaxa(temp);
		if(missednames.size() > 0){
			result = new pqlsymbol(ERROR, missedTaxaErrorMessage(missednames));
		}
		else{
			result = new pqlsymbol(smallest_clade(temp), ::NUM_TREES );
		}
	}
	else if (arglist[0]->is_vect() && arglist[0]->is_string()) {
		vector<string> missednames = ::biparttable.lm.catchDeclaredTaxa(arglist[0]->get_string_vect());
		if(missednames.size() > 0){
			result = new pqlsymbol(ERROR, missedTaxaErrorMessage(missednames));
		}
		else{
			result = new pqlsymbol(smallest_clade(arglist[0]->get_string_vect()), ::NUM_TREES );
		}
		
	}
	else if (arglist[0]->is_vect() && arglist[0]->is_int()) {
		result = new pqlsymbol(smallest_clade(arglist[0]->get_int_vect()), ::NUM_TREES );
	}
	else {
		cout << "smallest_clade expects a StringVec/Intvec. " << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

pqlsymbol * u_similarity_search(vector< pqlsymbol * > arglist) {

	string stree;

	if (arglist.size() == 1 && arglist[0]->is_string()){
		stree = arglist[0]->get_string();
	}
	else{
		cout << "similarity_search expects a single string representing a tree " << "Found " << get_arg_types(move(arglist)) << endl;
	}
	return new pqlsymbol(similarity_search(stree), ::NUM_TREES);
}

pqlsymbol * u_get_trees_with_taxa(vector< pqlsymbol * > arglist) {  
	pqlsymbol * result;

	//type check and catch errors and handle any method overloading. 
	if (arglist.size() != 1) {
		cout << "get_trees_with_taxa expects either 1 String, 1 StringVect or 1 IntVect argument. " << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR,"Type Error");
	}
	
	else if (arglist[0]->is_atom() && arglist[0]->is_string() ) {
		vector<string> temp;
		temp.push_back(arglist[0]->get_string());
		vector<string> missednames = ::biparttable.lm.catchDeclaredTaxa(temp);
		if(missednames.size() > 0){
			result = new pqlsymbol(ERROR, missedTaxaErrorMessage(missednames));
		}
		else{
			result = new pqlsymbol(get_trees_with_taxa(temp), ::NUM_TREES );
		}
	}
	else if (arglist[0]->is_vect() && arglist[0]->is_string() ) {
		vector<string> missednames = ::biparttable.lm.catchDeclaredTaxa(arglist[0]->get_string_vect());
		
		if(missednames.size() > 0){
			result = new pqlsymbol(ERROR, missedTaxaErrorMessage(missednames));
		}
		else{
			result = new pqlsymbol(get_trees_with_taxa(arglist[0]->get_string_vect() ), ::NUM_TREES );
		}
		
	}
	else if (arglist[0]->is_vect() && arglist[0]->is_int() ) {
		result = new pqlsymbol(get_trees_with_taxa(arglist[0]->get_int_vect() ), ::NUM_TREES );
	}
	else {
		cout << "get_trees_with_taxa expects either 1 StringVect or 1 IntVect argument. " << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}


pqlsymbol * u_get_trees_without_taxa(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	//type check and catch errors and handle any method overloading. 
	if (arglist.size() != 1){
		cout << "get_trees_without_taxa expects either 1 StringVect or 1 IntVect argument. " << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	else if (arglist[0]->is_string() && arglist[0]->is_vect()){
		result = new pqlsymbol(get_trees_without_taxa(arglist[0]->get_string_vect() ), NUM_TREES );
	}
	else if (arglist[0]->is_int() && arglist[0]->is_vect()){
		result = new pqlsymbol(get_trees_without_taxa(arglist[0]->get_int_vect() ), NUM_TREES );
	}
	else{
		cout << "get_trees_without_taxa expects either 1 StringVect or 1 IntVect argument. " << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

pqlsymbol * u_to_newick(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	//type check and catch errors and handle any method overloading. 
	if (arglist.size() != 1) {
		cout << "to_newick expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	else if (arglist[0]->is_treeset()){
		result = new pqlsymbol(to_newick(arglist[0]->get_treeset() ) );
	}
	
	else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == LIST) {
		result = new pqlsymbol(to_newick(arglist[0]->get_int_vect()) );
	}
    else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == ATOM) {
 		result = new pqlsymbol(to_newick(arglist[0]->get_int()) );
    }
	else{
		cout << "to_newick expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

pqlsymbol * u_least_conflict(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;
	if (arglist.size() == 1 && arglist[0]->is_treeset()){
		result = new pqlsymbol(least_conflict(arglist[0]->get_treeset()) );
		//cout << "consensus expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(arglist) << endl;
	}
	else {
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}   
  
      
pqlsymbol * u_greedy_consen(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;
	if (arglist.size() == 1 && arglist[0]->is_treeset()){
		result = new pqlsymbol(greedy_consen(arglist[0]->get_treeset(), 100 ) );
	}
	else if (arglist.size() == 2 && arglist[0]->is_treeset() && arglist[1]->is_double()){
	  result = new pqlsymbol(greedy_consen(arglist[0]->get_treeset(), (float)(arglist[1]->get_double()) ) );
	}
	else if (arglist.size() == 2 && arglist[0]->is_treeset() && arglist[1]->is_int()){
	  result = new pqlsymbol(greedy_consen(arglist[0]->get_treeset(), (float)(arglist[1]->get_int()) ) );
	}
	else {
		//cout << "consensus expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}    

pqlsymbol * u_strict_consen(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;
	
	if (arglist.size() > 0 && arglist[0]->is_treeset() && ::biparttable.is_taxa_homogenious(arglist[0]->get_treeset() ) ){
	
		if (arglist.size() == 1 && arglist[0]->is_treeset()){
			result = new pqlsymbol(consen(arglist[0]->get_treeset(), 100 ) );
		}
		else {
			//cout << "consensus expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(arglist) << endl;
			result = new pqlsymbol(ERROR, "Type Error");
		}
	}
	else{
		result = new pqlsymbol(ERROR, "Consensus can only handle taxa homogenious treesets. Taxa heterogenious trees were found");
	}
	return result;
}


pqlsymbol * u_majority_consen(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;
	
	if (arglist.size() > 0 && arglist[0]->is_treeset() && ::biparttable.is_taxa_homogenious(arglist[0]->get_treeset() ) ){
	
		if (arglist.size() == 1 && arglist[0]->is_treeset()){
			result = new pqlsymbol(consen(arglist[0]->get_treeset(), 50 ) );
		}
		else if (arglist.size() == 2 && arglist[0]->is_treeset() && arglist[1]->is_double()){
			result = new pqlsymbol(consen(arglist[0]->get_treeset(), (float)(arglist[1]->get_double()) ) );
		}
		else if (arglist.size() == 2 && arglist[0]->is_treeset() && arglist[1]->is_int()){
			result = new pqlsymbol(consen(arglist[0]->get_treeset(), (float)(arglist[1]->get_int()) ) );
		}
		else {
			//cout << "consensus expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(arglist) << endl;
			result = new pqlsymbol(ERROR, "Type Error");
		}
	}
	else{
		result = new pqlsymbol(ERROR, "Consensus can only handle taxa homogenious treesets. Taxa heterogenious trees were found");
	}
	return result;
}


pqlsymbol * u_consen(vector<pqlsymbol * > arglist) 
{  
	//pqlsymbol * result;
	
	set<unsigned int> tset;
	float threshold = 50.0;
	
	
	if (arglist.size() > 0 && arglist[0]->is_treeset()){
		tset =  arglist[0]->get_treeset();
	}
	else{
		return new pqlsymbol(ERROR, "Type Error, 1st value");
	}
	
	if (arglist.size() > 1 && arglist[1]->is_int()){
			threshold = (float)arglist[1]->get_int();
	}
	else if (arglist.size() > 1 && arglist[1]->is_double()){
			threshold = (float)arglist[1]->get_double();
	}
	else if (arglist.size() > 1 && arglist[1]->is_float()){
			threshold = (float)arglist[1]->get_float();
	}
	else if (arglist.size() > 1){
		return new pqlsymbol(ERROR, "Type Error, 2nd value");
	}
	
	if(!::biparttable.is_taxa_homogenious(tset) ){
		return new pqlsymbol(ERROR, "Consensus can only handle taxa homogenious treesets. Taxa heterogenious trees were found");
	}
	
	return new pqlsymbol(consen(tset, threshold) );
	
}

pqlsymbol * u_consensus_reso_rate(vector<pqlsymbol *> arglist)
{
	set<unsigned int> tset;
	float threshold = 50.0;

	if (arglist.size() > 0 && arglist[0]->is_treeset()){
		tset = arglist[0]->get_treeset();
	}
	else{

		cout << "u_consensus_reso_rate expects either 1 set of tree or a set of trees and a percentage " << "Found " << get_arg_types(move(arglist)) << endl;
		return new pqlsymbol(ERROR, "Type Error, 1st value");
	}
	
	if (arglist.size() > 1 && arglist[1]->is_int()){
			threshold = (float)arglist[1]->get_int();
	}
	else if (arglist.size() > 1 && arglist[1]->is_double()){
			threshold = (float)arglist[1]->get_double();
	}
	else if (arglist.size() > 1 && arglist[1]->is_float()){
			threshold = (float)arglist[1]->get_float();
	}
	else if (arglist.size() > 1){
		cout << "u_consensus_reso_rate expects either 1 set of tree or a set of trees and a percentage " << "Found " << get_arg_types(move(arglist)) << endl;
		return new pqlsymbol(ERROR, "Type Error, 2nd value");
	}

	if(!::biparttable.is_taxa_homogenious(tset) ){
		return new pqlsymbol(ERROR, "Can only handle taxa homogenious treesets. Taxa heterogenious trees were found");
	}

	return new pqlsymbol(consensus_reso_rate(tset, threshold));
}

pqlsymbol * u_reso_rate(vector<pqlsymbol *> arglist)
{
	string stree;

	if (arglist.size() == 1 && arglist[0]->is_string()){
		stree = arglist[0]->get_string();
	}
	else{
		cout << "reso_rate expects a single string representing a tree " << "Found " << get_arg_types(move(arglist)) << endl;
	}
	return new pqlsymbol(reso_rate(stree));

}

 pqlsymbol * u_HashCS(vector<pqlsymbol * > arglist) {  
	 pqlsymbol * result;
 
	 //type check and catch errors and handle any method overloading. 
	 if (arglist.size() != 1){
		 cout << "consensus expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(arglist) << endl;
		 result = new pqlsymbol(ERROR, "Type Error");
	 }
	 else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == LIST){
		 result = new pqlsymbol(HashCS(arglist[0]->get_int_vect()) );
	 }
     else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == ATOM) {
 		 result = new pqlsymbol(HashCS(arglist[0]->get_int()));
     }
	 else {
		 cout << "consensus expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(arglist) << endl;
		 result = new pqlsymbol(ERROR, "Type Error");
	 }
	 return result;
 }
/*
pqlsymbol * u_search_hashtable_strict(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;
	
	//type check and catch errors and handle any method overloading. 
	
	vector<int> left;
	vector<int> right;
	
	bool LeftIsGood = false;
	bool RightIsGood = false;
	int strict = 0;
	
	if (arglist.size() == 2 || arglist.size() == 3){
		if (arglist[0]->is_string() ){
			if (arglist[0]->is_vect() ){
				LeftIsGood = true;
				left = ::lm.lookUpLabels(arglist[0]->get_string_vect());
			}
			else if ( arglist[0]->is_atom() ){
				LeftIsGood = true;
				left = ::lm.lookUpLabels(arglist[0]->get_string());
			}
		}
		
		else if ( arglist[0]->is_emptylist() ){
			LeftIsGood = true;
		}
		
		if (arglist[1]->is_string() ){
			if (arglist[1]->is_vect() ){
				RightIsGood = true;
				right = ::lm.lookUpLabels(arglist[1]->get_string_vect());
			}
			else if ( arglist[1]->is_atom() ){
				RightIsGood = true;
				right = ::lm.lookUpLabels(arglist[1]->get_string());
			}
		}
		else if ( arglist[1]->is_emptylist() ){
			RightIsGood = true;
		}
	}
		
	if (arglist.size() == 3){
		if (arglist[2]->get_string() == "neither" || arglist[2]->get_string() == "Neither" || arglist[2]->get_string() == "n" || arglist[2]->get_string() == "N" ){
			strict = 0;		
		}
		else if (arglist[2]->get_string() == "left" || arglist[2]->get_string() == "Left" || arglist[2]->get_string() == "l" || arglist[2]->get_string() == "L"){
			strict = 1;	
		}
		else if (arglist[2]->get_string() == "right" || arglist[2]->get_string() == "Right" || arglist[2]->get_string() == "r" || arglist[2]->get_string() == "R"){
			strict = 2;		
		}
		else if (arglist[2]->get_string() == "both" || arglist[2]->get_string() == "Both" || arglist[2]->get_string() == "b" || arglist[2]->get_string() == "B"){
			strict = 3;	
		}
		
	}
	
	if(LeftIsGood && RightIsGood){
		result = new pqlsymbol(search_hashtable_strict_old( left, right, strict ), NUM_TREES );

	}
	else{
		cout << "search_hashtable_strict expects either 2 string vectors or 2 string vectors and an int. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	
	return result;
}
*/
pqlsymbol * u_random_search(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;
	
	int left = 2;
	int right = 2;
	int strict = 0;
	int iterations = 1;
	
	if (arglist.size() >= 2){
		left = arglist[0]->get_int();
		right = arglist[1]->get_int();
	}
	if (arglist.size() >= 3){
		strict = arglist[2]->get_int();
		
	}
	if (arglist.size() >= 4){
		iterations = arglist[3]->get_int();
	}
		
	result = new pqlsymbol(random_search( left, right, strict, iterations ) );
	
	return result;
}


pqlsymbol * u_random_search2(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;
	
	int left = 2;
	int right = 2;
	int strict = 0;
	int iterations = 1;
	
	if (arglist.size() >= 2){
		left = arglist[0]->get_int();
		right = arglist[1]->get_int();
	}
	if (arglist.size() >= 3){
		strict = arglist[2]->get_int();
		
	}
	if (arglist.size() >= 4){
		iterations = arglist[3]->get_int();
	}
		
	result = new pqlsymbol(random_search2( left, right, strict, iterations ) );
	
	return result;
}
/*
pqlsymbol * u_search_hashtable_auto_and_timed(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;
	
	//type check and catch errors and handle any method overloading. 
	
	vector<int> left;
	vector<int> right;
	
	bool LeftIsGood = false;
	bool RightIsGood = false;
	int strict = 0;
	
	if (arglist.size() == 2 || arglist.size() == 3){
		if (arglist[0]->is_int() ){
			if (arglist[0]->is_vect() ){
				LeftIsGood = true;
				left = arglist[0]->get_int_vect();
			}
			else if ( arglist[0]->is_atom() ){
				LeftIsGood = true;
				left.push_back(arglist[0]->get_int());
			}
		}
		
		else if ( arglist[0]->is_emptylist() ){
			LeftIsGood = true;
		}
		
		if (arglist[1]->is_int() ){
			if (arglist[1]->is_vect() ){
				RightIsGood = true;
				right = arglist[1]->get_int_vect();
			}
			else if ( arglist[1]->is_atom() ){
				RightIsGood = true;
				right.push_back(arglist[1]->get_int());
			}
		}
		else if ( arglist[1]->is_emptylist() ){
			RightIsGood = true;
		}
	}
		
	if (arglist.size() == 3){
		if (arglist[2]->get_string() == "neither" || arglist[2]->get_string() == "Neither" || arglist[2]->get_string() == "n" || arglist[2]->get_string() == "N" ){
			strict = 0;		
		}
		else if (arglist[2]->get_string() == "left" || arglist[2]->get_string() == "Left" || arglist[2]->get_string() == "l" || arglist[2]->get_string() == "L"){
			strict = 1;	
		}
		else if (arglist[2]->get_string() == "right" || arglist[2]->get_string() == "Right" || arglist[2]->get_string() == "r" || arglist[2]->get_string() == "R"){
			strict = 2;		
		}
		else if (arglist[2]->get_string() == "both" || arglist[2]->get_string() == "Both" || arglist[2]->get_string() == "b" || arglist[2]->get_string() == "B"){
			strict = 3;	
		}
	}
	
	if(LeftIsGood && RightIsGood){
		result = new pqlsymbol(search_hashtable_strict_and_timed( left, right, strict ), NUM_TREES );
	}
	else{
		cout << "consensus expects either 2 string vectors or 2 string vectors and an int. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	
	return result;
}
*/
pqlsymbol * u_search_hashtable_strict_and_timed(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;
	
	//type check and catch errors and handle any method overloading. 
	
	vector<int> left;
	vector<int> right;
	
	bool LeftIsGood = false;
	bool RightIsGood = false;
	int strict = 0;
	
	if (arglist.size() == 2 || arglist.size() == 3){
		if (arglist[0]->is_string() ){
			if (arglist[0]->is_vect() ){
				LeftIsGood = true;
				left = ::biparttable.lm.lookUpLabels(arglist[0]->get_string_vect());
			}
			else if ( arglist[0]->is_atom() ){
				LeftIsGood = true;
				left = ::biparttable.lm.lookUpLabels(arglist[0]->get_string());
			}
		}
		
		else if ( arglist[0]->is_emptylist() ){
			LeftIsGood = true;
		}
		
		if (arglist[1]->is_string() ){
			if (arglist[1]->is_vect() ){
				RightIsGood = true;
				right = ::biparttable.lm.lookUpLabels(arglist[1]->get_string_vect());
			}
			else if ( arglist[1]->is_atom() ){
				RightIsGood = true;
				right = ::biparttable.lm.lookUpLabels(arglist[1]->get_string());
			}
		}
		else if ( arglist[1]->is_emptylist() ){
			RightIsGood = true;
		}
	}
		
	if (arglist.size() == 3){
		if (arglist[2]->get_string() == "neither" || arglist[2]->get_string() == "Neither" || arglist[2]->get_string() == "n" || arglist[2]->get_string() == "N" ){
			strict = 0;		
		}
		else if (arglist[2]->get_string() == "left" || arglist[2]->get_string() == "Left" || arglist[2]->get_string() == "l" || arglist[2]->get_string() == "L"){
			strict = 1;	
		}
		else if (arglist[2]->get_string() == "right" || arglist[2]->get_string() == "Right" || arglist[2]->get_string() == "r" || arglist[2]->get_string() == "R"){
			strict = 2;		
		}
		else if (arglist[2]->get_string() == "both" || arglist[2]->get_string() == "Both" || arglist[2]->get_string() == "b" || arglist[2]->get_string() == "B"){
			strict = 3;	
		}
	}
	
	if(LeftIsGood && RightIsGood){
		result = new pqlsymbol(search_hashtable_strict_and_timed( left, right, strict ), NUM_TREES );
	}
	else{
		cout << "consensus expects either 2 string vectors or 2 string vectors and an int. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	
	return result;
}

pqlsymbol * u_new_search_hashtable_strict(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;
	
	//type check and catch errors and handle any method overloading. 
	
	vector<int> left;
	vector<int> right;
	
	bool LeftIsGood = false;
	bool RightIsGood = false;
	int strict = 0;
	
	if (arglist.size() == 2 || arglist.size() == 3){
		if (arglist[0]->is_string() ){
			if (arglist[0]->is_vect() ){
				LeftIsGood = true;
				left = ::biparttable.lm.lookUpLabels(arglist[0]->get_string_vect());
			}
			else if ( arglist[0]->is_atom() ){
				LeftIsGood = true;
				left = ::biparttable.lm.lookUpLabels(arglist[0]->get_string());
			}
		}
		
		else if ( arglist[0]->is_emptylist() ){
			LeftIsGood = true;
		}
		
		if (arglist[1]->is_string() ){
			if (arglist[1]->is_vect() ){
				RightIsGood = true;
				right = ::biparttable.lm.lookUpLabels(arglist[1]->get_string_vect());
			}
			else if ( arglist[1]->is_atom() ){
				RightIsGood = true;
				right = ::biparttable.lm.lookUpLabels(arglist[1]->get_string());
			}
		}
		else if ( arglist[1]->is_emptylist() ){
			RightIsGood = true;
		}
	}
		
	if (arglist.size() == 3){
		if (arglist[2]->get_string() == "neither" || arglist[2]->get_string() == "Neither" || arglist[2]->get_string() == "n" || arglist[2]->get_string() == "N" ){
			strict = 0;		
		}
		else if (arglist[2]->get_string() == "left" || arglist[2]->get_string() == "Left" || arglist[2]->get_string() == "l" || arglist[2]->get_string() == "L"){
			strict = 1;	
		}
		else if (arglist[2]->get_string() == "right" || arglist[2]->get_string() == "Right" || arglist[2]->get_string() == "r" || arglist[2]->get_string() == "R"){
			strict = 2;		
		}
		else if (arglist[2]->get_string() == "both" || arglist[2]->get_string() == "Both" || arglist[2]->get_string() == "b" || arglist[2]->get_string() == "B"){
			strict = 3;	
		}
	}
	
	if(LeftIsGood && RightIsGood){
		result = new pqlsymbol(search_hashtable_strict( left, right, strict ), NUM_TREES );
	}
	else{
		cout << "consensus expects either 2 string vectors or 2 string vectors and an int. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	
	return result;
}


pqlsymbol * u_quick_quartet(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	//type check and catch errors and handle any method overloading. 
	if (arglist.size() != 1){
		cout << "quick quartet expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == LIST){
		result = new pqlsymbol(quick_quartet(arglist[0]->get_int_vect()) );
	}
	else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == ATOM){
 		result = new pqlsymbol(quick_quartet(arglist[0]->get_int()));
 	}
	else if (arglist[0]->is_treeset() ){
 		result = new pqlsymbol(quick_quartet(arglist[0]->get_treeset() ) );
 	}
	else{
		cout << "quick quartet expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}


pqlsymbol * u_phlash(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	//type check and catch errors and handle any method overloading. 
	if (arglist.size() != 1){
		cout << "phlash expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == LIST){
		result = new pqlsymbol(phlash(arglist[0]->get_int_vect()) );
	}
    else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == ATOM){
 	result = new pqlsymbol(phlash(arglist[0]->get_int()));
    }
	else{
		cout << "phlash expects either 1 IntVect or 1 Int argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

int help( ) 
{
	cout << "add\n" << endl;
	cout << "get_trees_without_taxa\n" << endl;
	cout << "get_trees_with_taxa\n" << endl;
	cout << "get_trees_by_taxa\n" << endl;
	cout << "to_newick\n" << endl;
	cout << "get_trees_by_subtree\n" << endl;
	cout << "consensus\n" << endl;
	cout << "search_by_relationship\n" << endl;
	cout << "qq: Returns a quick quartet matrix of the selected trees "
			"and saves the output into qq.txt found in the temp directory.\n" << endl;
	cout << "phlash\n" << endl;
	
	return 0;
}

pqlsymbol * u_help(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;
	
	result = new pqlsymbol(help( ));

	return result;
}

pqlsymbol * u_unique(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	if (arglist.size() == 1 && arglist[0]->is_treeset() ){
		result = new pqlsymbol(unique(arglist[0]->get_treeset()), NUM_TREES );
	}
	else{
		cout << "unique expects 1 treeset argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}

	return result;
}

//unique_biparts(set< unsigned int > treesin)
pqlsymbol * u_unique_biparts(vector<pqlsymbol *> arglist)
{
	pqlsymbol * result;

	if(arglist.size() != 1){	
		cout << "Error: unique_biparts expects one argument of type tree set";
		result = new pqlsymbol(ERROR, "Type Error");
	}
	else if(!arglist[0]->is_treeset() ){
		cout << "Error: unique_biparts takes a tree set, Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	else{
		result = new pqlsymbol(unique_biparts(arglist[0]->get_treeset()));
		//result = new pqlsymbol(ERROR, "Type Error");
	}

	return result;
}

pqlsymbol * u_silhouette(vector<pqlsymbol * > arglist) {
	pqlsymbol * result = new pqlsymbol();

	if(arglist.size() == 1 && arglist[0]->is_vect()){
		result = new pqlsymbol(silhouette(arglist[0]->get_treeset_vect()));
	}
	else{
		cout << "silhouette expect a single treeset vector as input " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}	
	return result;

}

pqlsymbol * u_agglo_clust(vector<pqlsymbol * > arglist){
	pqlsymbol * result = new pqlsymbol();
	if (arglist.size() == 1 && arglist[0]->is_treeset()){
		result = new pqlsymbol(agglo_clust(arglist[0]->get_treeset()));
	}
	else{
		cout << "agglo_clust expects a single treeset as input " << " Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}

	return result;

	}		

// takes an int. Returns the ints which share the same topology. 
pqlsymbol * u_duplicates(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	if (arglist.size() == 1 && arglist[0]->is_int() &&  arglist[0]->is_atom()){
		result = new pqlsymbol(duplicates(arglist[0]->get_int() ),NUM_TREES );
	}
	else{
		cout << "unique expects 1 Int argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}

	return result;
}


pqlsymbol * u_count(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	if (arglist.size() == 1){
		result = new pqlsymbol(arglist[0]->get_size());
	}
	else{
		cout << "count expects 1 argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}

	return result;
}

pqlsymbol * u_union(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	if (arglist.size() == 2 && arglist[0]->is_treeset() &&  arglist[1]->is_treeset()){
		set<unsigned int> s1 = arglist[0]->get_treeset();
		set<unsigned int> s2 = arglist[1]->get_treeset();
		set<unsigned int> sunion; 
		
		std::set_union( s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter( sunion, sunion.begin() ) );
                    
        result = new pqlsymbol(sunion, NUM_TREES);
		
	}
	else{
		cout << "union expects 2 treeset arguments. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}

	return result;
}


pqlsymbol * u_difference(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	if (arglist.size() == 2 && arglist[0]->is_treeset() &&  arglist[1]->is_treeset()){
		set<unsigned int> s1 = arglist[0]->get_treeset();
		set<unsigned int> s2 = arglist[1]->get_treeset();
		set<unsigned int> sdiff; 
		
		std::set_difference( s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter( sdiff, sdiff.begin() ) );
                    
        result = new pqlsymbol(sdiff, NUM_TREES);
		
	}
	else{
		cout << "difference expects 2 treeset arguments. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}

	return result;
}


pqlsymbol * u_intersection(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	if (arglist.size() == 2 && arglist[0]->is_treeset() &&  arglist[1]->is_treeset()){
		set<unsigned int> s1 = arglist[0]->get_treeset();
		set<unsigned int> s2 = arglist[1]->get_treeset();
		set<unsigned int> sinter; 
		
		std::set_intersection( s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter( sinter, sinter.begin() ) );
                    
        result = new pqlsymbol(sinter, NUM_TREES);
		
	}
	else{
		cout << "interesection expects 2 treeset arguments. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}

	return result;
}

pqlsymbol * u_not(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	if (arglist.size() == 1 && arglist[0]->is_treeset()){
		set<unsigned int> s1 = arglist[0]->get_treeset();
		set<unsigned int> sdiff; 
		
		std::set_difference( all_trees.begin(), all_trees.end(), s1.begin(), s1.end(),  std::inserter( sdiff, sdiff.begin() ) );
                    
        result = new pqlsymbol(sdiff, NUM_TREES);
		
	}
	else{
		cout << "not expects 1 treeset argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}

	return result;
}

pqlsymbol * u_set(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;

	if (arglist.size() == 1 && arglist[0]->is_int() && arglist[0]->is_vect()){
		vector<int> tempvect = arglist[0]->get_int_vect(); 
		set<unsigned int> s1(tempvect.begin(), tempvect.end());

        result = new pqlsymbol(s1, NUM_TREES);
	}
	else{
		cout << "set expects 1 intvect argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

pqlsymbol * u_show(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result = new pqlsymbol();
  string mode = "";
  system("echo \"<table border=\\\"1\\\">\" > temp/svg.html");
  
  //type check and catch errors and handle any method overloading. 
  if (((arglist.size() != 1) && (arglist.size() != 2)) || (arglist.size() == 2 && !(arglist[1]->is_string()))) {
    cout << "show expects one IntVect or Int argument, and one optional string argument: \"text\", \"ortho\", or \"radial\". " << "Found " << get_arg_types(move(arglist)) << endl;
    return result = new pqlsymbol(ERROR, "Type Error");
  }
  if (arglist.size() == 1) {
    mode = "text";
  } else if (arglist[1]->get_string() != "text" && arglist[1]->get_string() != "ortho" && arglist[1]->get_string() != "radial") {
    cout << "Mode argument must be one of: \"text\", \"ortho\", or \"radial\"; found \"" << arglist[1]->get_string() << "\"" << endl;
    return result = new pqlsymbol(ERROR, "Invalid Argument");
  } else {
    mode = arglist[1]->get_string();
  }
  if (arglist[0]->is_treeset()) {
    show(arglist[0]->get_treeset(), mode);
  } else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == LIST) {
    show(arglist[0]->get_int_vect(), mode);
  } else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == ATOM) {
    show(arglist[0]->get_int(), mode);
  } else {
   cout << "show expects one IntVect or Int argument, and one optional string argument: \"text\", \"ortho\", or \"radial\". " << "Found " << get_arg_types(move(arglist)) << endl;
   result = new pqlsymbol(ERROR, "Type Error");
  }
  if (mode == "ortho" || mode == "radial") {
    system("xdg-open temp/svg.html");
  } else {
    system("rm temp/svg.html");
  }
  return result;
}

pqlsymbol * u_show_newick(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result = new pqlsymbol();
  string mode = "";
  system("echo \"<table border=\\\"1\\\">\" > temp/svg.html"); 
  if (arglist.size() != 1 && (arglist.size() != 2) || ((arglist.size() == 2) && !(arglist[1]->is_string()))) {
    cout << "show_newick expects one String or StringVect argument, and one optional string argument: \"text\", \"ortho\", or \"radial\". " << "Found " << get_arg_types(move(arglist)) << endl;
    return result = new pqlsymbol(ERROR, "Type Error");
  }
  if (arglist.size() == 1) {
    mode = "text";
  } else if (arglist[1]->get_string() != "text" && arglist[1]->get_string() != "ortho" && arglist[1]->get_string() != "radial") {
    cout << "Mode argument must be one of: \"text\", \"ortho\", or \"radial\"; found \"" << arglist[1]->get_string() << "\"" << endl;
    return result = new pqlsymbol(ERROR, "Invalid Argument");
  } else {
    mode = arglist[1]->get_string();
  }
  if (arglist[0]->get_data_type() == STRING && arglist[0]->get_object_type() == LIST) {
    show_newick(arglist[0]->get_string_vect(), mode);
  } else if (arglist[0]->get_data_type() == STRING && arglist[0]->get_object_type() == ATOM) {
    show_newick(arglist[0]->get_string(), "Newick Tree", "",  mode);
  } else {
   cout << "show_newick expects one String or StringVect argument, and one optional string argument: \"text\", \"ortho\", or \"radial\". " << "Found " << get_arg_types(move(arglist)) << endl;
   result = new pqlsymbol(ERROR, "Type Error");
  }
  if (mode == "ortho" || mode == "radial") {
    system("xdg-open temp/svg.html");
  } else {
    system("rm temp/svg.html");
  }
  return result;
}

pqlsymbol * u_export(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result = new pqlsymbol();
  string mode;
  //type check and catch errors and handle any method overloading. 
  if ((arglist.size() != 3 && arglist.size() != 2) || !(arglist[1]->is_string()) || (arglist.size() == 3 && !(arglist[2]->is_string()))) {
    cout << "export expects one IntVect or Int argument, one filesystem path, and one optional String argument: \"text\", \"ortho\", or \"radial\". " << "Found " << get_arg_types(move(arglist)) << endl;
    return result = new pqlsymbol(ERROR, "Type Error");
  }
  if (arglist.size() == 2) {
    mode = "text";
  } else if (arglist[2]->get_string() != "text" && arglist[2]->get_string() != "ortho" && arglist[2]->get_string() != "radial") {
    cout << "Mode argument must be one of: \"text\", \"ortho\", or \"radial\"; found \"" << arglist[2]->get_string() << "\"" << endl;
    return result = new pqlsymbol(ERROR, "Invalid Argument");
  } else {
    mode = arglist[2]->get_string();
  }
  if (arglist[0]->is_treeset()){
    exportt(arglist[0]->get_treeset(), arglist[1]->get_string(), mode);
  } else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == LIST) {
    exportt(arglist[0]->get_int_vect(), arglist[1]->get_string(), mode);
  } else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == ATOM) {
    exportt(arglist[0]->get_int(), arglist[1]->get_string(), mode);
  } else {
    cout << "export expects one IntVect or Int argument, one filesystem path, and one optional String argument: \"text\", \"ortho\", or \"radial\". " << "Found " << get_arg_types(move(arglist)) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}

pqlsymbol * u_export_newick(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result = new pqlsymbol();
  string mode;
  //type check and catch errors and handle any method overloading. 
  if ((arglist.size() != 3 && arglist.size() != 2) || !(arglist[1]->is_string()) || (arglist.size() == 3 && !(arglist[2]->is_string()))) {
    cout << "export_newick expects one String or StringVect argument, one filesystem path, and one optional String argument: \"text\", \"ortho\", or \"radial\". " << "Found " << get_arg_types(move(arglist)) << endl;
    return result = new pqlsymbol(ERROR, "Type Error");
  }
  if (arglist.size() == 2) {
    mode = "text";
  } else if (arglist[2]->get_string() != "text" && arglist[2]->get_string() != "ortho" && arglist[2]->get_string() != "radial") {
    cout << "Mode argument must be one of: \"text\", \"ortho\", or \"radial\"; found \"" << arglist[2]->get_string() << "\"" << endl;
    return result = new pqlsymbol(ERROR, "Invalid Argument");
  } else {
    mode = arglist[2]->get_string();
  }
  if (arglist[0]->get_data_type() == STRING && arglist[0]->get_object_type() == LIST) {
    export_newick(arglist[0]->get_string_vect(), arglist[1]->get_string(), mode);
  } else if (arglist[0]->get_data_type() == STRING && arglist[0]->get_object_type() == ATOM) {
    export_newick(arglist[0]->get_string(), "nwtree", arglist[1]->get_string(), "", mode);
  } else {
    cout << "export_newick expects one String or StringVect argument, one filesystem path, and one optional String argument: \"text\", \"ortho\", or \"radial\". " << "Found " << get_arg_types(move(arglist)) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}

pqlsymbol * u_classification(vector<pqlsymbol * > arglist) {
  pqlsymbol * result = new pqlsymbol();
  if (arglist.size() != 1) {
    cout << "classification expects either 1 String or 1 StringVect argument. " << "Found " << get_arg_types(move(arglist)) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  } else if (arglist[0]->get_data_type() == STRING && arglist[0]->get_object_type() == ATOM) {
    classification(arglist[0]->get_string() );
  } else if (arglist[0]->get_data_type() == STRING && arglist[0]->get_object_type() == LIST) {
    classification(arglist[0]->get_string_vect());
  } else {
    cout << "classification expects either 1 String or 1 StringVect argument. " << "Found " << get_arg_types(move(arglist)) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}

pqlsymbol * u_write_classifications(vector<pqlsymbol * > arglist) {
  pqlsymbol * result = new pqlsymbol();
  if (arglist.size() != 1) {
    cout << "write_classifications expects 1 String argument. " << "Found " << get_arg_types(move(arglist)) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  } else if (arglist[0]->is_string()) {
    write_classifications(arglist[0]->get_string() );
  } else {
    cout << "write_classifications expects 1 String argument. " << "Found " << get_arg_types(move(arglist)) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}

pqlsymbol * u_edit_classification(vector<pqlsymbol * > arglist) {
  pqlsymbol * result = new pqlsymbol();
  if (arglist.size() != 3 || !arglist[1]->is_string() || !arglist[2]->is_string()) {
    cout << "classification expects either 3 String arguments or 1 StringVect and 2 String arguments. " << "Found " << get_arg_types(move(arglist)) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  } else if (arglist[0]->get_data_type() == STRING && arglist[0]->get_object_type() == ATOM) {
    edit_classification(arglist[0]->get_string(), arglist[1]->get_string(), arglist[2]->get_string() );
  } else if (arglist[0]->get_data_type() == STRING && arglist[0]->get_object_type() == LIST) {
    edit_classification(arglist[0]->get_string_vect(), arglist[1]->get_string(), arglist[2]->get_string() );
  } else {
    cout << "classification expects either 3 String arguments or 1 StringVect and 2 String arguments. " << "Found " << get_arg_types(move(arglist)) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}

pqlsymbol * u_taxa_in_group(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result;
  if (arglist.size() == 1 && arglist[0]->is_string() ){
    result = new pqlsymbol(taxa_in_group(arglist[0]->get_string()));
  }
  else {
    cout << "taxa_in_group expects one String argument. " << "Found " << get_arg_types(arglist) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}

pqlsymbol * u_show_group_in_tree(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result = new pqlsymbol();
  string mode;
  system("echo \"<table border=\\\"1\\\">\" > temp/svg.html");
  if ((arglist.size() != 2 && arglist.size() != 3) || !(arglist[0]->is_string()) || !(arglist[1]->get_data_type() == INT || arglist[1]->is_treeset())) {
    cout << "show_group_in_tree expects one String argument (group), one Int argument (treenum), and one optional String argument (mode). " << "Found " << get_arg_types(move(arglist)) << endl;
    return result = new pqlsymbol(ERROR, "Type Error");
  } 
  if (arglist.size() == 2) {
    mode = "ortho";
  }
  if (arglist.size() == 3 && arglist[2]->is_string()) {
	  
    mode = arglist[2]->get_string();
    std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);

    if (mode != "ortho" && mode != "radial") {
      cout << "Invalid mode argument. Mode must be either 'ortho' or 'radial'" << endl;
      return result = new pqlsymbol(ERROR, "Invalid Argument");
    }
  }
  int ret; //return of highlight_group
  if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == ATOM) {
    vector<int> int_in_vect; 
    int_in_vect.push_back(arglist[1]->get_int()); //wrap the int argument in a vector (simplifies passing argument)
    ret = highlight_group(arglist[0]->get_string(), int_in_vect, mode);
  } else if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == LIST) {
    ret = highlight_group(arglist[0]->get_string(), arglist[1]->get_int_vect(), mode);
  } else if (arglist[1]->is_treeset()) {
    ret = highlight_group(arglist[0]->get_string(), arglist[1]->get_treeset(), mode);
  }
  if (ret == 0)
    system("xdg-open temp/svg.html");
  return result;
}

pqlsymbol * u_show_group(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result = new pqlsymbol();
  string mode;
  system("echo \"<table border=\\\"1\\\">\" > temp/svg.html");
  if ((arglist.size() != 2 && arglist.size() != 3) || !(arglist[0]->is_string()) || !(arglist[1]->get_data_type() == INT || arglist[1]->is_treeset())) {
    cout << "show_group expects one String argument (group), one Int argument (treenum), and one optional String argument (mode). " << "Found " << get_arg_types(move(arglist)) << endl;
    return result = new pqlsymbol(ERROR, "Type Error");
  } 
  if (arglist.size() == 2) {
    mode = "ortho";
  }
  if (arglist.size() == 3 && arglist[2]->is_string()) {
    mode = arglist[2]->get_string();
    std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
    if (mode != "ortho" && mode != "radial") {
      cout << "Invalid mode argument. Mode must be either 'ortho' or 'radial'" << endl;
      return result = new pqlsymbol(ERROR, "Invalid Argument");
    }
  }
  int ret; //return of show_group
  if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == ATOM) {
    vector<int> int_in_vect; 
    int_in_vect.push_back(arglist[1]->get_int()); //wrap the int argument in a vector (simplifies passing argument)
    ret = highlight_group(arglist[0]->get_string(), int_in_vect, mode, true);
  } else if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == LIST) {
    ret = highlight_group(arglist[0]->get_string(), arglist[1]->get_int_vect(), mode, true);
  } else if (arglist[1]->is_treeset()) {
    ret = highlight_group(arglist[0]->get_string(), arglist[1]->get_treeset(), mode, true);
  }
  if (ret == 0)
    system("xdg-open temp/svg.html");
  return result;
}

pqlsymbol * u_show_level(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result = new pqlsymbol();
  string mode;
  system("echo \"<table border=\\\"1\\\">\" > temp/svg.html");
  if ((arglist.size() != 2 && arglist.size() != 3) || !(arglist[0]->is_string()) || !(arglist[1]->get_data_type() == INT || arglist[1]->is_treeset())) {
    cout << "show_level expects one String argument (group), one Int argument (treenum), and one optional String argument (mode). " << "Found " << get_arg_types(move(arglist)) << endl;
    return result = new pqlsymbol(ERROR, "Type Error");
  } 
  if (arglist.size() == 2) {
    mode = "text";
  }
  if (arglist.size() == 3 && arglist[2]->is_string()) {
    mode = arglist[2]->get_string();
    std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);
    if (mode != "text" && mode != "ortho" && mode != "radial") {
      cout << "Invalid mode argument. Mode must be either 'text', 'ortho', or 'radial'" << endl;
      return result = new pqlsymbol(ERROR, "Invalid Argument");
    }
  }
  int ret; //return of show_level
  if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == ATOM) {
    ret = show_level(arglist[0]->get_string(), arglist[1]->get_int(), mode);
  } else if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == LIST) {
    ret = show_level(arglist[0]->get_string(), arglist[1]->get_int_vect(), mode);
  } else if (arglist[1]->is_treeset()) {
    ret = show_level(arglist[0]->get_string(), arglist[1]->get_treeset(), mode);
  }
  if (ret == 0 && mode != "text")
    system("xdg-open temp/svg.html");
  return result;
}

pqlsymbol * u_show_only(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result = new pqlsymbol();
  string mode;
  vector<string> groups;
  system("echo \"<table border=\\\"1\\\">\" > temp/svg.html");
  if ((arglist.size() != 2 && arglist.size() != 3) || !(arglist[0]->get_data_type() == STRING) || !(arglist[1]->get_data_type() == INT || arglist[1]->is_treeset())) {
    cout << "show_only expects one String or StringVect argument (group), one Int or IntVect argument (treenum), and one optional String argument (mode). " << "Found " << get_arg_types(move(arglist)) << endl;
    return result = new pqlsymbol(ERROR, "Type Error");
  } 
  if (arglist.size() == 2) {
    mode = "ortho";
  }
  if (arglist.size() == 3 && arglist[2]->is_string()) {
    mode = arglist[2]->get_string();
    std::transform(mode.begin(), mode.end(), mode.begin(), ::tolower);

    if (mode != "ortho" && mode != "radial") {
      cout << "Invalid mode argument. Mode must be either 'ortho' or 'radial'" << endl;
      return result = new pqlsymbol(ERROR, "Invalid Argument");
    }
  }
  int ret; //return of show_group
  if (arglist[0]->get_object_type() == ATOM) {
    groups.push_back(arglist[0]->get_string()); //wrap the string argument in a vector
  } else {
    groups = arglist[0]->get_string_vect();
  } 
  if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == ATOM) {
    vector<int> int_in_vect; 
    int_in_vect.push_back(arglist[1]->get_int()); //wrap the int argument in a vector
    ret = show_only(groups, int_in_vect, mode);
  } else if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == LIST) {
    ret = show_only(groups, arglist[1]->get_int_vect(), mode);
  } else if (arglist[1]->is_treeset()) {
    ret = show_only(groups, arglist[1]->get_treeset(), mode);
  }
  if (ret == 0)
    system("xdg-open temp/svg.html");
  return result;
}

pqlsymbol * u_load_trait_matrix(vector<pqlsymbol * > arglist) {
  pqlsymbol * result = new pqlsymbol();
  if (arglist.size() != 1) {
    cout << "load_trait_matrix expects 1 String argument. " << "Found " << get_arg_types(move(arglist)) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  } else if (arglist[0]->is_string()) {
    load_trait_matrix(arglist[0]->get_string() );
  } else {
    cout << "load_trait_matrix expects 1 String argument. " << "Found " << get_arg_types(move(arglist)) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}

pqlsymbol * u_test_trait_correlation(vector<pqlsymbol * > arglist) {
  pqlsymbol * result = new pqlsymbol();
  if (arglist.size() != 7) {
    cout << "test_trait_correlation expects 6 int arguments (t1ind, t1val, t2ind, t2val, tree, iterations) and 1 string argument (folder)." << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  } else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == ATOM
	     && arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == ATOM
	     && arglist[2]->get_data_type() == INT && arglist[2]->get_object_type() == ATOM
	     && arglist[3]->get_data_type() == INT && arglist[3]->get_object_type() == ATOM
	     && arglist[4]->get_data_type() == INT && arglist[4]->get_object_type() == ATOM
	     && arglist[5]->get_data_type() == INT && arglist[5]->get_object_type() == ATOM
	     && arglist[6]->is_string()) {
    test_trait_correlation(arglist[0]->get_int(), arglist[1]->get_int(), arglist[1]->get_int(), arglist[3]->get_int(), arglist[4]->get_int(), arglist[5]->get_int(), arglist[6]->get_string());
  } else {
    cout << "test_trait_correlation expects 6 int arguments (t1ind, t1val, t2ind, t2val, tree, iterations) and 1 string argument (folder)." << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}
/*
pqlsymbol * u_taxa_filter(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result = new pqlsymbol();
  vector<string> taxavect;
  if (arglist.size() != 2 || arglist[0]->get_data_type() != STRING) {
    cout << "taxa_filter expects one String or StringVect argument and one Int, Intvect, or Treeset argument. " << "Found " << get_arg_types(arglist) << endl;
   return result = new pqlsymbol(ERROR, "Type Error");
  }
  if (arglist[0]->get_object_type() == ATOM)
    taxavect.push_back(arglist[0]->get_string());
  else
    taxavect = arglist[0]->get_string_vect();
  if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == ATOM) {
    taxa_filter(taxavect, arglist[1]->get_int());
  }
  else if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == LIST) {
    taxa_filter(taxavect, arglist[1]->get_int_vect());
  }
  else if (arglist[1]->is_treeset()) {
    taxa_filter(taxavect, arglist[1]->get_treeset());
  }
  else {
    cout << "taxa_filter expects one String or StringVect argument and one Int, Intvect, or Treeset argument. " << "Found " << get_arg_types(arglist) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}
*/
/*
pqlsymbol * u_group_filter(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result = new pqlsymbol();
  vector<string> taxavect;
  if (arglist.size() != 2 || arglist[0]->get_data_type() != STRING) {
    cout << "group_filter expects one String or StringVect argument and one Int, Intvect, or Treeset argument. " << "Found " << get_arg_types(arglist) << endl;
   return result = new pqlsymbol(ERROR, "Type Error");
  }
  if (arglist[0]->get_object_type() == ATOM)
    taxavect.push_back(arglist[0]->get_string());
  else
    taxavect = arglist[0]->get_string_vect();
  if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == ATOM) {
    group_filter(taxavect, arglist[1]->get_int());
  }
  else if (arglist[1]->get_data_type() == INT && arglist[1]->get_object_type() == LIST) {
    group_filter(taxavect, arglist[1]->get_int_vect());
  }
  else if (arglist[1]->is_treeset()) {
    group_filter(taxavect, arglist[1]->get_treeset());
  }
  else {
    cout << "group_filter expects one String or StringVect argument and one Int, Intvect, or Treeset argument. " << "Found " << get_arg_types(arglist) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}
*/
/*
pqlsymbol * u_delete_tree(vector<pqlsymbol * > arglist) {
  pqlsymbol * result = new pqlsymbol();
  if (arglist.size() != 1) {
    cout << "delete_tree expects one Int, Intvect, or Treeset argument. " << "Found " << get_arg_types(arglist) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == ATOM) {
    delete_tree(arglist[0]->get_int());
  }
  else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == LIST) {
    delete_tree(arglist[0]->get_int_vect());
  }
  else if (arglist[0]->is_treeset()) {
    delete_tree(arglist[0]->get_treeset());
  }
  else {
    cout << "delete_tree expects one Int, Intvect, or Treeset argument. " << "Found " << get_arg_types(arglist) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}
*/
pqlsymbol * u_write_trz(vector<pqlsymbol * > arglist) {  
  pqlsymbol * result = new pqlsymbol();
  if (arglist.size() != 2 || !(arglist[1]->is_string())) {
    cout << "write_trz expects one IntVect, Int, or Treeset argument, and one file path. " << "Found " << get_arg_types(move(arglist)) << endl;
    result = new pqlsymbol(ERROR, "Type Error");
  } else if (arglist[0]->is_treeset()) {
    write_trz(arglist[0]->get_treeset(), arglist[1]->get_string());
  } else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == LIST) {
    write_trz(arglist[0]->get_int_vect(), arglist[1]->get_string());
  } else if (arglist[0]->get_data_type() == INT && arglist[0]->get_object_type() == ATOM) {
    write_trz(arglist[0]->get_int(), arglist[1]->get_string());
  } else {
   cout << "write_trz expects one IntVect, Int, or Treeset argument, and one file path. " << "Found "  << "Found " << get_arg_types(move(arglist)) << endl;
   result = new pqlsymbol(ERROR, "Type Error");
  }
  return result;
}


pqlsymbol * u_debug(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;
	if (arglist.size() == 1 && arglist[0]->is_bool() ){
		
        if (arglist[0]->get_bool()){
			::DEBUGMODE = true;
			result = new pqlsymbol("Debug set to true");
		}
		else{
			::DEBUGMODE = false;
			result = new pqlsymbol("Debug set to false");
		}		
	}
	else{
		cout << "debug expects 1 boolean argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}


pqlsymbol * u_set_hetero(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;
	if (arglist.size() == 1 && arglist[0]->is_bool() ){
		
        if (arglist[0]->get_bool()){
			::HETERO = true;
			result = new pqlsymbol("HETERO set to true");
		}
		else{
			::HETERO = false;
			result = new pqlsymbol("HETERO set to false");
		}		
	}
	else{
		cout << "set_hetero expects 1 boolean argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

pqlsymbol * u_print_taxa_in_trees(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;
	result = new pqlsymbol("Printing taxa In trees");
	::biparttable.print_taxa_in_trees();
	return result;
}

pqlsymbol * u_get_taxa_in_tree(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;
	if (arglist[0]->is_int()){
	result = new pqlsymbol(biparttable.get_taxa_in_tree(arglist[0]->get_int()));
	}
	else {
		cout << "get_taxa_in_tree expects one Int"  << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

pqlsymbol * u_print_biparttable(vector<pqlsymbol * > arglist) {
  pqlsymbol * result = new pqlsymbol();
  ::biparttable.print_biparttable();
  return result;
}


pqlsymbol * u_print_label_map(vector<pqlsymbol * > arglist){
  pqlsymbol * result = new pqlsymbol();
  ::biparttable.lm.printMap();
  return result;
}


pqlsymbol * u_tier(vector<pqlsymbol * > arglist){
	pqlsymbol * result = new pqlsymbol();
	if (arglist[0]->is_int()){
		printBipartTier(arglist[0]->get_int());
	}
    return result;
}

pqlsymbol * u_rate_tiers(vector<pqlsymbol * > arglist){
  pqlsymbol * result = new pqlsymbol();
	if (arglist[0]->is_int()){
		rateTierRogueness(arglist[0]->get_int());
	}
  return result;
}

pqlsymbol * u_tttier(vector<pqlsymbol * > arglist){
    pqlsymbol * result = new pqlsymbol();
	if (arglist[0]->is_int() && arglist[1]->is_int() ){
		printTaxaTierTuples(arglist[0]->get_int(), arglist[1]->get_int() );
	}
    return result;
}

pqlsymbol * u_prototype(vector<pqlsymbol * > arglist){
	generate_random_bt();
	pqlsymbol * result = new pqlsymbol();
	//pqlsymbol * result = new pqlsymbol(distinguishing_bipart(arglist[0]->get_treeset(), arglist[1]->get_treeset()));
	return result;
}



pqlsymbol * u_distinguishing_taxa(vector<pqlsymbol * > arglist){
	pqlsymbol * result;
	if (arglist[0]->is_treeset() && arglist[1]->is_treeset() ){
		result = new pqlsymbol(distinguishing_taxa(arglist[0]->get_treeset(), arglist[1]->get_treeset()));
	}
	else{
		cout << "u_distinguishing_taxa expects two treesets"  << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}



void add_function(string functname, afptr funct, string doc ){
	::argFunctMap.insert(std::make_pair(functname, funct));
	//This is so tab-autocomplete works
	::functionKeys.push_back(functname);	
}

//adds the functions to the their maps. 
void init_the_functs()
{
	//::functionKeys.push_back("");
	
	//search
	add_function("get_trees_without_taxa", &u_get_trees_without_taxa, "Returns the trees that do not have the input taxa");
	add_function("gtwot", &u_get_trees_without_taxa, "Returns the trees that do not have the input taxa");

	add_function("clade_size_search", &u_clade_size_search, "Returns trees with clade of given size and taxa");
	add_function("smallest_clade", &u_smallest_clade, "Returns trees with the smallest clade of given taxa");
	add_function("get_trees_with_taxa", &u_get_trees_with_taxa, "Returns the trees that have the input taxa");
	add_function("gtwt", &u_get_trees_with_taxa, "Returns the trees that have the input taxa");
	add_function("similarity_search", &u_similarity_search, "Returns the trees that are the most similar to the given tree in that they have the most shared bipartitions.");

	//add_function("get_trees_by_taxa", &u_get_trees_by_taxa, " ");
	add_function("get_trees_by_subtree", &u_get_trees_by_subtree, "Returns the trees that contain the input subtree");

	add_function("subtree_search", &u_get_trees_by_subtree, "Returns the trees that contain the input subtree");
	//add_function("search_by_relationship", &u_search_hashtable_strict, "Returns the trees that have the input bipartition");
	add_function("structural_search", &u_new_search_hashtable_strict, "Returns the trees that have the input bipartition");
	add_function("ss", &u_new_search_hashtable_strict, "Returns the trees that have the input bipartition");
	
	
	//add_function("ss", &u_search_hashtable_strict, "Returns the trees that have the input bipartition");
	
	add_function("timedsearch", &u_search_hashtable_strict_and_timed, "for testing only");


	//interface with outside programs
	add_function("phlash", &u_phlash, "Calls the program Phlash");
	add_function("qq", &u_quick_quartet, "Calls the program QuickQuartet");
 	add_function("HashCS", &u_HashCS, " ");
 		
	//set operations
	add_function("not", &u_not, "Returns the difference between all trees and the input treeset.");
	add_function("union", &u_union, "Returns the union of two treesets.");
	add_function("intersection", &u_intersection, "Returns the intersection of two treesets." );
	add_function("difference", &u_difference, "Returns the difference of two treesets.");

	//convert types
	add_function("set", &u_set, "Converts a list of ints into a treeset.");
	add_function("to_newick", &u_to_newick, "Returns the newick string for the tree matching the input index.");
 
	//analysis
	add_function("count", &u_count, "Returns the number of objects in the treeset or list");
	add_function("unique", &u_unique, "Returns the a subset of trees from a given treeset each with a unique topology.");
	add_function("unique_biparts", &u_unique_biparts, "Returns the number of all unique bipartitions given a treeset");
	add_function("silhouette", &u_silhouette, "Returns the silhouette distance between given clusters of trees");
	add_function("agglo_clust", &u_agglo_clust, "Returns the agglomerative clustering of the given input set of trees.");
	add_function("duplicates", &u_duplicates, "Returns the set of trees with are topologically equal to the input tree.");
		
		//consensus
		add_function("consensus", &u_consen, "Returns the newick string for the consensus tree for the input treeset.");
		add_function("consensus_strict", &u_strict_consen, "Returns the newick string for the strict consensus tree for the input treeset.");
		add_function("consensus_majority", &u_majority_consen, "Returns the newick string for the majority consensus tree for the input treeset.");
		add_function("consensus_greedy", &u_greedy_consen, " ");	
		add_function("consensus_least_conflict", &u_least_conflict, " ");
		//should move to a single function that calls consensus_reso_rate, or reso rate depending on input.
		add_function("consensus_reso_rate", &u_consensus_reso_rate, "Returns the consensus resolution rate for a set of trees and a given consensus strictness.");
		add_function("crr", &u_consensus_reso_rate, "Returns the consensus resolution rate for a  set of trees and a given consensus strictness.");
		add_function("reso_rate", &u_reso_rate, "Returns the resolution rate for a given tree.");

	
	
	//utilities
	add_function("help", &u_help, " ");

	//visualization
	add_function("show", &u_show, "Displays images of the specified tree or trees (Int, IntVect, or Treeset). Takes an optional String mode argument ('ortho' or 'radial') to display an SVG image. Default mode is text.");
	add_function("show_newick", &u_show_newick, "Displays images of the specified Newick strings (String or StringVect). Takes an optional String mode argument ('ortho' or 'radial') to display an SVG image. Default mode is text.");
	add_function("export", &u_export, "Exports images of the specified tree or trees (Int, IntVect, or Treeset) to the specified folder path (String). Takes an optional String mode argument ('ortho' or 'radial') to export an SVG image. Default mode is text.");
	add_function("export_newick", &u_export_newick, "Exports images of the specified Newick strings (String or StringVect) to the specified folder path (String). Takes an optional String mode argument ('ortho' or 'radial') to display an SVG image. Default mode is text.");

	//classification
	add_function("classification", &u_classification, "Prints the classification data for a given taxon or vector of taxa.");
	add_function("write_classifications", &u_write_classifications, "Saves the classification data currently in memory to the specified filename.");
	add_function("edit_classification", &u_edit_classification, "Takes (taxa, rank, info). Edits the classification data for the specified taxon or vector of taxa by overwriting the taxon's existing info in the specified rank with the given info.");
	add_function("taxa_in_group", &u_taxa_in_group, "Returns the taxa belonging to the specified taxonomic group.");
	add_function("show_group_in_tree", &u_show_group_in_tree, "Displays an SVG image of the given taxonomic group within the given tree(s). Highlights taxa in the clade that conflict with the expected grouping.");
	add_function("show_group", &u_show_group, "Displays an SVG image of the given taxonomic group isolated from its parent tree(s). Highlights taxa in the clade that conflict with the expected grouping.");
	add_function("show_level", &u_show_level, "If possible, displays an image of the given tree(s) abstracted to the level of the given taxonomic rank.");
	add_function("show_only", &u_show_only, "Computes a version of the given tree(s) including only taxa from the given group(s)");

	//traits
	add_function("load_trait_matrix", &u_load_trait_matrix, "Loads trait information from a NEXUS character matrix.");
	add_function("test_trait_correlation", &u_test_trait_correlation, "Runs algorithm to test correlated evolution of specified traits");

	//modify treeset
	//add_function("taxa_filter", &u_taxa_filter, "Creates a new tree or trees by removing all but the specified taxa from the specified tree(s). Original trees are kept intact.");
	//add_function("group_filter", &u_group_filter, "Creates a new tree or trees by removing all taxa except the ones in the specified taxonomic groups from the specified tree(s). Original trees are kept intact.");
	//add_function("delete_tree", &u_delete_tree, "Deletes the specified tree(s) from the working data set.");
	add_function("write_trz", &u_write_trz, "Writes specified Tree/TreeVect/Treeset to a .trz file with specified filename.");

	//Developer Functions
	add_function("debug", &u_debug, " ");
	add_function("print_taxa_in_trees", &u_print_taxa_in_trees, " ");
	add_function("print_biparttable", &u_print_biparttable, " ");
 	add_function("set_hetero", &u_set_hetero, " ");
 	add_function("rsearch", &u_random_search, " ");
 	add_function("rsearch2", &u_random_search2, " ");
 	add_function("tier", &u_tier, " ");
 	add_function("tttier", &u_tttier, " ");
 	add_function("rate_tiers", &u_rate_tiers, " ");
 	add_function("print_label_map", &u_print_label_map, " ");
 	add_function("get_taxa_in_tree", &u_get_taxa_in_tree, "Returns the taxa in a given tree. Takes a tree index");
 	add_function("distinguishing_taxa", &u_distinguishing_taxa, "Takes two Treesets. Returns the taxa that are in of one treeset and none of the other.");
 	add_function("proto", &u_prototype, "This is a protype function only use for the function in current developemnt.");
 	
 	cout << "Functions loaded into map" << endl;
 	
 	
 	
	
}

