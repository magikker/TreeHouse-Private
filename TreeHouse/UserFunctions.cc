#include "UserFunctions.h"
#include <boost/variant.hpp>
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


//needs type checking
pqlsymbol * u_get_trees_by_subtree(vector< pqlsymbol * > arglist) {  
	return new pqlsymbol(get_trees_by_subtree(arglist[0]->get_string() ) );
}

pqlsymbol * u_rf_dist_clade(vector< pqlsymbol * > arglist) { 
	if (arglist.size() == 2 ){
		if (arglist[0]->is_string() && arglist[1]->is_string()  ) {
			return new pqlsymbol(CRFStringDistance(arglist[0]->get_string(), arglist[1]->get_string() ) );
		}
		else if (arglist[0]->is_int() && arglist[1]->is_int()  ) {
			return new pqlsymbol(CRFDistance(arglist[0]->get_int(), arglist[1]->get_int()  ) );
		}
	}
	return new pqlsymbol(ERROR, "Function requires either two ints or two newick strings");
}

pqlsymbol * u_rf_dist_bipart(vector< pqlsymbol * > arglist) {  
	if (arglist.size() == 2 ){
		if (arglist[0]->is_string() && arglist[1]->is_string()  ) {
			return new pqlsymbol(BipartRFStringDistance(arglist[0]->get_string(), arglist[1]->get_string() ) );
		}
		else if (arglist[0]->is_int() && arglist[1]->is_int()  ) {
			return new pqlsymbol(BipartRFDistance(arglist[0]->get_int(), arglist[1]->get_int() ) );
		}
	}
	return new pqlsymbol(ERROR, "Function requires either two ints or two newick strings");
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
			result = new pqlsymbol(get_trees_by_taxa(arglist[0]->get_string_vect(), arglist[1]->get_string_vect()) );
		}
	}
	else {	
		cout << "get_trees_by_taxa expects 2 string vector arguments." << "Found " << get_arg_types(std::move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

pqlsymbol * u_shared_quartets_strict(vector< pqlsymbol * > arglist){  
	pqlsymbol * result;
	if(arglist.size()!=1 || !arglist[0]->is_treeset()){
		cout << "Error: this function expects one input- a treeSet. It was given " << get_arg_types(std::move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	else{
		shared_quartets_strict(arglist[0]->get_treeset());
		result = new pqlsymbol();
	}
	return result;
}

pqlsymbol * u_shared_quartets_majority(vector< pqlsymbol * > arglist){  
	pqlsymbol * result;
	if(arglist.size()!=1 || !arglist[0]->is_treeset()){
		cout << "Error: this function expects one input- a treeSet. It was given " << get_arg_types(std::move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	else{
		shared_quartets_majority(arglist[0]->get_treeset());
		result = new pqlsymbol();
	}
	return result;
}

pqlsymbol * u_print_quartets_from_tree(vector< pqlsymbol * > arglist){  
	pqlsymbol * result;
	if(arglist.size()!=1 || !arglist[0]->is_int()){
		cout << "Error: this function expects one input- an int. It was given " << get_arg_types(std::move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	else{
		printSet(generateQuartetsFromTree(arglist[0]->get_int()));
		result = new pqlsymbol();
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

/*
   pqlsymbol * u_get_subset_trees(vector< pqlsymbol * > arglist) {
   pqlsymbol * result;
   if(arglist.size()!=1 || !arglist[0]->is_int()){
   cout << "Error: get_subset_trees expects ONE argument of type INT!" << " Found " << get_arg_types(move(arglist)) << endl;
   result = new pqlsymbol(ERROR,"Type Error");
   }
   else{
   result = new pqlsymbol(get_subset_trees(arglist[0]->get_int()));
   }
   return result;
   }

   pqlsymbol * u_get_superset_trees(vector< pqlsymbol * > arglist) {
   pqlsymbol * result;
   if(arglist.size()!=1 || !arglist[0]->is_int()){
   cout << "Error: get_superset_trees expects ONE argument of type INT!" << " Found " << get_arg_types(move(arglist)) << endl;
   result = new pqlsymbol(ERROR,"Type Error");
   }
   else{
   result = new pqlsymbol(get_superset_trees(arglist[0]->get_int()));
   }
   return result;
   }
   */
pqlsymbol * u_average_ancestral_distance(vector< pqlsymbol * > arglist) {
	pqlsymbol * result;
	bool error = false;
	if(arglist.size()<2 || arglist.size()>3){
		error = true;
	}
	else if(arglist.size()==2){
		if(arglist[0]->is_int() && arglist[1]->is_int()){
			result = new pqlsymbol(average_ancestral_distance(arglist[0]->get_int(), arglist[1]->get_int()));
		}
		else{
			error = true;
		}
	}
	else if(arglist.size()==3){
		if(arglist[0]->is_int() && arglist[1]->is_int() && arglist[2]->is_treeset()){
			result = new pqlsymbol(average_ancestral_distance(arglist[0]->get_int(), arglist[1]->get_int(), arglist[2]->get_treeset()));
		}
		else{
			error = true;
		}
	}	
	if(error){
		cout << "Error: average_ancestral_distance expects two INTS and an OPTIONAL treeset (leave empty for all trees)!" << " Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR,"Type Error");
	}

	return result;
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
			result = new pqlsymbol(clade_size_search(temp, arglist[1]->get_int() ) );
		}
	}
	else if (arglist[0]->is_vect() && arglist[0]->is_string() && arglist[1]->is_int()) {
		vector<string> missednames = ::biparttable.lm.catchDeclaredTaxa(arglist[0]->get_string_vect());
		if(missednames.size() > 0){
			result = new pqlsymbol(ERROR, missedTaxaErrorMessage(missednames));
		}
		else{
			result = new pqlsymbol(clade_size_search(arglist[0]->get_string_vect(), arglist[1]->get_int() ) );
		}

	}
	else if (arglist[0]->is_vect() && arglist[0]->is_int() && arglist[1]->is_int()) {
		result = new pqlsymbol(clade_size_search(arglist[0]->get_int_vect(), arglist[1]->get_int() ) );
	}
	else {
		cout << "clade_size_search expects a StringVec/Intvec and an Integer" << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

pqlsymbol * u_smallest_clade_search(vector< pqlsymbol * > arglist) {
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
			result = new pqlsymbol(smallest_clade_search(temp) );
		}
	}
	else if (arglist[0]->is_vect() && arglist[0]->is_string()) {
		vector<string> missednames = ::biparttable.lm.catchDeclaredTaxa(arglist[0]->get_string_vect());
		if(missednames.size() > 0){
			result = new pqlsymbol(ERROR, missedTaxaErrorMessage(missednames));
		}
		else{
			result = new pqlsymbol(smallest_clade_search(arglist[0]->get_string_vect()) );
		}

	}
	else if (arglist[0]->is_vect() && arglist[0]->is_int()) {
		result = new pqlsymbol(smallest_clade_search(arglist[0]->get_int_vect()) );
	}
	else {
		cout << "smallest_clade expects a StringVec/Intvec. " << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

pqlsymbol * u_similarity_search(vector< pqlsymbol * > arglist) {
	pqlsymbol * result;
	string stree;

	stree = arglist[0]->get_string();
	result = new pqlsymbol(similarity_search(stree));

	return result;
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
			result = new pqlsymbol(get_trees_with_taxa(temp) );
		}
	}
	else if (arglist[0]->is_vect() && arglist[0]->is_string() ) {
		vector<string> missednames = ::biparttable.lm.catchDeclaredTaxa(arglist[0]->get_string_vect());

		if(missednames.size() > 0){
			result = new pqlsymbol(ERROR, missedTaxaErrorMessage(missednames));
		}
		else{
			result = new pqlsymbol(get_trees_with_taxa(arglist[0]->get_string_vect() ) );
		}

	}
	else if (arglist[0]->is_vect() && arglist[0]->is_int() ) {
		result = new pqlsymbol(get_trees_with_taxa(arglist[0]->get_int_vect() ) );
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
		result = new pqlsymbol(get_trees_without_taxa(arglist[0]->get_string_vect() ) );
	}
	else if (arglist[0]->is_int() && arglist[0]->is_vect()){
		result = new pqlsymbol(get_trees_without_taxa(arglist[0]->get_int_vect() ) );
	}
	else{
		cout << "get_trees_without_taxa expects either 1 StringVect or 1 IntVect argument. " << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	return result;
}

/*
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
*/

pqlsymbol * u_to_newick(vector<pqlsymbol *> arglist){

	to_newick(arglist[0]->get_int());

	return new pqlsymbol(to_newick(arglist[0]->get_int()));
}

pqlsymbol * u_least_conflict(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	result = new pqlsymbol(least_conflict(arglist[0]->get_treeset()) );

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
	else if(arglist.size()==1 && arglist[0]->is_int()){	
		int tree = arglist[0]->get_int();
		stree = to_newick(tree);
		//append a semicolon at the end of stree because that is the format it needs to be in
		stree.push_back(';');
	}	
	else{
		cout << "reso_rate expects a single string representing a tree " << "Found " << get_arg_types(move(arglist)) << endl;
		return new pqlsymbol(ERROR, "Type error");
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
result = new pqlsymbol(search_hashtable_strict_old( left, right, strict ) );

}
else{
cout << "search_hashtable_strict expects either 2 string vectors or 2 string vectors and an int. " << "Found " << get_arg_types(arglist) << endl;
result = new pqlsymbol(ERROR, "Type Error");
}

return result;
}
*/

pqlsymbol * u_rsk(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;

	int left = 2;
	int right = 2;
	int iterations = 1;

	if (arglist.size() >= 2){
		left = arglist[0]->get_int();
		right = arglist[1]->get_int();
	}
	if (arglist.size() >= 3){
		iterations = arglist[2]->get_int();
	}

	result = new pqlsymbol(random_search_ktet( left, right, iterations ) );

	return result;
}

pqlsymbol * u_ssk(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;

	int left = 2;
	int right = 2;
	int iterations = 1;

	if (arglist.size() >= 2){
		left = arglist[0]->get_int();
		right = arglist[1]->get_int();
	}
	if (arglist.size() >= 3){
		iterations = arglist[2]->get_int();
	}

	result = new pqlsymbol(sucess_search_ktet( left, right, iterations ) );

	return result;
}

pqlsymbol * u_fsk(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;

	int left = 2;
	int right = 2;
	int iterations = 1;

	if (arglist.size() >= 2){
		left = arglist[0]->get_int();
		right = arglist[1]->get_int();
	}
	if (arglist.size() >= 3){
		iterations = arglist[2]->get_int();
	}

	result = new pqlsymbol(fail_search_ktet( left, right, iterations ) );

	return result;
}

pqlsymbol * u_rssc(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;

	int left = 2;
	int iterations = 1;

	if (arglist.size() >= 1){
		left = arglist[0]->get_int();
	}
	if (arglist.size() >= 2){
		iterations = arglist[1]->get_int();
	}

	result = new pqlsymbol(random_search_subclade( left, iterations ) );

	return result;
}

pqlsymbol * u_sssc(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;

	int left = 2;
	int iterations = 1;

	if (arglist.size() >= 1){
		left = arglist[0]->get_int();
	}
	if (arglist.size() >= 2){
		iterations = arglist[1]->get_int();
	}

	result = new pqlsymbol(sucess_search_subclade( left, iterations ) );

	return result;
}

pqlsymbol * u_fssc(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;

	int left = 2;
	int iterations = 1;

	if (arglist.size() >= 1){
		left = arglist[0]->get_int();
	}
	if (arglist.size() >= 2){
		iterations = arglist[1]->get_int();
	}
	result = new pqlsymbol(fail_search_subclade( left, iterations ) );

	return result;
}

pqlsymbol * u_ssc(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;

	int left = 2;
	int iterations = 1;

	if (arglist.size() >= 1){
		left = arglist[0]->get_int();
	}
	if (arglist.size() >= 2){
		iterations = arglist[1]->get_int();
	}

	result = new pqlsymbol(sucess_search_clade( left, iterations ) );

	return result;
}

pqlsymbol * u_fsc(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;

	int left = 2;
	int iterations = 1;

	if (arglist.size() >= 1){
		left = arglist[0]->get_int();
	}
	if (arglist.size() >= 2){
		iterations = arglist[1]->get_int();
	}
	result = new pqlsymbol(fail_search_clade( left, iterations ) );

	return result;
}



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
result = new pqlsymbol(search_hashtable_strict_and_timed( left, right, strict ) );
}
else{
cout << "consensus expects either 2 string vectors or 2 string vectors and an int. " << "Found " << get_arg_types(arglist) << endl;
result = new pqlsymbol(ERROR, "Type Error");
}

return result;
}
*/
pqlsymbol * u_search_ktet(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;

	//type check and catch errors and handle any method overloading. 

	vector<int> left;
	vector<int> right;

	bool LeftIsGood = false;
	bool RightIsGood = false;
	int strict = 0;

	if (arglist.size() == 2){
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

	if(LeftIsGood && RightIsGood){
		result = new pqlsymbol(search_ktet( left, right ) );
	}
	else{
		cout << "consensus expects either 2 string vectors. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}

	return result;
}

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
		result = new pqlsymbol(search_hashtable_strict_and_timed( left, right, strict ) );
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
		result = new pqlsymbol(search_hashtable_strict( left, right, strict ) );
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

/*int help( ) 
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
  */

pqlsymbol * u_help(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;	
	if(arglist.size()!=1 || !arglist[0]->is_string()){
		cout << "help expects one string argument. " << "Found " << get_arg_types(arglist) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
	else{
		help(arglist[0]->get_string());
		result = new pqlsymbol();
	}


	return result;
}

pqlsymbol * u_unique(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	result = new pqlsymbol(unique(arglist[0]->get_treeset()) );

	return result;
}

//unique_biparts(set< unsigned int > treesin)
pqlsymbol * u_unique_biparts(vector<pqlsymbol *> arglist)
{
	pqlsymbol * result;

	result = new pqlsymbol(unique_biparts(arglist[0]->get_treeset()));

	return result;
}

pqlsymbol * u_silhouette(vector<pqlsymbol * > arglist) {
	pqlsymbol * result = new pqlsymbol();

	result = new pqlsymbol(silhouette(arglist[0]->get_treeset_vect(), arglist[1]->get_string()));

	return result;

}

pqlsymbol * u_rand_index(vector <pqlsymbol *> arglist){
	pqlsymbol * result = new pqlsymbol();

	result = new pqlsymbol(rand_index(arglist[0]->get_treeset_vect(), arglist[1]->get_treeset_vect()));

	return result;
}

pqlsymbol * u_adjusted_rand_index(vector <pqlsymbol *> arglist){
	pqlsymbol * result = new pqlsymbol();

	result = new pqlsymbol(adjusted_rand_index(arglist[0]->get_treeset_vect(), arglist[1]->get_treeset_vect()));

	return result;
}

pqlsymbol * u_agglo_clust(vector<pqlsymbol * > arglist){
	pqlsymbol * result = new pqlsymbol();

	result = new pqlsymbol(agglo_clust(arglist[0]->get_treeset(), arglist[1]->get_int(), arglist[2]->get_string()));

	return result;

}

pqlsymbol * u_kmeans_clust(vector<pqlsymbol *> arglist){
	pqlsymbol * result = new pqlsymbol();

	result = new pqlsymbol(kmeans_clust(arglist[0]->get_treeset(),arglist[1]->get_int(), arglist[2]->get_string()));

	return result;
}

pqlsymbol * u_dbscan_clust(vector <pqlsymbol *> arglist){
	pqlsymbol * result = new pqlsymbol();

	result = new pqlsymbol(dbscan_clust(arglist[0]->get_treeset(),arglist[1]->get_int(),arglist[2]->get_int(),arglist[3]->get_string()));

	return result;
}

pqlsymbol * u_burnin_clust(vector <pqlsymbol *> arglist){
	pqlsymbol * result = new pqlsymbol();

	result = new pqlsymbol(burnin_clust(arglist[0]->get_treeset(), arglist[1]->get_string()));
	return result;
}

pqlsymbol * u_search_clade(vector <pqlsymbol *> arglist){
	pqlsymbol * result = new pqlsymbol();

	//bool strict = false;
	//if(arglist.size() == 1){
		 //string st = arglist[1]->get_string();
		 //if (boost::iequals(st, "strict")){
			// strict = true;
		// }
	//}

	result = new pqlsymbol(search_clade(arglist[0]->get_string_vect()));
	return result;
}

// takes an int. Returns the ints which share the same topology. 
pqlsymbol * u_duplicates(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	result = new pqlsymbol(biparttable.duplicates(arglist[0]->get_int() ) );


	return result;
}

pqlsymbol * u_sample_trees(vector<pqlsymbol *> arglist)
{
	pqlsymbol * result;

	result = new pqlsymbol(sample_trees(arglist[0]->get_treeset(), arglist[1]->get_int()));

	return result;
}


pqlsymbol * u_count(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;

	result = new pqlsymbol(arglist[0]->get_size());


	return result;
}

pqlsymbol * u_union(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;


	set<unsigned int> s1 = arglist[0]->get_treeset();
	set<unsigned int> s2 = arglist[1]->get_treeset();
	set<unsigned int> sunion; 

	std::set_union( s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter( sunion, sunion.begin() ) );

	result = new pqlsymbol(sunion);


	return result;
}


pqlsymbol * u_difference(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;


	set<unsigned int> s1 = arglist[0]->get_treeset();
	set<unsigned int> s2 = arglist[1]->get_treeset();
	set<unsigned int> sdiff; 

	std::set_difference( s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter( sdiff, sdiff.begin() ) );

	result = new pqlsymbol(sdiff);


	return result;
}


pqlsymbol * u_intersection(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;


	set<unsigned int> s1 = arglist[0]->get_treeset();
	set<unsigned int> s2 = arglist[1]->get_treeset();
	set<unsigned int> sinter; 

	std::set_intersection( s1.begin(), s1.end(), s2.begin(), s2.end(), std::inserter( sinter, sinter.begin() ) );

	result = new pqlsymbol(sinter);


	return result;
}

pqlsymbol * u_not(vector<pqlsymbol * > arglist) 
{  
	pqlsymbol * result;


	set<unsigned int> s1 = arglist[0]->get_treeset();
	set<unsigned int> sdiff; 

	std::set_difference( all_trees.begin(), all_trees.end(), s1.begin(), s1.end(),  std::inserter( sdiff, sdiff.begin() ) );

	result = new pqlsymbol(sdiff);


	return result;
}

pqlsymbol * u_set(vector<pqlsymbol * > arglist) {  
	pqlsymbol * result;


	vector<int> tempvect = arglist[0]->get_int_vect(); 
	set<unsigned int> s1(tempvect.begin(), tempvect.end());

	result = new pqlsymbol(s1);
	return result;
}

//User function for display_clusters, opens gnuplot with a graph of the trees color coded by cluster
pqlsymbol * u_display_clusters(vector<pqlsymbol *> arglist){
	pqlsymbol * result = new pqlsymbol();

	display_clusters(arglist[0]->get_string(), arglist[1]->get_string(), arglist[2]->get_treeset_vect());
	return result;
}

pqlsymbol * u_display_heatmap(vector<pqlsymbol *> arglist){
	pqlsymbol * result = new pqlsymbol();
	display_heatmap(arglist[0]->get_treeset(), arglist[1]->get_string(), arglist[2]->get_string());

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

pqlsymbol * u_conflicting_quartet_distance(vector<pqlsymbol * > arglist) {
	pqlsymbol * result = new pqlsymbol();
	if (arglist.size() == 2) {
		if(arglist[0]->is_int() && arglist[1]->is_int()){
			return new pqlsymbol(conflictingQuartetDistance(arglist[0]->get_int(), arglist[1]->get_int()));
			}	
		else{
			cout << "conflicting_quartet_distance expect two ints (tree indices)." << endl;
			result = new pqlsymbol(ERROR, "Type Error");
			}
	}
	 else {
		cout << "conflicting_quartet_distance expect two ints (tree indices)." << endl;
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
			::biparttable.hetero = true;
			cout << "HETERO set to true!\n";
			result = new pqlsymbol("HETERO set to true");
		}
		else{
			::biparttable.hetero = false;
			cout << "HETERO set to false!\n";
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

	result = new pqlsymbol(biparttable.get_taxa_in_tree(arglist[0]->get_int()));


	return result;
}

pqlsymbol * u_print_biparttable(vector<pqlsymbol * > arglist) {
	pqlsymbol * result = new pqlsymbol();
	::biparttable.print_biparttable();
	return result;
}

pqlsymbol * u_print_cladetable(vector<pqlsymbol * > arglist) {
	pqlsymbol * result = new pqlsymbol();
	::biparttable.print_clade_table();
	return result;
}

pqlsymbol * u_print_bipartitions(vector<pqlsymbol * > arglist) {
	pqlsymbol * result = new pqlsymbol();
	::biparttable.print_bipartitions(arglist[0]->get_int());
	return result;
}

pqlsymbol * u_print_clades(vector<pqlsymbol * > arglist) {
	pqlsymbol * result = new pqlsymbol();
	::biparttable.print_clades(arglist[0]->get_int());
	return result;
}

pqlsymbol * u_print_inverted_index(vector<pqlsymbol * > arglist) {
	pqlsymbol * result = new pqlsymbol();
	::biparttable.print_inverted_index();
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
	cout << "entering the prototype function" << endl;
	pqlsymbol * result = new pqlsymbol();
	::biparttable.print_CladeMap();
	//dTree(arglist[0]->get_treeset_vect());
	return result;
}

/*
pqlsymbol * u_prototype(vector<pqlsymbol * > arglist){
	cout << "entering the prototype function" << endl;
	pqlsymbol * result = new pqlsymbol();
	psupport(arglist[0]->get_treeset_vect());
	return result;
}
*/


pqlsymbol * u_homogenize(vector<pqlsymbol * > arglist){
	pqlsymbol * result = new pqlsymbol();
	vector<unsigned int> mask;
	mask = ::biparttable.homog_taxa();
	::biparttable.apply_taxa_mask(mask);
	return result;
}


pqlsymbol * u_distinguishing_taxa(vector<pqlsymbol * > arglist){
	pqlsymbol * result;

	result = new pqlsymbol(distinguishing_taxa(arglist[0]->get_treeset(), arglist[1]->get_treeset()));


	return result;
}

pqlsymbol * u_distance_between_taxa(vector<pqlsymbol * > arglist){
	pqlsymbol * result;

	result = new pqlsymbol(distance_between_taxa(arglist[0]->get_int(), arglist[1]->get_int(), arglist[2]->get_int()));

}

pqlsymbol * u_distance_to_common_ancestor(vector<pqlsymbol * > arglist){
	pqlsymbol * result;


	result = new pqlsymbol(distance_to_common_ancestor(arglist[0]->get_int(), arglist[1]->get_int(), arglist[2]->get_int()));

	return result;
}

pqlsymbol * u_distance_to_root(vector<pqlsymbol * > arglist){
	pqlsymbol * result;

	if (arglist.size()==2)

		result = new pqlsymbol(distance_to_root(arglist[0]->get_int(), arglist[1]->get_int()));


	return result;
}


pqlsymbol * u_print_conflicting_quartets(vector<pqlsymbol * > arglist){
  pqlsymbol * result;

  if (arglist.size()==2)
  {
	  if(arglist[0]->is_int() && arglist[1]->is_int()){
	     set<quartet> conflictingQuartets = generateConflictingQuartets(arglist[0]->get_int(), arglist[1]->get_int());
   	     printSet(conflictingQuartets);
	     result = new pqlsymbol();
	  }
  } 
  else {
	cout << "print_conflicting_quartets expects two INTs, i.e. bipartition indices. "  << "Found " << get_arg_types(move(arglist)) << endl;
	result = new pqlsymbol(ERROR, "Type Error");
  	}
}

pqlsymbol * u_average_depth(vector<pqlsymbol * > arglist){
	pqlsymbol * result;


	result = new pqlsymbol(average_depth(arglist[0]->get_int()));


	return result;
}
/*
pqlsymbol * u_print_conflicting_quartets(vector<pqlsymbol * > arglist){
	pqlsymbol * result;

	if (arglist.size()==2)
	{
		if(arglist[0]->is_int() && arglist[1]->is_int()){
			set<quartet> conflictingQuartets = generateConflictingQuartets(arglist[0]->get_int(), arglist[1]->get_int());
			printSet(conflictingQuartets);
			result = new pqlsymbol();
		}
	} 
	else {
		cout << "print_conflicting_quartets expects two INTs, i.e. bipartition indices. "  << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}
>>>>>>> a9281369f04d5a01ccb75c1a457ae2070e310c1e

	return result;
}
*/

pqlsymbol * u_get_num_quartets(vector<pqlsymbol * > arglist){
  pqlsymbol * result;

  if (arglist.size()==1)
  {
	  if(arglist[0]->is_int()){
   	     
	     result = new pqlsymbol(getNumQuartets(arglist[0]->get_int()));
	}
	else{
	cout << "get_num_quartets expects one argument of type INT.. "  << "Found " << get_arg_types(move(arglist)) << endl;
	result = new pqlsymbol(ERROR, "Type Error");
  	}
  } 
  else {
	cout << "get_num_quartets expects one argument of type INT.. "  << "Found " << get_arg_types(move(arglist)) << endl;
	result = new pqlsymbol(ERROR, "Type Error");
  	}

  return result;
}





pqlsymbol * u_num_conflicting_quartets(vector<pqlsymbol * > arglist){
	pqlsymbol * result;

	if (arglist.size()==2)
	{
		if(arglist[0]->is_int() && arglist[1]->is_int()){
			result = new pqlsymbol(numConflictingQuartets(arglist[0]->get_int(), arglist[1]->get_int()));
		}
	} 
	else {
		cout << "num_conflicting_quartets expects two INTs, i.e. bipartition indices. "  << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
	}

	return result;
}

pqlsymbol * u_num_shared_quartets(vector<pqlsymbol * > arglist){
  pqlsymbol * result;

  if (arglist.size()==2)
  {
	  if(arglist[0]->is_int() && arglist[1]->is_int()){
	     result = new pqlsymbol(getNumSameQuartets(arglist[0]->get_int(), arglist[1]->get_int()));
	  }
  } 
  else {
	cout << "num_conflicting_quartets expects two INTs, i.e. bipartition indices. "  << "Found " << get_arg_types(move(arglist)) << endl;
	result = new pqlsymbol(ERROR, "Type Error");
  	}

  return result;
}


pqlsymbol * u_num_overlapping_quartets(vector<pqlsymbol * > arglist){
  pqlsymbol * result;

  if (arglist.size()==2)
  {
	  if(arglist[0]->is_int() && arglist[1]->is_int()){
	     result = new pqlsymbol(numOverlappingQuartets(arglist[0]->get_int(), arglist[1]->get_int()));
	  }
  } 
  else {
	cout << "num_overlapping_quartets expects two INTs, i.e. bipartition indices. "  << "Found " << get_arg_types(move(arglist)) << endl;
	result = new pqlsymbol(ERROR, "Type Error");
  	}

  return result;
}

pqlsymbol * u_quartet_distance(vector<pqlsymbol * > arglist){
  pqlsymbol * result;

  if (arglist.size()==2)
  {
	  if(arglist[0]->is_int() && arglist[1]->is_int()){
 	     set<quartet> qDist = generateDifferentQuartetsFromTrees3(arglist[0]->get_int(), arglist[1]->get_int());
	     unsigned int size = qDist.size();
	     result = new pqlsymbol(size);
	  }
  } 
  else {
	cout << "quartet_distance expects two INTs, i.e. tree indices. "  << "Found " << get_arg_types(move(arglist)) << endl;
	result = new pqlsymbol(ERROR, "Type Error");
  	}

  return result;
}

pqlsymbol * u_print_quartet_distance(vector<pqlsymbol * > arglist){
  pqlsymbol * result;

  if (arglist.size()==2)
  {
	  if(arglist[0]->is_int() && arglist[1]->is_int()){
 	     set<quartet> qDist = generateDifferentQuartetsFromTrees3(arglist[0]->get_int(), arglist[1]->get_int());
	     printSet(qDist);
	     result = new pqlsymbol();
	  }
  } 
  else {
	cout << "print_quartet_distance expects two INTs, i.e. tree indices. "  << "Found " << get_arg_types(move(arglist)) << endl;
	result = new pqlsymbol(ERROR, "Type Error");
  	}

  return result;
}

pqlsymbol * u_print_shared_quartets(vector<pqlsymbol * > arglist){
  pqlsymbol * result;

  if (arglist.size()==2)
  {
	  if(arglist[0]->is_int() && arglist[1]->is_int()){
		printSameQuartets(arglist[0]->get_int(), arglist[1]->get_int());	     
		result = new pqlsymbol();
	  }
  } 
  else {
	cout << "num_conflicting_quartets expects two INTs, i.e. bipartition indices. "  << "Found " << get_arg_types(move(arglist)) << endl;
	result = new pqlsymbol(ERROR, "Type Error");
  	}

  return result;
}

pqlsymbol * u_is_bifurcating(vector<pqlsymbol * > arglist){
  pqlsymbol * result;

  if (arglist.size()==1)
  {
	  if(arglist[0]->is_int()){
		result = new pqlsymbol(isBifurcating(to_newick(arglist[0]->get_int())));
	  }
	else if(arglist[0]->is_string()){
		result = new pqlsymbol(isBifurcating(arglist[0]->get_string()));
		}
	else{
	cout << "is_bifurcating expects either an int or a string (representing a tree). "  << "Found " << get_arg_types(move(arglist)) << endl;
	result = new pqlsymbol(ERROR, "Type Error");
		}
  } 
  else {
	cout << "is_bifurcating expects either an int or a string (representing a tree). "  << "Found " << get_arg_types(move(arglist)) << endl;
	result = new pqlsymbol(ERROR, "Type Error");
  	}

  return result;
}
/*
pqlsymbol * u_average_depth(vector<pqlsymbol * > arglist){
  pqlsymbol * result;

  if (arglist.size()==1)
  {
	  if(arglist[0]->is_int()){
		result = new pqlsymbol(average_depth(arglist[0]->get_int()));
	  }
	else if(arglist[0]->is_treeset()){
		result = new pqlsymbol(average_depth(arglist[0]->get_treeset()));
		}
	else{
		cout << "average_depth expects an int (i.e. tree index) or a treeset."  << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");
		}
		
  } 
  else {
	cout << "average_depth expects an int (i.e. tree index) or a treeset."  << "Found " << get_arg_types(move(arglist)) << endl;
	result = new pqlsymbol(ERROR, "Type Error");
  	}

  return result;
}
*/

pqlsymbol * u_average_distance_between_taxa(vector<pqlsymbol * > arglist){
	pqlsymbol * result;

	if (arglist.size()==2)
	{
		if(arglist[0]->is_int() && arglist[1]->is_int()){
			result = new pqlsymbol(average_distance_between_taxa(arglist[0]->get_int(), arglist[1]->get_int()));
		}
	} 
	else if (arglist.size()==3) {
		if(arglist[0]->is_int() && arglist[1]->is_int() && arglist[2]->is_treeset()){
			result = new pqlsymbol(average_distance_between_taxa(arglist[0]->get_int(), arglist[1]->get_int(), arglist[2]->get_treeset()));
		}
		else{ cout << "average_distance_between_taxa expects two ints, i.e. two taxa, and an optional treeset. "  << "Found " << get_arg_types(move(arglist)) << endl;
			result = new pqlsymbol(ERROR, "Type Error");
		}
	}
	else{
		cout << "average_distance_between_taxa expects two ints, i.e. two taxa, and an optional treeset. "  << "Found " << get_arg_types(move(arglist)) << endl;
		result = new pqlsymbol(ERROR, "Type Error");

	}

	return result;
}

pqlsymbol * u_expected_depth(vector<pqlsymbol * > arglist){
	pqlsymbol * result;
	
	if (arglist.size()==1)
	{
		if(arglist[0]->is_int() && arglist[1]->is_int()){
			result = new pqlsymbol(expected_average_depth(arglist[0]->get_int()));
		}
		else{
			cout << "expected_average_depth expects an int. "  << "Found " << get_arg_types(move(arglist)) << endl;
  			result = new pqlsymbol(ERROR, "Type Error");
			}
	} 
	
  	else{
  		cout << "expected_average_depth expects an int. "  << "Found " << get_arg_types(move(arglist)) << endl;
  		result = new pqlsymbol(ERROR, "Type Error");
  		
	}
	
	return result;
}

pqlsymbol * u_calculate_C(vector<pqlsymbol * > arglist){
	pqlsymbol * result;
	
	if (arglist.size()==1)
	{
		if(arglist[0]->is_int()){
			result = new pqlsymbol(calculate_C(arglist[0]->get_int()));
		}
		else if(arglist[0]->is_string()){
			result = new pqlsymbol(calculate_C(arglist[0]->get_string()));
			}
		else{
  			cout << "calculate_C expects either an int or a Newick string. "  << "Found " << get_arg_types(move(arglist)) << endl;
  			result = new pqlsymbol(ERROR, "Type Error");
			}
	} 
	
  	else{
  		cout << "calculate_C expects either an int or a Newick string. "  << "Found " << get_arg_types(move(arglist)) << endl;
  		result = new pqlsymbol(ERROR, "Type Error");
  		
	}
	
	return result;
}

pqlsymbol * u_edit_distance_greedy(vector<pqlsymbol * > arglist){
	pqlsymbol * result;
	
	if (arglist.size()==2)
	{
		if(arglist[0]->is_int() && arglist[1]->is_int()){
			result = new pqlsymbol(edit_distance_greedy(arglist[0]->get_int(), arglist[1]->get_int()));
		}
		else{
  			cout << "edit_distance_greedy expects two trees. "  << "Found " << get_arg_types(move(arglist)) << endl;
  			result = new pqlsymbol(ERROR, "Type Error");
			}
	} 
	
  	else{
  		cout << "edit_distance_greedy expects two trees. "  << "Found " << get_arg_types(move(arglist)) << endl;
  		result = new pqlsymbol(ERROR, "Type Error");
 
	}
	return result;
}

pqlsymbol * u_edit_distance_total(vector<pqlsymbol * > arglist){
	pqlsymbol * result;
	
	if (arglist.size()==2)
	{
		if(arglist[0]->is_int() && arglist[1]->is_int()){
			result = new pqlsymbol(edit_distance_total(arglist[0]->get_int(), arglist[1]->get_int()));
		}
		else{
  			cout << "edit_distance_total expects two trees. "  << "Found " << get_arg_types(move(arglist)) << endl;
  			result = new pqlsymbol(ERROR, "Type Error");
			}
	} 
	
  	else{
  		cout << "edit_distance_total expects two trees. "  << "Found " << get_arg_types(move(arglist)) << endl;
  		result = new pqlsymbol(ERROR, "Type Error");
 
	}
	return result;
}

pqlsymbol * u_edit_distance_average(vector<pqlsymbol * > arglist){
	pqlsymbol * result;
	
	if (arglist.size()==2)
	{
		if(arglist[0]->is_int() && arglist[1]->is_int()){
			result = new pqlsymbol(edit_distance_average(arglist[0]->get_int(), arglist[1]->get_int()));
		}
		else{
  			cout << "edit_distance_average expects two trees. "  << "Found " << get_arg_types(move(arglist)) << endl;
  			result = new pqlsymbol(ERROR, "Type Error");
			}
	} 
	
  	else{
  		cout << "edit_distance_average expects two trees. "  << "Found " << get_arg_types(move(arglist)) << endl;
  		result = new pqlsymbol(ERROR, "Type Error");
 
	}
	return result;
}

pqlsymbol * u_edit_distance_minimum(vector<pqlsymbol * > arglist){
	pqlsymbol * result;
	
	if (arglist.size()==2)
	{
		if(arglist[0]->is_int() && arglist[1]->is_int()){
			result = new pqlsymbol(edit_distance_minimum(arglist[0]->get_int(), arglist[1]->get_int()));
		}
		else{
  			cout << "edit_distance_minimum expects two trees. "  << "Found " << get_arg_types(move(arglist)) << endl;
  			result = new pqlsymbol(ERROR, "Type Error");
			}
	} 
	
  	else{
  		cout << "edit_distance_minimum expects two trees. "  << "Found " << get_arg_types(move(arglist)) << endl;
  		result = new pqlsymbol(ERROR, "Type Error");
 
	}
	return result;
}


void add_function(string functname, afptr funct, string doc ){
	::argFunctMap.insert(std::make_pair(functname, funct));
	//This is so tab-autocomplete works
	::functionKeys.push_back(functname);
	::helpRef.insert(std::make_pair(functname, doc));

}

//A function to directly check if the given type matches the expected type
//returns a string, changes a boolean variable given
string type_check(dataType expected, pqlsymbol * given, bool &error){
	string ret_string;
	switch(expected){

		case TYPE_ATOM: //Atom
			if(!given->is_atom()){
				ret_string = "Expects an atom as input ";
				error = true;
			}
			break;
		case TYPE_INT: //int
			if(!given->is_int()){
				ret_string =  "Expects an integer as input ";
				error = true;
			}
			break;
		case TYPE_STRING: //String
			if(!given->is_string()){
				ret_string =  "Expects a string as input ";
				error = true;
			}
			break;
		case TYPE_TREESET: //Treeset
			if(!given->is_treeset()){
				ret_string =  "Expects a treeset as input ";
				error = true;
			}
			break;
		case TYPE_VECT: //Vect
			if(!given->is_vect()){
				ret_string =  "Expects a vector as input ";
				error = true;
			}
			break;
		case TYPE_INTVECT://Vector of integers
			if(!given->is_int() && given->is_vect()){
				ret_string =  "Expects an integer vector as input ";
				error = true;
			}
			break;
		case TYPE_STRINGVECT://Vector of strings
			if(!given->is_string() && given->is_vect()){
				ret_string =  "Expects a string vector as input ";
				error = true;
			}
			break;
		case TYPE_TREESETVECT://Vector of treesets
			if(!given->is_treeset() && !given->is_vect()){
				ret_string =  "Expects a vector of treesets as input ";
				error = true;
			}
			break;
		case TYPE_BOOL: //Bool
			if(!given->is_bool()){
				ret_string =  "Expects a boolean as input ";
				error = true;
			}
			break;
		case TYPE_CHAR: //Char
			if(!given->is_char()){
				ret_string =  "Expects a character as input ";
				error = true;
			}
			break;
		case TYPE_FLOAT: //Float
			if(!given->is_float()){
				ret_string =  "Expects a float as input ";
				error = true;
			}
			break;
		case TYPE_DOUBLE: //Double
			if(!given->is_double()){
				ret_string =  "Expects a double as input ";
				error = true;
			}
			break;
		case TYPE_SYMBOL: //Symbol
			if(!given->is_symbol()){
				ret_string =  "Expects a symbol as input ";
				error = true;
			}
			break;
		case TYPE_FUNCTION://Function
			if(!given->is_funct()){
				ret_string =  "Expects a function as input ";
				error = true;
			}
			break;
		default: //Not a known type given
			ret_string = "No known type given.";
			error = true;
			break;
	}
	return ret_string;
}

//Does typechecking for all functions input to add_function with a list of argument types
//simplifies type checking for all proper inputs (if you add a second add_function for
//the same afptr and different dataTypes type checking will be done with both in mind
//for the overloaded function. They will still have to be checked again in the user function
//for casting purposes as of this comment.)
pqlsymbol * u_template(vector <pqlsymbol *> arglist, string functName){

	pqlsymbol * result;

	//Checks to see if the function exists
	std::map<std::string, afptr>::iterator it = ptrMap.find(functName);
	if (it == ptrMap.end()){
		cout << "u_template error: Could not find function pointer\n";
		return new pqlsymbol(ERROR, "Invalid function");
	}
	//Checks that there are the proper number of inputs
	for(unsigned int i = 0; i < argMap[functName].size(); i ++){//for each potential set of arguments
		if (argMap[functName][i].size() != arglist.size()){
			cout << "Expects " << argMap[functName].size() << " inputs, Found: " << arglist.size() << endl;
			//			return new pqlsymbol(ERROR, "Input size Error");
		}
	}

	string errString = functName + " ";
	int inputNum = 0;

	//Checks for the proper types
	for(unsigned int j = 0; j < argMap[functName].size(); j++){//for each set of potential arguments
		bool error = false;

		for(unsigned int i = 0; i < argMap[functName][j].size(); i++){//For each arg type in argMap
			errString += type_check(argMap[functName][j][i], arglist[i], error);
			if(error){
				inputNum = i + 1;
				break;
			}
		}
		if(error){
			continue;
		}

		/* Leftover from attempts at working with overloaded functions
		   if(argMap[functName].size() > 1){
		   return result = new pqlsymbol( ( *it).second(arglist));
		   }	
		   */

		return ptrMap[functName](arglist);

	}
	//Outpus the error string accumulated
	cout << errString << inputNum <<  " Found : " << get_arg_types(arglist) << endl;
	return new pqlsymbol(ERROR, "Type Error");

}

//Overload of the add_function which accepts a variable number of inputs after the documentation string
//representing the types for each arguement accepted by the function (arguements should be of enumerated type
//dataType)
template <class... Args>
void add_function(string functname, afptr funct, string doc, Args... sepArgs){
	std::vector<dataType> args = {sepArgs...};
	::argFunctMap.insert(std::make_pair(functname, funct));
	//This is so tab-autocomplete works
	::functionKeys.push_back(functname);
	::helpRef.insert(std::make_pair(functname, doc));

	std::map<string, afptr>::iterator it;
	it = ::argFunctMap.find(functname);
	if(it == ::argFunctMap.end()){
		vector < vector < dataType > > vectargs;
		vectargs.push_back(args);
		::argMap.insert(std::make_pair(functname, vectargs));
	}
	else{
		::argMap[functname].push_back(args);
	}
	::ptrMap.insert(std::make_pair(functname, funct));
}

//adds the functions to the their maps. 
void init_the_functs()
{


	//Leftovers from an attempt to reduce the number of user functions which have
	//to be written
	/*
	   boost::function<void(void)> voidHelp;
	   voidHelp = boost::bind(::help, "void");
	   vector <boost::function<void(void)>> functVect;
	   functVect.push_back(voidHelp);
	   boost::bind(functVect[0], "rand_index");
	   voidHelp();
	   */

	/*
	   boost::variant<boost::function<void(std::string)> > varTest;
	   boost::variant< int, double > varInt;
	   varTest = boost::bind(help, _1);
	   std::string input = "show";
	   boost::get<boost::function<void(std::string)>>(varTest)(input);
	   */

	//::functionKeys.push_back("");

	//search

	//The New way
	add_function("search_ktet", &u_search_ktet, "Returns trees which contain the searched ktet");
	add_function("search_clade", &u_search_clade, "Returns trees which contain a clade containing the given taxa. This function has an optional strictness value. When turned on the search will only return trees with exact clades. ");




	add_function("average_ancestral_distance", &u_average_ancestral_distance, "Returns average distance to common ancestor given two taxa and an optional set of trees (no third input = all trees");
	add_function("get_trees_without_taxa", &u_get_trees_without_taxa, "Returns the trees that do not have the input taxa", TYPE_VECT);
	add_function("gtwot", &u_get_trees_without_taxa, "Returns the trees that do not have the input taxa", TYPE_VECT);

	add_function("clade_size_search", &u_clade_size_search, "Returns trees with clade of given size and taxa");
	add_function("clade_size_search", &u_clade_size_search, "Returns trees with clade of given size and taxa", TYPE_VECT, TYPE_INT);
	add_function("smallest_clade_search", &u_smallest_clade_search, "Returns trees with the smallest clade of given taxa", TYPE_VECT);
	
	add_function("get_trees_with_taxa", &u_get_trees_with_taxa, "Returns the trees that have the input taxa", TYPE_VECT);
	add_function("gtwt", &u_get_trees_with_taxa, "Returns the trees that have the input taxa", TYPE_VECT);
	add_function("similarity_search", &u_similarity_search, "Returns the trees that are the most similar to the given tree in that they have the most shared bipartitions.", TYPE_STRING);

	//add_function("get_trees_by_taxa", &u_get_trees_by_taxa, " ");
	add_function("get_trees_by_subtree", &u_get_trees_by_subtree, "Returns the trees that contain the input subtree", TYPE_STRING);

	add_function("subtree_search", &u_get_trees_by_subtree, "Returns the trees that contain the input subtree", TYPE_STRING);
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
	add_function("not", &u_not, "Returns the difference between all trees and the input treeset.", TYPE_TREESET);
	add_function("union", &u_union, "Returns the union of two treesets.", TYPE_TREESET, TYPE_TREESET);
	add_function("intersection", &u_intersection, "Returns the intersection of two treesets.", TYPE_TREESET, TYPE_TREESET );
	add_function("difference", &u_difference, "Returns the difference of two treesets.", TYPE_TREESET, TYPE_TREESET);

	//convert types
	add_function("set", &u_set, "Converts a list of ints into a treeset.");
	add_function("to_newick", &u_to_newick, "Returns the newick string for the tree matching the input index.", TYPE_INT);

	//analysis
	add_function("count", &u_count, "Returns the number of objects in the treeset or list", TYPE_TREESET);
	add_function("unique", &u_unique, "Returns the a subset of trees from a given treeset each with a unique topology.", TYPE_TREESET);
	add_function("unique_biparts", &u_unique_biparts, "Returns the number of all unique bipartitions given a treeset", TYPE_TREESET);
	add_function("duplicates", &u_duplicates, "Returns the set of trees with are topologically equal to the input tree.", TYPE_TREESET);

	add_function("sample_trees", &u_sample_trees, "Returns a random sampling of the input treeset of the requested size", TYPE_TREESET, TYPE_INT);

	//clustering
	add_function("silhouette", &u_silhouette, "Returns the silhouette distance between given clusters of trees", TYPE_TREESETVECT, TYPE_STRING);
	add_function("rand_index", &u_rand_index, "Returns the rand index of two clusterings for a set of trees", TYPE_TREESETVECT, TYPE_TREESETVECT);
	add_function("adjusted_rand_index", &u_adjusted_rand_index, "Returns the adjusted rand index of two clusterings for a set of trees", TYPE_TREESETVECT, TYPE_TREESETVECT);
	add_function("agglo_clust", &u_agglo_clust, "Returns the agglomerative clustering of the given input set of trees.", TYPE_TREESET, TYPE_INT, TYPE_STRING);
	add_function("kmeans_clust", &u_kmeans_clust, "returns the k means clustering of the given input set of trees.", TYPE_TREESET, TYPE_INT, TYPE_STRING);
	add_function("dbscan_clust", &u_dbscan_clust, "Return the dbscan clustering of the given input set of trees. Last cluster is always noise (trees that don't fit a cluster based on current criteria", TYPE_TREESET, TYPE_INT, TYPE_INT, TYPE_STRING);
	add_function("burnin_clust", &u_burnin_clust, "User Function for testing a clustering method - Currently unfinished", TYPE_TREESET, TYPE_STRING);

	//quartets
	add_function("shared_quartets_strict", &u_shared_quartets_strict, "Returns quartets present in every tree of treeset",TYPE_TREESET);
	add_function("shared_quartets_majority", &u_shared_quartets_majority, "Returns quartets present in a majority of trees", TYPE_TREESET);
	add_function("print_quartets_from_tree", &u_print_quartets_from_tree, "Prints all quartets of a tree given its index as an int", TYPE_INT);
	add_function("print_conflicting_quartets", &u_print_conflicting_quartets, "Given two bipartitions, prints all conflicting quartets");
	add_function("num_conflicting_quartets", &u_num_conflicting_quartets, "Returns number of conflicting quartets across two given bipartitions");
	//add_function("get_num_quartets", &u_get_num_quartets, "Returns the number of quartets implied by a given bipartition.", TYPE_INT);
	add_function("num_quartets", &u_get_num_quartets, "Returns the number of quartets implied by a given bipartition.", TYPE_INT);
	add_function("print_shared_quartets", &u_print_shared_quartets, "Prints the shared quartets between two given bipartitions");
	add_function("num_shared_quartets", &u_num_shared_quartets, "Returns number of shared quartets between two biparts");
	add_function("num_overlapping_quartets", &u_num_overlapping_quartets, "Returns number of overlapping quartets between two 	biparts");
	add_function("quartet_distance", &u_quartet_distance, "Returns quartet distance between two trees");
	add_function("print_quartet_distance", &u_print_quartet_distance, "Prints differing quartets between two trees.");
	add_function("conflicting_quartet_distance", &u_conflicting_quartet_distance, "Returns conflicting quartet distance between two trees");
	//distance
	add_function("distance_between_taxa", &u_distance_between_taxa, "Returns distance between two specified taxa in a specified tree");
	add_function("distance_to_common_ancestor", &u_distance_to_common_ancestor, "Returns distance between two specified taxa in a specified tree");
	add_function("distance_to_root", &u_distance_to_root, "Returns distance to root of specified taxa in specified tree");

	add_function("average_depth", &u_average_depth, "Returns the average depth of a specified taxa among all trees");
	add_function("expected_average_depth", &u_average_depth, "Returns the theoretical expected average depth given a number of taxa");
	add_function("average_distance_between_taxa", &u_average_distance_between_taxa, "Returns the average distance between two taxa among all trees");
	add_function("is_bifurcating", &u_is_bifurcating, "Returns boolean representing whether inputted tree is bifurcating");
	add_function("calculate_C", &u_calculate_C, "Returns the value of C (a symmetry measure) for a tree. Input is either an int or a string. Input tree MUST be bifurcating!");
	add_function("edit_distance_greedy", &u_edit_distance_greedy, "Given two trees, returns greedy edit distance");
	add_function("edit_distance_minimum", &u_edit_distance_minimum, "Given two trees, returns minimum edit distance");
	add_function("edit_distance_total", &u_edit_distance_total, "Given two trees, returns total edit distance");
	add_function("edit_distance_average", &u_edit_distance_average, "Given two trees, returns average edit distance");

	//add_function("average_depth", &u_average_depth, "Returns the average depth of a specified taxa among all trees",TYPE_TREESET);
	//add_function("average_distance_between_taxa", &u_average_distance_between_taxa, "Returns the average distance between two taxa among all trees",TYPE_INT, TYPE_INT, TYPE_TREESET);

	//consensus
	add_function("consensus", &u_consen, "Returns the newick string for the consensus tree for the input treeset.", TYPE_TREESET);
	add_function("consensus_strict", &u_strict_consen, "Returns the newick string for the strict consensus tree for the input treeset.", TYPE_TREESET);
	add_function("consensus_majority", &u_majority_consen, "Returns the newick string for the majority consensus tree for the input treeset.", TYPE_TREESET, TYPE_DOUBLE);
	add_function("consensus_greedy", &u_greedy_consen, " ");	
	add_function("consensus_least_conflict", &u_least_conflict, " ");
	//should move to a single function that calls consensus_reso_rate, or reso rate depending on input.
	add_function("consensus_reso_rate", &u_consensus_reso_rate, "Returns the consensus resolution rate for a set of trees and a given consensus strictness.");
	add_function("crr", &u_consensus_reso_rate, "Returns the consensus resolution rate for a  set of trees and a given consensus strictness.");
	add_function("reso_rate", &u_reso_rate, "Returns the resolution rate for a given tree.");



	//utilities

	add_function("help", &u_help, "Prints function description given a function name", TYPE_STRING);
	add_function("help", &u_help, "Prints function description given a function name");
	add_function("print_inverted_index", &u_print_inverted_index, "prints the inverted index, which is num trees long by num biparts wide.");


	//visualization
	add_function("display_clusters", &u_display_clusters, "Displays a gnuplot graph of the trees color coded by cluster", TYPE_STRING, TYPE_STRING, TYPE_TREESETVECT);
	add_function("display_heatmap", &u_display_heatmap, "Displays a gnuplot heatmap of the input trees using the given distance measure", TYPE_TREESET, TYPE_STRING, TYPE_STRING);
	add_function("show", &u_show, "Displays images of the specified tree or trees (Int, IntVect, or Treeset). Takes an optional String mode argument ('ortho' or 'radial') to display an SVG image. Default mode is text.");
	add_function("show_newick", &u_show_newick, "Displays images of the specified Newick strings (String or StringVect). Takes an optional String mode argument ('ortho' or 'radial') to display an SVG image. Default mode is text.");
	add_function("export", &u_export, "Exports images of the specified tree or trees (Int, IntVect, or Treeset) to the specified folder path (String). Takes an optional String mode argument ('ortho' or 'radial') to export an SVG image. Default mode is text.");
	add_function("export_newick", &u_export_newick, "Exports images of the specified Newick strings (String or StringVect) to the specified folder path (String). Takes an optional String mode argument ('ortho' or 'radial') to display an SVG image. Default mode is text.");

	//classification
	add_function("classification", &u_classification, "Prints the classification data for a given taxon or vector of taxa.", TYPE_STRINGVECT);
	add_function("write_classifications", &u_write_classifications, "Saves the classification data currently in memory to the specified filename.");
	add_function("edit_classification", &u_edit_classification, "Takes (taxa, rank, info). Edits the classification data for the specified taxon or vector of taxa by overwriting the taxon's existing info in the specified rank with the given info.");
	add_function("taxa_in_group", &u_taxa_in_group, "Returns the taxa belonging to the specified taxonomic group.", TYPE_STRING);
	add_function("show_group_in_tree", &u_show_group_in_tree, "Displays an SVG image of the given taxonomic group within the given tree(s). Highlights taxa in the clade that conflict with the expected grouping.");
	add_function("show_group", &u_show_group, "Displays an SVG image of the given taxonomic group isolated from its parent tree(s). Highlights taxa in the clade that conflict with the expected grouping.");
	add_function("show_level", &u_show_level, "If possible, displays an image of the given tree(s) abstracted to the level of the given taxonomic rank.");
	add_function("show_only", &u_show_only, "Computes a version of the given tree(s) including only taxa from the given group(s)");

	//traits
	add_function("load_trait_matrix", &u_load_trait_matrix, "Loads trait information from a NEXUS character matrix.", TYPE_STRING);
	add_function("test_trait_correlation", &u_test_trait_correlation, "Runs algorithm to test correlated evolution of specified traits", TYPE_INT, TYPE_INT, TYPE_INT,TYPE_INT,TYPE_INT, TYPE_INT,TYPE_STRING);

	//modify treeset
	//add_function("taxa_filter", &u_taxa_filter, "Creates a new tree or trees by removing all but the specified taxa from the specified tree(s). Original trees are kept intact.");
	//add_function("group_filter", &u_group_filter, "Creates a new tree or trees by removing all taxa except the ones in the specified taxonomic groups from the specified tree(s). Original trees are kept intact.");
	//add_function("delete_tree", &u_delete_tree, "Deletes the specified tree(s) from the working data set.");
	add_function("write_trz", &u_write_trz, "Writes specified Tree/TreeVect/Treeset to a .trz file with specified filename.");

	//Developer Functions
	add_function("debug", &u_debug, " ");
	add_function("print_taxa_in_trees", &u_print_taxa_in_trees, " ");
	add_function("print_biparttable", &u_print_biparttable, " ");
	add_function("set_hetero", &u_set_hetero, " ", TYPE_BOOL);
	add_function("rsearch", &u_random_search, " ");
	add_function("rsearch2", &u_random_search2, " ");
	
	
	
	add_function("tier", &u_tier, " ", TYPE_INT);
	add_function("tttier", &u_tttier, " ", TYPE_INT);
	add_function("rate_tiers", &u_rate_tiers, " ", TYPE_INT);
	add_function("print_label_map", &u_print_label_map, " ");
	add_function("get_taxa_in_tree", &u_get_taxa_in_tree, "Returns the taxa in a given tree. Takes a tree index", TYPE_INT);
	add_function("distinguishing_taxa", &u_distinguishing_taxa, "Takes two Treesets. Returns the taxa that are in of one treeset and none of the other.", TYPE_TREESET, TYPE_TREESET);
	add_function("proto", &u_prototype, "This is a protype function only use for the function in current developemnt.");
	add_function("homog", &u_homogenize, "This is a protype function which homogenizes a tree set removing all taxa which do not appear in all of the trees. Extra branches are collasped. This is currently a destructive function and non-reversable.");
	add_function("print_cladetable", &u_print_cladetable, " ");
	add_function("print_bipartitions", &u_print_bipartitions, " ");
	add_function("print_clades", &u_print_clades, " ");

	
	
	//GRB NEW Distances
	add_function("rf_dist_clade", &u_rf_dist_clade, " ");
	add_function("rf_dist_bipart", &u_rf_dist_bipart, " ");

	//GRB NEW TESTING FUNCTIONS!
	add_function("rsk", &u_rsk, " ");
	add_function("ssk", &u_ssk, " ");
	add_function("fsk", &u_fsk, " ");
	
	add_function("rssc", &u_rssc, " ");
	add_function("sssc", &u_sssc, " ");
	add_function("fssc", &u_fssc, " ");
	
	add_function("ssc", &u_ssc, " ");
	add_function("fsc", &u_fsc, " ");	

	cout << "Functions loaded into map" << endl;




}
