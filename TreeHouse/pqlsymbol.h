#ifndef _PQLSYMBOL_H
#define _PQLSYMBOL_H

// Symbols Symboltables and Bears oh my. 
#include <iostream>
#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <sstream>
#include <string.h>

using namespace std;

typedef enum {ATOM, LIST, FUNCT, ERROR, SYMBOL} object_type;
typedef enum {THE_EMPTY_LIST, BOOLEAN, CHAR, INT, FLOAT, DOUBLE, STRING, VECTOR, TREESET} data_type;

class pqlsymbol
{
	private:
	
	object_type otype;
	data_type dtype;
		
	vector<pqlsymbol * > list;
	union {
		bool boolvalue;	
		char charvalue;
		int intvalue;
		float floatvalue;
		double doublevalue;
		std::string stringvalue;
		std::string symbolname; 
	};
		unsigned int ntrees;
		std::set<unsigned int> treeset;		


	public:
	pqlsymbol(){
		dtype = THE_EMPTY_LIST;
		otype = LIST;
	}

	pqlsymbol(std::set<unsigned int> inval, unsigned int num_trees){
		dtype = TREESET;
		otype = ATOM;
		ntrees = num_trees;
		
		for(set<unsigned int>::const_iterator pos = inval.begin(); pos != inval.end(); ++pos){
			if(*pos < ntrees && *pos >= 0){
				treeset.insert(*pos);
			}
		}
	}

	pqlsymbol(bool inval){
		dtype = BOOLEAN;
		otype = ATOM;
		boolvalue = inval;
	}

	pqlsymbol(vector<bool> inval){
		dtype = BOOLEAN;
		otype = LIST;
		for(std::vector<bool>::size_type i = 0; i != inval.size(); i++){
			list.push_back(new pqlsymbol(inval [i]));
		} 
	}

	pqlsymbol(char inval){
		dtype = CHAR;
		otype = ATOM;
		charvalue = inval;
	}

	pqlsymbol(vector<char> inval){
		dtype = CHAR;
		otype = LIST;
		for(std::vector<char>::size_type i = 0; i != inval.size(); i++){
			list.push_back(new pqlsymbol(inval [i]));
		} 
	}

	pqlsymbol(int inval){
		dtype = INT;
		otype = ATOM;
		intvalue = inval;
	}

	pqlsymbol(vector<int> inval){
		dtype = INT;
		otype = LIST;
		for(std::vector<int>::size_type i = 0; i != inval.size(); i++){
			list.push_back(new pqlsymbol(inval [i]));
		} 
	}

	pqlsymbol(float inval){
		dtype = FLOAT;
		otype = ATOM;
		floatvalue = inval;
	}

	pqlsymbol(vector<float> inval){
		dtype = FLOAT;
		otype = LIST;
		for(std::vector<float>::size_type i = 0; i != inval.size(); i++){
			list.push_back(new pqlsymbol(inval [i]));
		} 
	}

	pqlsymbol(double inval){
		dtype = DOUBLE;
		otype = ATOM;
		doublevalue = inval;
	}

	pqlsymbol(vector<double> inval){
		dtype = DOUBLE;
		otype = LIST;
		for(std::vector<double>::size_type i = 0; i != inval.size(); i++){
			list.push_back(new pqlsymbol(inval [i]));
		} 
	}

	pqlsymbol(vector<string> inval){
		dtype = STRING;
		otype = LIST;
		for(std::vector<int>::size_type i = 0; i != inval.size(); i++){
			list.push_back(new pqlsymbol(inval [i]));
		} 
	}

	pqlsymbol(vector< set <unsigned int>> inval){
		dtype = TREESET;
		otype = LIST;
			for(std::vector<int>::size_type i = 0; i != inval.size(); i++){
				list.push_back(new pqlsymbol(inval[i], inval[i].size()));
			}
	}

	pqlsymbol(const std::string& instring){
		dtype = STRING;
		otype = ATOM;
		new (&stringvalue) string(instring);
	}

	pqlsymbol(const std::string&& instring){
		dtype = STRING;
		otype = ATOM;
		new (&stringvalue) string(std::move(instring));
	}

	pqlsymbol(object_type intype, const string& instring){
		dtype = STRING;
		otype = intype;
		new (&stringvalue) string(instring);
	}

	pqlsymbol(object_type intype, const string&& instring){
		dtype = STRING;
		otype = intype;
		new (&stringvalue) string(std::move(instring));
	}
	
	pqlsymbol(vector< pqlsymbol * > inval){	
		if (inval.size() > 0){
			dtype = inval[0]->dtype;
			otype = LIST;
		}
		else{	
			dtype = THE_EMPTY_LIST;
			otype = LIST;
		}
		list = std::move(inval);
	}

	~pqlsymbol(){
        // call destructor of string if needed
		
		if (dtype == STRING && otype == ATOM){
			 stringvalue.~string();
		}
		 
		if (otype == LIST){
			while(!list.empty()){
				delete list.back();
				list.pop_back();
			}
		}
    }

//Gets

	pqlsymbol * get_value_as_new_pqlsymbol(){
		
		if (otype == ERROR){
			return new pqlsymbol(ERROR, get_string() );
		}
		
		//if (otype == TREESET){
		//	return new pqlsymbol(get_treeset(), ntrees);
		//}
		
		if (otype == ATOM){
			switch (dtype){
				case BOOLEAN:
	    			return new pqlsymbol(get_bool() );

				case CHAR:
					return new pqlsymbol(get_char() );
	
				case INT:
	 		   		return new pqlsymbol(get_int() );
	
				case FLOAT:
	 		   		return new pqlsymbol(get_float() );

				case DOUBLE:
		    		return new pqlsymbol(get_double() );

		  		case STRING:
		    		return new pqlsymbol(get_string() );

				case TREESET:
				return new pqlsymbol(get_treeset(), ntrees );

		    	case THE_EMPTY_LIST:
		    		return new pqlsymbol();

				default:
	    			return new pqlsymbol(ERROR, "Shouldn't hit this in [] return") ;	
			}
			
		}
		
		if (otype == LIST){
			switch(dtype){
				
				case BOOLEAN:
	    			return new pqlsymbol(get_bool_vect() );

				case CHAR:
					return new pqlsymbol(get_char_vect() );
	
				case INT:
	 		   		return new pqlsymbol(get_int_vect() );
	
				case FLOAT:
	 		   		return new pqlsymbol(get_float_vect() );

				case DOUBLE:
		    		return new pqlsymbol(get_double_vect() );

		  		case STRING:
		    		return new pqlsymbol(get_string_vect() );

				case TREESET:
				return new pqlsymbol(get_treeset_vect() );
		    		
		    	case THE_EMPTY_LIST:
		    		return new pqlsymbol();
	//			I need to implement a deep copy of the vector. 			
	//			case VECTOR:
	//				return shared_ptr<pqlsymbol>( new pqlsymbol() );

				default:
	    			return new pqlsymbol(ERROR, "Shouldn't hit this in [] return") ;	
				
			}
		}
	}
	
	int get_size(){
		if (dtype == TREESET){
			return (int)treeset.size();
		}
		
		if (otype == ATOM){
			return 1;
		}
		
		if (otype == LIST){
			return list.size();
			}
		
		return -1;
	}

	std::set<unsigned int> get_treeset(){
		return treeset;
	}
	
	bool get_bool(){
		return boolvalue;
	}

	bool get_bool_from_vect(int i){
		return list[i]->boolvalue;
	}

	vector<bool> get_bool_vect(){
		vector<bool> rval;
		for(std::vector<pqlsymbol>::size_type i = 0; i != list.size(); i++){
			rval.push_back(list[i]->boolvalue);
		} 
		return rval;
	}
		
	char get_char(){
		return charvalue;
	}	

	char get_char_from_vect(int i){
		return list[i]->charvalue;
	}	

	vector<char> get_char_vect(){
		vector<char> rval;
		for(std::vector<pqlsymbol>::size_type i = 0; i != list.size(); i++){
			rval.push_back(list[i]->charvalue);
		} 
		return rval;
	}

	int get_int(){
		return intvalue;
	}

	int get_int_from_vect(int i){
		return list[i]->intvalue;
	}	

	vector<int> get_int_vect(){
		vector<int> rval;
		for(std::vector<pqlsymbol>::size_type i = 0; i != list.size(); i++){
			rval.push_back(list[i]->intvalue);
		} 
		return rval;
	}

	vector <set <unsigned int > > get_treeset_vect(){
		vector<set < unsigned int > > rval;
		for(std::vector<pqlsymbol>::size_type i = 0; i != list.size(); i++){
			rval.push_back(list[i]->get_treeset());
		}
		return rval;
	}


	float get_float(){
		return floatvalue;
	}
	
	float get_float_from_vect(int i){
		return list[i]->floatvalue;
	}	

	vector<float> get_float_vect(){
		vector<float> rval;
		for(std::vector<pqlsymbol>::size_type i = 0; i != list.size(); i++){
			rval.push_back(list[i]->floatvalue);
		} 
		return rval;
	}

	double get_double(){
		return doublevalue;
	}

	double get_double_from_vect(int i){
		return list[i]->doublevalue;
	}	

	vector<double> get_double_vect(){
		vector<double> rval;
		for(std::vector<pqlsymbol>::size_type i = 0; i != list.size(); i++){
			rval.push_back(list[i]->doublevalue);
		} 
		return rval;
	}

	std::string get_string(){
		return stringvalue;
	}
	
	string get_string_from_vect(int i){
		return list[i]->stringvalue;
	}	
	
	vector<string> get_string_vect(){
		vector<string> rval;
		for(std::vector<pqlsymbol>::size_type i = 0; i != list.size(); i++){
			rval.push_back(list[i]->stringvalue);
		}
		return rval;
	}
	
	std::string get_symbolname(){
		return symbolname;
	}

	int get_vect_size(){
		return list.size();
	}

	int get_treeset_size(){
		return treeset.size();
	}
//~ 
	//~ set<int> * get_int_vect_as_set_pointer(){
		//~ vector<int> tempvect = get_int_vect();
		//~ 
		//~ 
		//~ return new set<int>(tempvect.begin(), tempvect.end());
	//~ }

	pqlsymbol * get_list_symbol(int i){
		
		if (otype ==  LIST){
			switch (dtype){
				case BOOLEAN:
	    			return new pqlsymbol(list[i]->get_bool() );

				case CHAR:
					return new pqlsymbol(list[i]->get_char() );
	
				case INT:
	 		   		return new pqlsymbol(list[i]->get_int() );
	
				case FLOAT:
	 		   		return new pqlsymbol(list[i]->get_float() );

				case DOUBLE:
		    		return new pqlsymbol(list[i]->get_double() );

		  		case STRING:
		    		return new pqlsymbol(list[i]->get_string() );

				case TREESET:
				return new pqlsymbol(list[i]->get_treeset(), ntrees );
	
	//			I need to implement a deep copy of the vector. 			
	//			case VECTOR:
	//				return shared_ptr<pqlsymbol>( new pqlsymbol() );

				default:
	    			return new pqlsymbol(ERROR, "Shouldn't hit this in [] return") ;	
			}
		}	
	}

//	vector< shared_ptr< pqlsymbol > > get_list_symbol(){
//		return list;
//	}


	bool is_treeset(){
		return dtype == TREESET;
	}

	bool is_bool(){
		return dtype == BOOLEAN;
	}

	bool is_char(){
		return dtype == CHAR;
	}
	
	bool is_int(){
		return dtype == INT;
	}

	bool is_float(){
		return dtype == FLOAT;
	}

	bool is_double(){
		return dtype == DOUBLE;
	}	

	bool is_string(){
		return dtype == STRING;
	}

	bool is_symbol(){
		return otype == SYMBOL;
	}

	bool is_emptylist(){
		return dtype == THE_EMPTY_LIST;
	}

	bool is_vect(){
		return otype == LIST;
	}

	bool is_atom(){
		return otype == ATOM;
	}

	bool is_funct(){
		return otype == FUNCT;
	}

	bool is_error(){
		return otype == ERROR;
	}

	object_type get_object_type(){
		return otype;
	}
	
	data_type get_data_type(){
		return dtype;
	}

	string value_to_string(){
		std::stringstream out;
		
		if (otype ==  ERROR){
			out << get_string();
	    	return out.str();
		}
				
		if (otype ==  ATOM){

			switch (dtype){
				case BOOLEAN:
					out << get_bool();
	    			return out.str();

				case CHAR:
				{
					out << "\'"	;
					out << get_char();
					out << "\'"	;
	 		   		return out.str();
				}
				case INT:
				{
					out << get_int();
	 		   		return out.str();
				}
				case FLOAT:
				{	
					out << get_float();
	 		   		return out.str();
				}
				case DOUBLE:
				{
					out << get_double();
		    		return out.str();
				}
		  		case STRING:
				{
					out << "\"";
					out << get_string();
					out << "\"";
					return out.str();
				}
				case TREESET:
				{
					//~ cout << "what's here? = ";
					//~ std::copy(treeset.begin(), treeset.end(), std::ostream_iterator<int>(std::cout, ", ") );
					//~ cout << endl;
			
					out << "{";
					std::copy(treeset.begin(), treeset.end(), std::ostream_iterator<int>(out, ", "));
				
					//~ // pre-increment and pre-decrement are faster than post-increment and post-decrement...
					//~ for(set<int>::const_iterator pos = treeset.begin(); pos != treeset.end(); ++pos){
					//~ int temp = *pos;
					//~ out << (int)*pos << ', ';
					//~ cout << *pos << endl;
					//~ }
					string temp = out.str();
					if (temp.size () > 1)  temp.resize (temp.size () - 2);
			
						temp.append("}");
			
						//~ out << "}";		
						//~ string temp = out.str();
						//~ if (temp.size () > 0)  temp.resize (temp.size () - 2);
						return temp;
					//	return out.str();

				}
				case THE_EMPTY_LIST:
					return "";

				default:
	    			return "****shouldn't hit this thing ever!!! ****** ATOM";	
					
			}
		}
		if (otype == LIST){
			out << "[";
			for(std::vector<pqlsymbol * >::size_type i = 0; i != list.size(); i++){								
					out << list[i]->value_to_string();
					if (i+1 < list.size()){
							out << ", ";
					}
			}

			out << "]";		
			return out.str();
		}
	print_dtype();
	print_otype();
	return "*************NO REALLY THIS IS BAD";
	}


	void print_dtype(){
		cout << "dtype = ";
		switch (dtype){
			
				case THE_EMPTY_LIST: 
					cout << "THE_EMPTY_LIST" << endl;
					break;
					
				case BOOLEAN:
					cout << "BOOLEAN" << endl;
					break;
					
				case CHAR:
					cout << "CHAR" << endl;
					break;
						
				case INT:
					cout << "INT" << endl;
					break;
					
				case FLOAT:
					cout << "FLOAT" << endl;
					break;
						
				case DOUBLE:
					cout << "DOUBLE" << endl;
					break;
					
		  		case STRING:
					cout << "STRING" << endl;
					break;
					
				case VECTOR:
					cout << "VECTOR" << endl;
					break;

				case TREESET:
					cout << "TREESET" << endl;
					break;
		}
	}

	void print_otype(){
		cout << "otype = ";
		switch (otype){
			
				case ATOM: 
					cout << "ATOM" << endl;
					break;
									
				case LIST:
					cout << "LIST" << endl;
					break;
					
				case FUNCT:
					cout << "FUNCT" << endl;
					break;
						
				case ERROR:
					cout << "ERROR" << endl;
					break;
					
				case SYMBOL:
					cout << "SYMBOL" << endl;
					break;
		}
	}


};


#endif

