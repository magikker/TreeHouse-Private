#ifndef _THGLOBALS_H_
#define _THGLOBALS_H_

#include <map>
#include <ctime>
#include <limits>

#include "label-map.hh"
#include "pqlsymbol.h"
#include "buildtree.h"
#include "BipartitionTable.h"
#include "Taxon.h"
#include "TreeTable.h"
#include "global.h"
#include "UtilityFunctions.h"

using namespace std;

//class pqlsymbol;

typedef pqlsymbol * (*vfptr)();
typedef pqlsymbol * (*afptr)(vector < pqlsymbol * >);
typedef std::map<std::string, vfptr> voidFuncts;
typedef std::map<std::string, afptr> argFuncts;

//~ typedef struct {
  //~ std::string name;			/* User printable name of the function. */
  //~ afptr *func;				/* Function to call to do the job. */
  //~ std::string doc;			/* Documentation for this function.  */
//~ } command;
//~ 
//~ typedef std::map<std::string, command> Functs;



extern bool DEBUGMODE;

//extern vector<int> shuffledTaxa;

//extern bool hetero; // = 0; // If the trees are taxa heterogenious
//extern unsigned num_trees; // = 0; // number of trees
//extern unsigned num_taxa; //  = 0; // number of taxa

extern BipartitionTable biparttable;
//extern TreeTable treetable;

extern LabelMap lm;										// A mapping between taxa name and placement in bitstring

extern vector<Taxon *> taxa_info;
//These objects are number of total trees long
extern vector< bool* > taxa_in_trees;					//Which taxa are in which trees // bitstrings 1's are in the tree. place in the bitsting maps to the Lable map.

extern vector< vector< unsigned int> > inverted_index;	//???

extern set< unsigned int > all_trees;							//stores the number for each tree. This is the world as to be used to determine the result of not operations. // This should be moved to being returned from a function call.
extern set< unsigned int > original_trees;
extern int NUM_TREES_INIT;  

extern std::map<int, vector<int> > tree_dups;

extern vector< pqlsymbol* > query_results;

//Function maps. Use void if the function takes no args and arg if it does. 
extern voidFuncts voidFunctMap;
extern argFuncts argFunctMap;
extern std::vector<std::string> functionKeys;

//map of the symbol in the language. These are were declared varibles are stored. If a user enters r = 7; r is set to 7 in the map. 
extern std::map<string, pqlsymbol * > symbol_table;
extern std::map<string, bool > constant_table;

extern ofstream output_file;


extern vector<std::clock_t> clocks;


extern unsigned int SetOps;
extern unsigned int SetInsertions;
extern double SetTime;
extern double HetTime;
extern double AHetTime;
extern double SearchTime;


bool write_to_output(string input);
bool change_output_pointer(string input);
bool init_output();

bool are_taxa_in_tree(int treeindex, vector<int> setoftaxa);

int num_taxa_in_tree(int treeindex);
vector<string> get_taxa_in_tree(unsigned int treeindex);
vector<string> get_all_taxa_vect();
set<string> get_all_taxa_set();
vector<int> get_taxa_with_trait(unsigned int trait_id, int trait_value);
vector<int> get_taxa_without_trait(unsigned int trait_id, int trait_value);
vector<int> get_taxa_in_clade(vector<int> taxa, unsigned int tree);
int distance_between_taxa(unsigned int taxon1, unsigned int taxon2, unsigned int tree);
int distance_to_common_ancestor(unsigned int taxon1, unsigned int taxon2, unsigned int tree);

string to_lower(string str);
int index_in_labelmap(string label);
vector< bool *> get_tree_bipartitions(unsigned int id);
vector<unsigned int> get_tree_bs_sizes(unsigned int id);
vector<float> get_tree_branches(unsigned int id);
vector<unsigned int> get_tree_data(unsigned int id, vector < bool *>& tree_bipartitions, vector <unsigned int>& tree_bs_sizes, vector<float>& tree_branches);
string th_compute_tree(BipartitionTable& bpt, unsigned id, bool branch);
void recompute_tree_dups();

vector<unsigned int> which_trees_single_strict(int bitstringindex, vector<int> positions);
vector<unsigned int> which_trees_double_strict(int bitstringindex, int numTaxaSearched);
vector<unsigned int> which_strict_hetero(int bitstringindex, vector<int> positions1, vector<int> positions2);
vector<unsigned int> which_strict_hetero(int bitstringindex, vector<int> positions1, vector<int> positions2, int &bitcomps);
vector<unsigned int> which_strict_hetero(int bitstringindex, vector<int> positions);
vector<unsigned int> which_strict_hetero(int bitstringindex, vector<int> positions, int &bitcomps);

vector<unsigned int> which_hetero(int bitstringindex, vector<int> positions);
vector<unsigned int> which_hetero(int bitstringindex, vector<int> positions, int &bitcomps);

bool is_strict_homog(int bitstringindex, vector<int> positions);
bool is_strict_homog(int bitstringindex, vector<int> positions, int &bitcomps);

bool is_taxa_homogenious(set<unsigned int> treeset);
vector<string> get_common_taxa();
vector<string> get_mask_for_homog();

bool is_compat(bool *bitstring1, int length1, bool *bitstring2, int length2);
bool is_compat(int bitstringindex1, int bitstringindex2);
bool is_compat(bool *bitstring1, int length1, int bitstringindex2);
vector<unsigned int> find_incompat(bool *bitstring, int length, set<unsigned int> inputtrees);
vector<unsigned int> find_incompat_old(bool *bitstring, int length, set<unsigned int> inputtrees);
BipartitionTable get_consen_bt(set<unsigned int> inputtrees, float percent);


//void set_hetero(bool input);
//bool get_hetero();
//unsigned int get_num_trees();
//void set_num_trees(unsigned int input);
//unsigned int get_num_taxa();
//void set_num_taxa(unsigned int input);

void print_taxa_in_trees();

void debugstatement(string input);

string get_time_stamp();
void start_clock();
double stop_clockbp();
double stop_clockb();
double stop_clockf();

#endif
