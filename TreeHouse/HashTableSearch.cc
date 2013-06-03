#include <sstream>
#include <iterator>
#include <iostream>
#include <string>
#include <set>
#include <ctime>
#include <algorithm>

#include <stdlib.h>
#include <string.h>

#include "compressfunc.h"
#include "global.h"
#include "THGlobals.h"
#include "label-map.hh"
#include "pql.h"

#include "newick.h"

#include "HashTableSearch.h"


using namespace std;

//I need to call a function on a newick string
	// This function will be called in the parser
	// Returns a set of bipartitions 

	//I need to turn the newick string into a a NEWICK tree to be processed. 
	// Write a fuction for newick.c that takes a string. 

	//newick tree needs to go into a function that returns the bipartitons
	//or the subtrees... Or the Ktets
	//or something usefull.... 

	// I need a function up front that will check if all the labels are in the tree set... If not I can quit right there. 

void handle_newick_error(int err){
  switch (err) {
  case -1:
    printf("Out of memory\n");
    break;
  case -2:
    printf("parse error\n");
    break;
  case -3:
    printf("Can't load file\n");
    break;
  default:
    printf("Error %d\n", err);
  }
  exit(0);
}

bool * dfs_compute_bitstrings(NEWICKNODE* startNode, NEWICKNODE* parent, vector< vector < unsigned > > &solution )
{
  if (HETERO && !HCHECK)
	{
      cout << "if (HETERO && !HCHECK)" << endl;
      return NULL;
	}

  vector <unsigned int> subtree;

  // If the node is leaf node, just set the place of the taxon name in the bit string to '1' 
  //bitstrings here are not used to compute hash values -- however, we do collect them.
  startNode->myparent = parent;
  ///fprintf(stderr, "AT NODE: ");
  if (startNode->Nchildren == 0) 
  { // leaf node
    string temp(startNode->label);
    unsigned int idx = ::lm[temp];
    bool * bs = new bool[::NUM_TAXA];
	
    for (unsigned int i = 0; i < ::NUM_TAXA; i++)
	{
      bs[i] = 0;
	}
    bs[idx] = true;

	

	//string attempt(bs);	
	
	//Error Checking
    int numOnes = 0;
    for (unsigned int i = 0; i < ::NUM_TAXA; i++)
	{
      numOnes+=bs[i];
	  if (bs[i] ==  true)
	  {
        subtree.push_back(i);
      }
	}
    if (numOnes > 1 || numOnes == 0) 
	{
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring creation error\n";
      exit(0);
    }
	//End Error Checking

	//print_bitstring(bs, ::NUM_TAXA);
	solution.push_back(subtree);
	//cout << "Here we go" << endl;
	//cout << solution.size() << endl;
	//print_vector_of_bs(solution, ::NUM_TAXA);
    return bs;
  }

  else 
  {
    bool * ebs[startNode->Nchildren];
    //fprintf(stderr, "I am an internal node!\n");
    for (int i=0; i<startNode->Nchildren; ++i) 
    {
      ebs[i] = dfs_compute_bitstrings(startNode->child[i], startNode, solution);
    }
    
    bool * bs = new bool[::NUM_TAXA];
    for (unsigned int i  = 0; i < ::NUM_TAXA; i++)
      bs[i] = 0;

    ///fprintf(stderr, "the bitstrings I'm going to be comparing are:\n");
    /*for (int i=0; i<startNode->Nchildren; ++i) {
      fprintf(stderr, "%s\n", (startNode->child[i])->label);
      }
      fprintf(stderr, "\n");*/

    for (int i=0; i<startNode->Nchildren; ++i) 
	{
      //cerr << "bitstring: " << *(ebs[i]) << endl;
	  if (ebs[i]) 
	  {
	    for (unsigned int j = 0; j < ::NUM_TAXA; j++)
	      bs[j] |= ebs[i][j];
	    delete [] ebs[i];
	    ebs[i] = NULL;
	  }
	  else 
	  {
	    if (!HETERO)
		{
	      cout << "ERROR: null bitstring\n";
	      exit(0);
        }
	  }
    }

   int numOnes = 0;
   for (unsigned int i= 0; i < ::NUM_TAXA; i++)
   {
     numOnes+=bs[i];
	 if (bs[i] ==  true)
	 {
       subtree.push_back(i);
     }
   }
    if (numOnes < 1) {
      cout << "numOnes = " << numOnes << endl;
      cerr << "ERROR: bitstring OR error\n";
      exit(0);
    }

	//print_bitstring(bs, ::NUM_TAXA);
	solution.push_back(subtree);
	//cout << "Here we go" << endl;
	//cout << solution.size() << endl;
	//print_vector_of_bs(solution, ::NUM_TAXA);
    return bs;
  } //end else
}


set<unsigned int> search_hashtable_ktets(vector< bool * > list_bs, vector < vector < unsigned int > > subtrees)
{
  //cout << "We show the bipartitions found by the search and corresponding trees" << endl;
  //keep in mind that hashtable and hash_lengths are global variables

  vector<unsigned int> all_trees = returnAllTrees(NUM_TREES);
  set<unsigned int> total_trees;
  
  copy(all_trees.begin(), all_trees.end(), inserter(total_trees, total_trees.end()));  


  for(unsigned int i = 0; i < subtrees.size(); i++) // for each subtree
  {
	 set<unsigned int> trees;
     
     for (unsigned int j = 0; j < list_bs.size(); j++) // for each bipartition in the hashtable
	 {
	   bool foundFlag = true;
       for(unsigned int k = 0; k < subtrees[i].size(); k++) // for each piece of the subtree
	   {
         if(subtrees[i][k] >= ::length_of_bitstrings[j] || list_bs[j][subtrees[i][k]] == 0 )
		 {
		   cout << "Hit it" << endl;
           foundFlag= false;
		   //break;
		 }

       }

       if (foundFlag == true)
	   {
         for (unsigned int l = 0; l < hash_lengths[j]; l++)
	     {	
	       trees.insert(hashtable[j][l]);
	     }
       }

     }

    //Printing
    cout << "the search is ";
	for(unsigned int k = 0; k < subtrees[i].size(); k++) // for each piece of the subtree
		cout <<  subtrees[i][k] << " "; 
	cout << endl;
	cout << "the trees I found are : ";
	//print_set(trees);
	//end printing

    std::set<unsigned int>::iterator it1 = total_trees.begin();
    std::set<unsigned int>::iterator it2 = trees.begin();
    
    while ( (it1 != total_trees.end()) && (it2 != trees.end()) ) 
      {
        if (*it1 < *it2) 
          {
            total_trees.erase(it1++);
          } 
          else if (*it2 < *it1) 
          {
            ++it2;
    	  } 
          else 
          { 
		  // *it1 == *it2
            ++it1;
            ++it2;
    	  }
	  }
	  // Anything left in set_1 from here on did not appear in set_2,
	  // so we remove it.
	  total_trees.erase(it1, total_trees.end());
	  //$trees=$a.trees;
	  cout << "len of total_trees " << total_trees.size() << endl;
	  cout << endl;
  }

  return total_trees;




/*
  for (unsigned int i = 0; i < list_bs.size(); i++)
  { //Each bipartition in the treeset

    foundFlag = true; //assume true and check each case for a contridiction to set it as false. If it passes the tests then it is true. 
    
    for(unsigned int j = 0; j < querys.size(); j++)
	{
	  for(unsigned int k = 0; k < querys[j].size(); k++)
		{
		  if( querys[j][k] >= ::length_of_bitstrings[i] || list_bs[querys[j][k]] == 0 )
		  {
            foundFlag= false;
			break;
		  }
		}
	}

	//if we got here and it's true it passed all the checks. 
	if (foundFlag == true)
	{
      for (unsigned int l = 0; l < hash_lengths[i]; l++)
	  {
		
		trees.insert(hashtable[i][l]);
	  }
      //cout << "]" << endl;
    }
  }
  return trees;

*/

}



set<unsigned int> subtree_to_ktets(string subtree, vector< bool * > list_bs)
{
	set<unsigned int> trees;
	NEWICKTREE *newickTree;
	int err;
	char * cs = strdup(subtree.c_str());
	newickTree = stringloadnewicktree(cs, &err);

	vector< vector < unsigned int > > subtrees; 

	if (!newickTree) 
	{
		handle_newick_error(err);
    }

    else 
	{
		//unsigned int numBitstr=0;
		dfs_compute_bitstrings(newickTree->root, NULL, subtrees);
		
  		cout << "We show bitstrings:" << endl;
  		//keep in mind that hashtable and hash_lengths are global variables
  		for (unsigned int i = 0; i < subtrees.size(); i++)
		{
    		for (unsigned int j = 0; j < subtrees[i].size(); j++)
			{
      			cout << subtrees[i][j] << " " ;
    		}
    		cout << endl;
  		}
		trees = search_hashtable_ktets(list_bs, subtrees);
		//print_vector_of_bs(bitstrings, ::NUM_TAXA);
		//if (MAXVAL == 0)
		//{
		//	MAXVAL = numBitstr -1;
		//}		
		killnewicktree(newickTree);
    }
	return trees;
}


/*
void search_hashtable(vector< bool * > list_bs, vector< bool * > &search_bs){
  cout << "We show the hash table below, with the corresponding bitstring reps of the bipartitions:" << endl << endl;
  //keep in mind that hashtable and hash_lengths are global variables
  bool foundFlag = true;
	
  for (unsigned int i = 0; i < search_bs.size(); i++){ //each search string
    for (unsigned int j = 0; j < NUMBIPART; j++){ // Each bipartition
      foundFlag = true;
      for (unsigned int k = 0; k < NUM_TAXA; k++){ // each bit
        //cout << search_bs[i][k] << " " << list_bs[j][k] << endl;
	if (search_bs[i][k] != list_bs[j][k]){
	  
	  foundFlag = false;
	  break;
	}
	
      }
      if (foundFlag == true){
        cout << "FOUND ONE" << endl;
        unsigned int mycurrsize = 0;
	for (unsigned int l = 0; l < NUM_TAXA; l++){
	  cout << list_bs[j][l];
	}
	cout << " --> [ ";
	mycurrsize = hash_lengths[j];
	for (unsigned int l = 0; l < mycurrsize; l++){
	  cout << hashtable[j][l] << " ";
	}
	cout << "]" << endl;
      }
    }
  }
}
*/
/*
set<int> search_hashtable2(vector< bool * > list_bs, vector<int> leftside, vector<int> rightside){
  //cout << "We show the bipartitions found by the search and corresponding trees" << endl;
  //keep in mind that hashtable and hash_lengths are global variables
  bool foundFlag = true;
  
  set<int> trees;
	
  for (unsigned int i = 0; i < list_bs.size(); i++){ //Each bipartition in the treeset
   
	foundFlag = true;
	
	// if the first taxa on the left and the first taxa on the right side of the search bipart are on the same side of the bipart
	// we know we can bail out before checking anything else. 
	if (list_bs[i][leftside[0]] == list_bs[i][rightside[0]]){ 
      foundFlag = false;
      continue;
    }
	
	//  
	//    
    for (unsigned int j = 0; j < leftside.size(); j++){ // Each taxa on the left side
		if (list_bs[i][leftside[j]] == list_bs[i][leftside[0]])
		  continue;
		else{
	  		foundFlag = false;
	  		break;
		}
    }
    
	for (unsigned int j = 0; j < rightside.size(); j++){ // Each taxa on the right side 
		if (list_bs[i][rightside[j]] == list_bs[i][rightside[0]])
	  		continue;
		else{
	  		foundFlag = false;
	  		break;
		}
    }
    
	if (foundFlag == true){
      //cout << "FOUND ONE" << endl;
      unsigned int mycurrsize = 0;
      //for (unsigned int l = 0; l < NUM_TAXA; l++){
	  //cout << list_bs[i][l];
      //}
      //cout << " --> [ ";
      mycurrsize = hash_lengths[i];
      for (unsigned int l = 0; l < mycurrsize; l++){
		trees.insert(hashtable[i][l]);
      //  cout << hashtable[i][l] << " ";
      }
      //cout << "]" << endl;
    }
  }
  return trees;
}
*/

set<unsigned int> search_hashtable_strict(vector< bool * > list_bs, vector<unsigned int> leftside, vector<unsigned int> rightside, int side){
  //cout << "We show the bipartitions found by the search and corresponding trees" << endl;
  //keep in mind that hashtable and hash_lengths are global variables
  bool foundFlag = true;

  set<unsigned int> trees;
  bool lvalue = 0;
  bool rvalue = 0;
  
  for (unsigned int i = 0; i < list_bs.size(); i++){ //Each bipartition in the treeset


//      for (unsigned int l = 0; l < length_of_bitstrings[i]; l++){
//	  cout << list_bs[i][l];
//      }
//      cout << " --> [ ";
//      for (unsigned int l = 0; l < hash_lengths[i]; l++){
//        cout << hashtable[i][l] << " ";
//      }
//      cout << "]" << endl;

    foundFlag = true; //assume true and check each case for a contridiction to set it as false. If it passes the tests then it is true. 
	
	// if the first taxa on the left and the first taxa on the right side of the search bipart are on the same side of the bipart
	// we know we can bail out before checking anything else. But first we need to make sure there are taxa on each side of the | or we
	// get a segfault. 


	if (leftside.size() > 0)
	{
		if (list_bs[i][leftside[0]] > ::length_of_bitstrings[i]) 
		{
			lvalue = 0;
		}
		else
		{
			lvalue = list_bs[i][leftside[0]];
		}

		for (unsigned int j = 0; j < leftside.size(); j++)// Each taxa on the left side
		{ 
			//cout << length_of_bitstrings[i] <<" " << leftside[j] << " " << lvalue << endl;

			if(::length_of_bitstrings[i] <= leftside[j])
			{	
				if(lvalue == 1)
				{	
					//cout<< "if(length_of_bitstrings[i] < leftside[j]  && lvalue == 1)" << endl;
					foundFlag = false;
			  		break; //Break at any mistakes. 
				}
			}
			else if (list_bs[i][leftside[j]] != lvalue)
			{
				//cout<< "broke on if (list_bs[i][leftside[j]] != lvalue)" << endl;
				foundFlag = false;
		  		break; //Break at any mistakes. 
			}
		}

	}

	if (rightside.size() > 0)
	{
		if (list_bs[i][rightside[0]] > ::length_of_bitstrings[i]) 
		{
			rvalue = 0;
		}
		else
		{
			rvalue = list_bs[i][rightside[0]];
		}

		for (unsigned int j = 0; j < rightside.size(); j++)// Each taxa on the left side
		{ 
			//cout << length_of_bitstrings[i] <<" " << rightside[j] << " " << rvalue << endl;

			if(::length_of_bitstrings[i] <= rightside[j])
			{	
				if(rvalue == 1)
				{	
				//	cout<< "if(length_of_bitstrings[i] < rightside[j]  && rvalue == 1)" << endl;
					foundFlag = false;
			  		break; //Break at any mistakes. 
				}
			}
			else if (list_bs[i][rightside[j]] != rvalue)
			{
				//cout<< "broke on if (list_bs[i][rightside[j]] != rvalue)" << endl;
				foundFlag = false;
		  		break; //Break at any mistakes. 
			}
		}
	}

	if (leftside.size() > 0 && rightside.size() > 0 && lvalue == rvalue)
	{
		//cout<< "broke on if (leftside.size() > 0 && rightside.size() > 0 && lvalue == rvalue)" << endl;
		foundFlag = false;
		continue;	
	}

	if (side > 0)
	{
		unsigned int mycount = count (list_bs[i], list_bs[i]+::length_of_bitstrings[i], 1);
		if (side == 1 || side == 3)
		{

			//cout << foundFlag << " "  << side << " " <<  lvalue << " " << leftside.size() << " " << (NUM_TAXA-mycount) << " " << mycount << endl;
			if(lvalue == 1 && leftside.size() != mycount)
			{
				//cout << "Hit lvalue == 1 && leftside.size() != mycount)" << endl;
				foundFlag = false;
	  			continue;
			}
			else if(lvalue == 0 && leftside.size() != (NUM_TAXA-mycount))
			{
				//cout << "(lvalue == 0 && leftside.size() != (NUM_TAXA-mycount))" << endl;
				foundFlag = false;
	  			continue;
			}
		}
		
		if (side == 2 || side == 3)
		{
			if(rvalue == 1 && rightside.size() != mycount)
			{
				//cout << "(rvalue == 1 && rightside.size() != mycount)" << endl;
				foundFlag = false;
	  			continue;
			}
			else if(rvalue == 0 && rightside.size() != (NUM_TAXA-mycount))
			{
				//cout << "(rvalue == 0 && rightside.size() != (NUM_TAXA-mycount))" << endl;
				foundFlag = false;
	  			continue;
			}
		}

	}
	//If it is true I don't know if 0's mean that it's on the otherside of the bipartition or not in the tree at all. So check. 

	//if we got here and it's true it passed all the checks. 
	if (foundFlag == true)
	{
      //cout << "FOUND ONE" << endl;
      //for (unsigned int l = 0; l < ::length_of_bitstrings[i]; l++)
	  //{
	  //  cout << list_bs[i][l];
      //}
      //cout << " --> [ ";
     
      for (unsigned int l = 0; l < hash_lengths[i]; l++)
	  {
		if (HETERO == true)
		{
			bool taxaFlag = true;
			for (unsigned int j = 0; j < leftside.size(); j++)// Each taxa on the left side
			{
				if (taxa_in_trees[hashtable[i][l]][leftside[j]] == 0)
				{
					//cout << "broke at if (taxa_in_trees[l][leftside[j]] == 0)" << endl; cout << hashtable[i][l] << " "<< leftside[j] << endl;
					taxaFlag = false;
					break;
				}
			}
			for (unsigned int j = 0; j < rightside.size(); j++)// Each taxa on the left side
			{ 
				if (taxa_in_trees[hashtable[i][l]][rightside[j]] == 0)
				{
					//cout << "broke at if (taxa_in_trees[l][rightside[j]] == 0)"<< endl; cout << hashtable[i][l] << " " << rightside[j] << endl;
					taxaFlag = false;
					break;
				}
			}
			if (taxaFlag == true)
			{
				trees.insert(hashtable[i][l]);
			}
		}
		else
		{		
			trees.insert(hashtable[i][l]);
		}
        //cout << hashtable[i][l] << " "; Probably not the print statement I need. 
      }
      //cout << "]" << endl;
    }
	}
  return trees;
}




void LookUpLabels(vector<string> leftside, vector<string> rightside, vector<unsigned int> &lside, vector<unsigned int> &rside)
{
	for (int i = 0; i < (int)leftside.size(); i++){
    //cout << i << ltokens[i] << endl;
	   lside.push_back(::lm[leftside[i]]);
    //cout << ltokens[i] << endl;
    }
	
	for (int i = 0; i < (int)rightside.size(); i++){
    //cout << i << rtokens[i] << endl;
        rside.push_back(::lm[rightside[i]]);
    //cout << rtokens[i] << endl;
    }
}



vector<unsigned int> returnAllTrees(unsigned int numberOfTrees)
{	
	vector<unsigned int> allTrees;

	for(unsigned int i = 0; i < numberOfTrees; i++)
	{
		allTrees.push_back(i);
	}

	return allTrees;
}

set<unsigned int> taxaRequirements(vector<string> RequiredTaxa, vector<string> ExcludedTaxa)
{
    set<unsigned int> trees;
    set<unsigned int>::iterator it;
	
	vector<unsigned int> required;
	vector<unsigned int> excluded;

    //vector< bool* > taxa_in_trees;				// which taxa are in which trees					NEED
	//vector<unsigned int> length_of_bitstrings;

	LookUpLabels(RequiredTaxa, ExcludedTaxa, required, excluded);
	
    for (unsigned int i = 0; i < NUM_TREES; i++)
    {
		int flag = 0;

        for (unsigned int j = 0; j < required.size(); j++)
		{
			if(::taxa_in_trees[i][required[j]] == 0)
			{
				flag = 1;
				break;
			}
		} 
		
		for (unsigned int j = 0; j < excluded.size(); j++)
		{
			if(::taxa_in_trees[i][excluded[j]] == 1)
			{
				flag = 1;
				break;
			}
		} 

		if (flag == 0)
		{
			trees.insert(i);
		}
    }
    
    return trees;
}

set<unsigned int> getTrees(vector<string> leftside, vector<string> rightside, vector< bool * > &list_bs, int strictFlag)
{
	vector<unsigned int> lside;
	vector<unsigned int> rside;
	set<unsigned int> trees;
	set<unsigned int>::iterator it;
	
	LookUpLabels(leftside, rightside, lside, rside);
	
	//printBipartition(lside, rside);
	//printBipartition(leftside, rightside);

	trees = search_hashtable_strict(list_bs, lside, rside, strictFlag);
	
	//cout << "returned the following trees" << endl;
	//printSetOfTrees(trees);
	//cout << endl;
	//cout << endl;

	return trees;
}



