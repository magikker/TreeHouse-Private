//*****************************************************************/
/*This is HashSim, a fast topological similarity matrix algorithm. It is based on 
Fast HashRF, which is in turn based on HashRF, by SeungJin Sul.

(c) 2010 HashSim, Fast HashRF: Suzanne Matthews
(c) 2009 HashRF: SeungJin Sul
*/
/*****************************************************/
#ifndef _COMPUTE_HH
#define _COMPUTE_HH

#include <iostream>
using namespace std;

//for unweighted, binary trees
void calculate(int a, int b, int d, int distance_measure, unsigned int &value, float &value2);
//for unweighted, multifurcating trees
void calculatem(int a, int b, int c, int d, int distance_measure, float & value);
//for weighted trees
void calculatew(float a, float ab, float ac, int distance_measure, float &value);
#endif
