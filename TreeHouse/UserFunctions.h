
#ifndef _USER_FUNCTIONS_HH_
#define _USER_FUNCTIONS_HH_

//Global Vars like the hashtable.
#include "THGlobals.h"

//The label map class. 
#include "label-map.hh"

//Function Files 
#include "AnalysisFunctions.h"
#include "SearchFunctions.h"
#include "UtilityFunctions.h"
#include "InterfaceFunctions.h"
#include "quartet.h"
#include "ConsensusFunctions.h"
#include "clustering.h"
#include "distance.h"
#include "TreeMeasures.h"
#include "NewickParser.h"



#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/variant.hpp>

pqlsymbol * u_template(vector<pqlsymbol * > arglist, string functName);

void init_the_functs();

#endif
