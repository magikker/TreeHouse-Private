grammar pql;

options
{
  language=C;
}

@parser::includes
{
  #include <iostream> 
  #include <string>
  #include <vector>
  #include <sstream>
  #include <set>
  #include <algorithm>
  #include <memory>	  
  
  //#include "HashTableSearch.h"
  #include "UtilityFunctions.h"
  #include "UserFunctions.h"
  #include "label-map.hh"
  #include "global.h"
  #include "pqlsymbol.h"
  #include "THGlobals.h"
  using namespace std;
}

@members 
{
   // Define some general things. 
   // vector< bool * > list_bs; // Move to global?
    //vector<unsigned int> allTrees; // Move to global?
     //std::vector<std::shared_ptr<pqlsymbol> > symbol_table;
}

//prog[vector< bool * > list_bs2] 
//  returns [pqlsymbol *result]
//  @init 
//  {
//      ::list_bs = $list_bs2; // Move to global?
//      ::all_trees = returnAllTrees(NUM_TREES);  // Move to global? 
//  }
//  : a=query+ {$result = $a.result;} // I need to handle things that ask for multiple commands in one run... That would need to be handled right here. 
//  ;
  
  
  prog
  
  : a=query+ // I need to handle things that ask for multiple commands in one run... That would need to be handled right here. 
  ;
  
 
query 
  : (a=assignment_expression 
    {
      pqlsymbol *  result = new pqlsymbol();
       
       if( $a.result->is_symbol() )
      {
          std::map<std::string, pqlsymbol * >::const_iterator iter = symbol_table.find( $a.result->get_symbolname() );
          if (iter == symbol_table.end())
          {
            // not found
            result = new pqlsymbol(ERROR, "symbol undefined, which shoudln't ever happen");
            delete $a.result;
          }
          else
          {
             //probably need a copy of the value here. 
            result = ( (*iter).second->get_value_as_new_pqlsymbol() );
            cout << $a.result->get_symbolname() << " = " << result->value_to_string() << endl;
            delete $a.result;
          }
          /*cout << "Query result " << result->value_to_string() <<  endl;*/
          query_results.push_back( result ); 
      }
      else
      {
         /*cout << "Query result " << $a.result->value_to_string() <<  endl;*/
         cout << $a.result->value_to_string() << endl;
         query_results.push_back( $a.result ); 
         delete result;
      } 
    } 
    )?
    NEWLINE
  ;
  
// E x p r e s s i o n s

//I  probably don't want to allow comma seperated expressions... 
//that stuff can get complicated, and I'm not trying to allow crazy c-like
//side effects

//expression
//  returns [pqlsymbol result]
//   : assignment_expression (',' assignment_expression)*
//    ;

assignment_expression
  returns [pqlsymbol * result]
    : id=IDENTIFIER '=' a=assignment_expression 
    {
     /*cout <<  pANTLR3_COMMON_TOKEN_to_string($id) << "=" << $a.result->value_to_string() << endl;*/
     /*cout << "setting symbol " << pANTLR3_COMMON_TOKEN_to_string($id) << " to " <<  $a.result->value_to_string() << endl;*/

      //At this point in the code I have two copies of the value. 
      //One is stored in the symbol table and the other returned up the interperter tree. 
      
       if( $a.result->is_symbol() )
      {
      
          std::map<std::string, pqlsymbol * >::const_iterator iter = symbol_table.find( $a.result->get_symbolname() );
          if (iter == symbol_table.end())
          {
            // not found
            //notFound = true;
            $result = new pqlsymbol(ERROR, "symbol undefined");
          }
          else
          {
             //probably need a copy of the value here. 
            $result = ( (*iter).second->get_value_as_new_pqlsymbol() );           
          }
      }
      
      else
      {
        //check if the symbol is a constant before we do any work. 
        
         std::map<std::string, bool >::const_iterator iter1 = constant_table.find( pANTLR3_COMMON_TOKEN_to_string($id) );
         if (iter1 == constant_table.end())
         {
          ::symbol_table[pANTLR3_COMMON_TOKEN_to_string($id)]= $a.result->get_value_as_new_pqlsymbol();
          $result = new pqlsymbol(SYMBOL, pANTLR3_COMMON_TOKEN_to_string($id) );
         }
        else{
        $result = new pqlsymbol(ERROR, "Can not assign a value to a static varible");
        }
       
        
      }
      delete $a.result;
      
    }
    |
      a = not_expression 
      {
        $result = $a.result;
      }
    ;

//AND, OR, and NOT should really just be functions. 

//logical_not_expression
//  returns [pqlsymbol *result]
//    : logical_and_expression ('||' logical_and_expression)*
//    ;



not_expression
  returns [pqlsymbol *result]
    @init
    {
      bool diffFlag = false;
      bool errorFlag = false;
      set<unsigned int> s1;
    }
       : ('!' { diffFlag = true; } )? a = difference_expression
         
    {
      if (diffFlag == false){
        $result = $a.result;
      }
      else{
        if( $a.result->is_treeset()  ){
          set<unsigned int> s1 = $a.result->get_treeset();
           
          set<unsigned int> sdiff;
         
          std::set_difference( all_trees.begin(), all_trees.end(), s1.begin(), s1.end(),  std::inserter( sdiff, sdiff.begin() ) );
          
          $result = new pqlsymbol(sdiff);
        }
        else{
          $result = new pqlsymbol(ERROR, "'Not' is only avaible for a treeset." );
        }
      }
    }
    ;

difference_expression
  returns [pqlsymbol *result]
    @init
    {
      bool diffFlag = false;
      bool errorFlag = false;
      set<unsigned int> s1;
    }
       : a = union_expression{ if( $a.result->is_treeset() ){s1 = $a.result->get_treeset();  }}
       ( '-' b = union_expression { diffFlag = true; }  )*
    {
      if (diffFlag == false){
        $result = $a.result;
      }
      else{
        if( $a.result->is_treeset() && $b.result->is_treeset()  ){
          set<unsigned int> s2 = $b.result->get_treeset();
           
          set<unsigned int> diff;
          
          std::set_difference( s1.begin(), s1.end(), s2.begin(), s2.end(),
    	std::inserter( diff, diff.begin() ) );
           s1 = diff;
          $result = new pqlsymbol(diff);
        }
        else{
          $result = new pqlsymbol(ERROR, "'Difference' is currently only avaible between two int vector types." );
        }
      }
    }
    ;

union_expression
  returns [pqlsymbol *result]
    @init
    {
      bool orFlag = false;
      bool errorFlag = false;
      set<unsigned int> s1;
    }
    : a = intersection_expression{    
      if( $a.result->is_treeset() ){s1 = $a.result->get_treeset();  }
      }  ( '+' b = intersection_expression 
    {
     orFlag = true; 
     if(  $b.result->is_treeset()  ){
          /*set<int> s1 = $a.result->get_treeset();*/
          set<unsigned int> s2 = $b.result->get_treeset();
          set<unsigned int> uni;
          std::set_union( s1.begin(), s1.end(), s2.begin(), s2.end(),
    	std::inserter( uni, uni.begin() ) );
           s1 = uni;
           
          $result = new pqlsymbol(uni);
     }
     else{errorFlag = true;}
        
    }  )*
    {
      if (orFlag == false){
        $result = $a.result;
      }
      else if (errorFlag){
          $result = new pqlsymbol(ERROR, "'Union' is currently only avaible between two int vector types." );
      }
    }
    ;

intersection_expression
  returns [pqlsymbol *result]
    @init
    {
      bool andFlag = false;
      bool errorFlag = false;
      set<unsigned int> s1;
    }
    : a = equality_expression {if( $a.result->is_treeset() ){s1 = $a.result->get_treeset();  }}
     ('^' b = equality_expression
    {
          andFlag = true;
          if( $a.result->is_treeset() && $b.result->is_treeset()  ){
          
          set<unsigned int> s2 = $b.result->get_treeset();
           
          set<unsigned int> intersection;
          
          std::set_intersection( s1.begin(), s1.end(), s2.begin(), s2.end(),
    	std::inserter( intersection, intersection.begin() ) );
          
          s1 = intersection;
                    
          $result = new pqlsymbol(intersection);
        }
        else{
          errorFlag = true;
        }
        
    })*
    {
      if (andFlag == false){
        $result = $a.result;
      }
      else if (errorFlag){
        $result = new pqlsymbol(ERROR, "'intersection' is currently only avaible between two treesets." );
      }
    }
   
    ;

equality_expression
  returns [pqlsymbol * result]
  @init
  {
    bool equalityFlag = false;
    bool eqFlag = true; 
  }
  
    : a = postfix_expression (('==' {eqFlag = true; equalityFlag = true;}|'!=' {eqFlag = false; equalityFlag = true;}) b = postfix_expression)*
      {
        if (equalityFlag)
        {
           if(eqFlag)
           {
             if( $a.result == $b.result )
             {
               $result = new  pqlsymbol( true);
             }
             else
             {
                $result = new pqlsymbol( false);
             }
           }
           else
           {
              if( $a.result != $b.result )
              {
                 $result = new pqlsymbol( true);
              }
              else
              {
                $result = new pqlsymbol( false);
              }
           }
         }
         else
         {
           $result = $a.result;
         }
      }
    ;
        
postfix_expression
  returns [pqlsymbol * result]
  @init 
  {
      bool arrayorfunct = false;
  }
    :  a = primary_expression {
    /*cout << "In postfix_expression a = " <<  $a.result->value_to_string() << endl; */
    }
    (  '[' b = assignment_expression ']'  //I still have to work this part out!
      {
         // you didn't apply [] to a vector
         if (! $a.result->is_vect() )
         {
            //cout << "Oh, it's not an vect that you are indexing.. that's really bad" << endl;
            $result = new pqlsymbol(ERROR, "Tried to index a vector without using an int");
         }
         
         //if you don't give [] an int
         else if (! $b.result->is_int() ) 
         {
            //cout << "Oh, it's not an Int.. that's really bad" << endl;
            $result = new pqlsymbol(ERROR, "Tried to index a vector without using an int");
         }
          
         // if you give [] an int, but out of range 
         else if (! ($a.result->get_vect_size() >= $b.result->get_int() ) )
         {
           //cout << "index out of range for vector" << endl;
           $result = new pqlsymbol(ERROR, "index out of range");
         }
          
         else if (! ( $b.result->get_int() >= 0) )
         {
           //cout << "index out of range for vector, we don't support negative numbers. " << endl;
           $result = new pqlsymbol(ERROR, "index out of range");
         } 
         // you gave it all correct things, so now lets return your symbol. 
         else
         {
            $result = $a.result->get_list_symbol($b.result->get_int() );
         }
         arrayorfunct = true;
         
      }
      
    | '(' ')' 
      {
        arrayorfunct = true;
        /*cout << "testing if " << $a.result->get_string()<< " is a function" <<endl;*/
        std::map<std::string, vfptr>::const_iterator iter = voidFunctMap.find($a.result->get_string());
        if (iter == voidFunctMap.end())
        {
        // not found
          $result = new pqlsymbol(ERROR, "not a legal function");
          //cout << "not a legal function" << endl;
        }
        
        else
        {
          $result = (*iter).second();
        }
        
       // delete $a.result;  
          
      }
      
    | '(' c=argument_expression_list ')'
      {
        arrayorfunct = true;
        /*cout << "testing if " << $a.result->get_string()<< " is a function" <<endl;*/
        std::map<std::string, afptr>::const_iterator iter = argFunctMap.find($a.result->get_string());
        if (iter == argFunctMap.end())
        {
        // not found         
        //cout << "not a legal function" << endl;
        $result = new pqlsymbol(ERROR, "not a legal function");
        }
        else
        {
        //Extra if statement and iterator to check if the template for type checking is used
        std::map<std::string, vector <vector < dataType> > >::iterator it = argMap.find($a.result->get_string());
        if(it==argMap.end()){
          	//Gotta have a value copy thing here. 
          	$result =  (*iter).second($c.result) ;
          	while(!$c.result.empty()){
                	  delete $c.result.back();
                	  $c.result.pop_back();
                	  }
          }
          else{
          	$result = u_template($c.result, $a.result->get_string());
          	}
        }
    }
    
    
    )*
    {
      if (arrayorfunct == false)
      {
        $result = $a.result;
      }
    }
    ;
			
argument_expression_list
  returns [vector<pqlsymbol * > result]
    :   a = assignment_expression 
    {
      result.push_back($a.result); 
    }
     
    (',' b = assignment_expression 
      {
        result.push_back($b.result); 
      }
    )*
    ;

primary_expression 
   returns [pqlsymbol * result]
    : id = IDENTIFIER 
    {
          std::map<std::string, pqlsymbol * >::const_iterator iter1 = symbol_table.find( pANTLR3_COMMON_TOKEN_string_lit_to_string( $id)  );
          if (iter1 == symbol_table.end())
         {
            // not found
            // TODO: I should check if it's a function name, and if so pass up a function type.
            $result = new pqlsymbol(pANTLR3_COMMON_TOKEN_string_lit_to_string( $id) );
            //$result = new pqlsymbol(ERROR, "Varible undefined.");
         }
         else
         {
             //probably need a copy of the value here. 
             $result = ( (*iter1).second->get_value_as_new_pqlsymbol() );
          }
      

        /*
        std::map<std::string, pqlsymbol * >::const_iterator iter = symbol_table.find(pANTLR3_COMMON_TOKEN_to_string($id));
        if (iter == symbol_table.end())
        {
        // not found
          $result = new pqlsymbol(ERROR, "symbol undefined");
         }
        else
        {
        //probably need a copy of the value here. 
          $result = ( (*iter).second->get_value_as_new_pqlsymbol() );
        }
       */
    }
    
    
    | a = atom {$result = $a.result;}
    | '(' a = assignment_expression {$result = $a.result;} ')'
    ;
    
			
    
    
    
    
    
//look
constant_expression_list
  returns [pqlsymbol * result]
   @init 
  {
      bool errored = false;
      vector <pqlsymbol * > thelist;
  }
    :   (
          a = assignment_expression { thelist.push_back($a.result); }  
         | 
            ra = RANGE_LITERAL{ 
                vector<int> range = pANTLR3_COMMON_TOKEN_to_intvect($ra);
                for(int i = 0; i < range.size(); i++ ){
                  thelist.push_back(new pqlsymbol(range[i]));
                } 
            }   
        )

        (',' (b = assignment_expression{
                  if ( $b.result->get_object_type() == thelist[0]->get_object_type() && $b.result->get_data_type() == thelist[0]->get_data_type() )
                  {
                       thelist.push_back($b.result); 
                  } 
                  else{
                  	//cout << "HIT AN ERROR IN A [LIST]" << endl;
                  	errored = true; 
                  }
                 
             } 
             | 
                rb = RANGE_LITERAL{ 
                vector<int> range = pANTLR3_COMMON_TOKEN_to_intvect($rb);
                for(int i = 0; i < range.size(); i++ ){
                  thelist.push_back(new pqlsymbol(range[i]));
                } 
            }
             )
             
       )*
       {
         if(! errored)
         {
           $result = new pqlsymbol( thelist );
           /*cout << "did not hit an error" << endl;*/
           /*cout << "The value of the expression list = " <<$result->value_to_string() << endl;*/
         }
         else 
           $result = new pqlsymbol(ERROR, "type mismatch in list");
       }
    ;
    
atom
returns [pqlsymbol * result]
  @init 
  {
      bool emptylist = true;
      bool errored = false;
  }
: a = constant {
/*cout << "DECIMAL_LITERAL = " << $a.result->value_to_string() << endl;*/
$result = $a.result;
/*cout << "DECIMAL_LITERAL = " << $result->value_to_string() << endl;*/

} 
|'[' 
  (
      b = constant_expression_list 
      {
        $result = $b.result;
        emptylist = false; 
      } 
  )? 
 ']'
 {
   if (emptylist == true)
   {
     $result = new pqlsymbol();
   }
 }
 | '{' (
 b = constant_expression_list
{
        //set<int> *tempset  = new set<int>;
        vector<int> tempvect = $b.result->get_int_vect();
        std::set<unsigned int> tempset( tempvect.begin(), tempvect.end() );
        $result = new pqlsymbol(tempset);
         emptylist = false;
}
)?  '}'
{
   if (emptylist == true)
   {
     set<unsigned int> tempset;
     $result = new pqlsymbol(tempset);
   }
}
 
 
;


constant
    returns [pqlsymbol * result]
      :   ( a = DECIMAL_LITERAL {
    	/*cout << "DECIMAL_LITERAL = " << pANTLR3_COMMON_TOKEN_to_int($a) << endl;*/
    	$result =  new pqlsymbol( pANTLR3_COMMON_TOKEN_to_int($a)  );
    } )
    
    |  ( a = CHARACTER_LITERAL {
    	/*cout <<"CHARACTER_LITERAL = " <<  pANTLR3_COMMON_TOKEN_to_string($a) << endl;*/
    	$result = new pqlsymbol( pANTLR3_COMMON_TOKEN_to_string($a) );
    } ) 
    
    |  ( a = STRING_LITERAL {
    	/*cout << "STRING_LITERAL = " << pANTLR3_COMMON_TOKEN_to_string($a) << endl;*/
    	$result =  new pqlsymbol( pANTLR3_COMMON_TOKEN_string_lit_to_string( $a)  );
    } ) 
    
    |  ( a = FLOATING_POINT_LITERAL {
    	/*cout << "FLOATING_POINT_LITERAL = " << pANTLR3_COMMON_TOKEN_to_double($a) << endl;*/
    	 $result =  new pqlsymbol( pANTLR3_COMMON_TOKEN_to_double($a) );
    } ) 
    ;
   
IDENTIFIER
    :  LETTER (LETTER|'0'..'9')*
    ;
	
fragment
LETTER
    :  '$'
    |  'A'..'Z'
    |  'a'..'z'
    |  '_'
    ;

CHARACTER_LITERAL
    :   '\'' ( EscapeSequence | ~('\''|'\\') ) '\''
    ;

STRING_LITERAL
    :  '"' STRING_GUTS '"'
    ;
    
fragment
STRING_GUTS 
    :  ( EscapeSequence | ~('\\'|'"') )* 
    ;
    
RANGE_LITERAL
    : ('-')? ('0'..'9')+ '..' ('-')?('0'..'9')+ 
    ;

DECIMAL_LITERAL 
    : ('0' | ('-')? '1'..'9' '0'..'9'*)  
    ;

FLOATING_POINT_LITERAL
    :   ('-')?('0'..'9')+ '.' ('0'..'9')* 
    |   '.' ('0'..'9')+ 
    //|  ('-')? ('0'..'9')+ 
    ;

fragment
EscapeSequence
    :   '\\' ('b'|'t'|'n'|'f'|'r'|'\"'|'\''|'\\')
    ;

COMMENT
    :   '/*' ( options {greedy=false;} : . )* '*/' {$channel=HIDDEN;}
    ;

LINE_COMMENT
    : '//' ~('\n'|'\r')* {$channel=HIDDEN;}
    ;

NEWLINE: ('\r\n'|'\n'|'\r');

WS  :  (' '|'\r'|'\t'|'\u000C'|'\n') {$channel=HIDDEN;}
    ;
