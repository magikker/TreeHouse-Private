//****************************************************************/
/*This is HashSim, a fast topological similarity matrix algorithm. It is based on 
Fast HashRF, which is in turn based on HashRF, by SeungJin Sul.

(c) 2010 HashSim, Fast HashRF: Suzanne Matthews
(c) 2009 HashRF: SeungJin Sul
*/
/*****************************************************/

#include <sys/time.h>
#include <sys/resource.h>
#include <cassert>
#include <valarray>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <string>
#include <cmath>
#include <cfloat>

#include "compute.hh"
#include "global.h"

using namespace std;

bool intcomp (unsigned int i,unsigned int j) { return (i<j); }

int max(int x, int y){ 
  if (x>y)
    return x;
  else
    return y;
} 

void calculate(int a, int b, int d, int distance_measure, unsigned int &value, float &value2){
  int sig1, sig2;
  float temp, temp2;
  //cout << "----" << endl;
  value = 0;
  value2 = 0;
    switch (distance_measure) {
    case 0: 
    case 24: 
      value = b; //RF
      break;
    case 1:
      value = b+b; //Hamming
      break;
    case 2: 
      value2 =(float)a/(a+2*b); //Jaccard/Tanimoto
      break;
    case 3:
      value2 = (float)a/(a+b); //Dice
      break;
    case 4:
      value2 = (float)a/NUMBIPART; //Russel/Rao
      break;
    case 5:
      value2 = (float)a/(a+4*b); //Sokal/Sneath
      break;
    case 6:
      d = NUMBIPART - (a+b+b);
      value2 = (float)(a+d-2*b)/NUMBIPART; //Hamann
      break;
    case 7: 
      d = NUMBIPART - (a+b+b);
      value2 = (float)(a*d - b*b)/(a*d + b*b); //Yule
      break;
    case 8:
      value2 = (float)a/(sqrt((a+b)*(a+b))); //Ochai/Cosine
      break;
    case 9:
      d = NUMBIPART - (a+b+b);
      sig1 = max(a,b)+max(b,d)+max(a,b)+max(b,d);
      sig2 = max(a+b,b+d) + max(a+b,b+d);
      value2 = (float)(sig1-sig2)/(2*NUMBIPART); //Anderberg
      break;
    case 10: 
      value2 = (float)(NUMBIPART*a)/((a+b)*(a+b)); //Forbes
      break;
    case 11:
      temp = (float)(0.5*(a*b + a*b) + b*b); //Mountford
      if (temp == 0)
	value2 = 1;
      else
	value2 = (float)a/(0.5*(a*b + a*b) + b*b); //Mountford
      break;
    case 12:
      value2 = (float)(pow((float)a,2))/((a+b)*(a+b)); //Sorgenfrei
      break;
    case 13:
      temp = (float)(a+b)/NUMBIPART;
      if (a != 0)
	value2 = (float)log10(a) - log10(NUMBIPART) - log10(temp) - log10(temp); //Gilbert/Wells
      else 
	value2 = FLT_MIN -  (float)log10(NUMBIPART) - log10(temp) - log10(temp); //Gilbert/Wells
      break;
    case 14:
      temp = (float) (NUMBIPART*a);
      temp = temp - ((a+b)*(a+b));
      temp2 =(float) (NUMBIPART*a + (a+b)*(a+b)); //Tarwid
      value2 = (float)temp/temp2;
      break;
    case 15:
      value2 = (float)a/(a+b) + (float)a/(a+b); //Johnson
      break;
    case 16:
      d = NUMBIPART - (a+b+b);
      temp = (float)(a*b + 2*b*b + b*d); //Peirce
      if (temp == 0)
	value2 = 1;
      else
	value2 = (float)(a*b + b*b)/(a*b + 2*b*b + b*d); //Peirce
      break;
    case 17:
      temp = (float)1/(a+b);
      value2 = (float)(a/2)*(temp + temp); //driver/kroeber
      break;
    case 18:
      value2 = (float)a/(sqrt((a+b)*(a+b))) - (float)max(a+b,a+b)/2; //fager/mcgowan 
      break;
    case 19:
      value2 = (float)sqrt(b+b); //euclidean
      break;
    case 20:
      d = NUMBIPART - (a+b+b);
      temp2 = (float)fabs(a*d - b*b);
      temp2 -= (float)(NUMBIPART/2);
      temp = (float)pow(temp2, 2);
      temp = (float)NUMBIPART*temp;
      temp2 = (a+b)*(a+b); 
      temp2 *= (b+d)*(b+d);
      temp2 = (float)temp/temp2;
      value2 = (float)log10(temp2); //Stiles
      break;
    case 21:
      d = NUMBIPART - (a+b+b);
      temp = a+d;
      value2 = (float)temp/(2*(b+b)+a+d); //Rogers/Tanimoto
      break;
    case 22:
      d = NUMBIPART - (a+b+b);
      temp2 = (float)NUMBIPART*a;
      temp2 -= ((a+b)*(a+b));
      temp = (float)(NUMBIPART*NUMBIPART)*temp2;
      temp2 =(float)(a+b)*(a+b)*(b+d)*(b+d);
      value2 = (float)temp/temp2; //Eyraud
      break;
    case 23:
      temp = b;
      temp2 = b+a;
      value2 = (float)temp/temp2; //weighted RF
      break;
    case 25:
      value2 = (float)(2*b)/NUMBIPART; //mean manhattan
      break;
    case 26:
      d = NUMBIPART - (a+b+b);
      value2 = (float)(a+d)/NUMBIPART; //simple matching
      break;
    case 27:
      value2 = (float)(a*d - b*b)/(sqrt((NUMBIPART*(a+b)*(a+b)))); //Dennis
      break;
    case 28:
      d = NUMBIPART - (a+b+b);
      sig1 = max(a,b)+max(b,d)+max(a,b)+max(b,d);
      sig2 = max(a+b,b+d) + max(a+b,b+d);
      value2 = (float)(sig1-sig2)/(2*NUMBIPART - sig2); //Goodman/Kruskal
      break;
    case 29:
      value2 = (float)(NUMBIPART*pow((a-0.5),2))/((a+b)*(a+b)); //Fossum
      break;
    case 30:
      value2 = (float)(b+b)/(2*a + b + b); //bray/curtis
      break;
    case 31:
      d = NUMBIPART - (a+b+b);
      temp = (float)(fabs(a*(b+d)));
      temp2 = (float)(fabs(b*(a+b)));
      if (temp2 == 0)
	value2 = NAN;
      else
	value2 = (float)temp/temp2; //AMPLE
      break;
    case 32:
      d = NUMBIPART - (a+b+b);
      temp = (float)(sqrt(a*d) - b);
      temp2 = (float)(sqrt(a*d) + b);
      value2 = (float)temp/temp2; //Yule-W
      break;
    case 33:
      value2 = (float)(a - (b+b))/NUMBIPART; //FAITH
      break;
    case 34:
      value2 = (float)(a - (b+b))/(a+b+b); //FAITH-2
      break;
    default:
      cerr << "Invalid option!" << endl;
      exit(0);
    }
}

void calculatem(int a, int b, int c, int d, int distance_measure, float & value){ //for multifurcating
  int sig1,sig2;
  float temp, temp2;

  //cout << "distance measure is: " << distance_measure << endl;
  //cout << "b is: " << b << endl;
  //cout << "c is: " << c << endl;
  //cout << "d is: " << d << endl;
  //cout << "a is: " << a << endl;
    switch (distance_measure) {
    case 0: 
    case 24: 
      value = (float)(b+c)/2; //RF
      //cout << "value is: " << value << endl;
      break;
    case 1: 
      value =  b+c; //Hamming
      break;
    case 2: 
      value =  (float)a/(a+b+c); //Jaccard/Tanimoto
      break;
    case 3: 
      value = (float)(2*a)/(2*a + b + c); //Dice
      break;
    case 4: 
      value = (float)a/NUMBIPART; //Russel/Rao
      break;
    case 5: 
      value = (float)a/(a + 2*b + 2*c); //Sokal/Sneath
      break; 
    case 6:
      d = NUMBIPART - (a+b+c);
      value = (float)(a+d-b-c)/NUMBIPART; //Hamann
      break;
    case 7:
      d = NUMBIPART - (a+b+c);
      value = (float)(a*d - b*c)/(a*d + b*c); //Yule
      break;
    case 8:
      value = (float)a/(sqrt((a+b)*(a+c))); //Ochai/Cosine
      break;
    case 9:
      d = NUMBIPART - (a+b+c);
      sig1 = max(a,b)+max(c,d)+max(a,c)+max(b,d);
      sig2 = max(a+c,b+d) + max(a+b,c+d);     
      value = (float)(sig1-sig2)/(2*NUMBIPART); //Anderberg
      break;
    case 10: 
      value = (float)(NUMBIPART*a)/((a+b)*(a+c)); //Forbes
      break;
    case 11:
      temp = (float)(0.5*(a*b + a*c)) + b*c;
      if (temp == 0)
	value = 1;
      else
      value = (float)a/(0.5*(a*b + a*c) + b*c); //Mountford
      break;
    case 12:
      value = (float)(pow((float)a,2))/((a+b)*(a+c)); //Sorgenfrei
      break;
    case 13:
      temp = (float)(a+b)/NUMBIPART;
      temp2 = (float)(a+c)/NUMBIPART;
      if (a!= 0)
	value = (float)log10(a) - log10(NUMBIPART) - log10(temp) - log10(temp2); //gilbert/wells
      else
	value = FLT_MIN - (float)log10(NUMBIPART) - log10(temp) - log10(temp2); //gilbert/wells
      break;
    case 14:
      temp = (float)(NUMBIPART*a);
      temp = temp - ((a+b)*(a+c));
      temp2 = (float)(NUMBIPART*a + (a+b)*(a+c)); 
      value = (float)temp/temp2; //Tarwid
      break;
    case 15:
      value = (float)(a/(a+b)) + (float)(a/(a+c)); //Johnson
      break;
    case 16:
      d = NUMBIPART - (a+b+c);
      temp = (float)(a*b + 2*b*c + c*d);
      if (temp == 0)
	value = 1;
      else
	value = (float)(a*b + b*c)/(a*b + 2*b*c + c*d); //Peirce
      break;
    case 17:
      temp = (float)1/(a+b);
      temp2 = (float)1/(a+c);
      value = (float)(a/2)*(temp + temp2); //Driver/Kroeber
      break;
    case 18:
      value = (float)a/(sqrt((a+b)*(a+c))) - (float)max(a+b,a+c)/2; //Fager/McGowan 
      break;
    case 19:
      value = (float)sqrt(b+c); //Euclid
      break;
    case 20:
      d = NUMBIPART - (a+b+c);
      temp = (float)NUMBIPART*pow((float)(fabs(a*d - b*c)-(NUMBIPART/2)),2);
      temp2 = (a+b)*(a+c);
      temp2*= (b+d)*(c+d);
      temp2 = (float)temp/temp2;
      value = (float)log10(temp2); //Stiles
      break;
    case 21:
      d = NUMBIPART - (a+b+c);
      temp = a+d;
      value = (float)temp/(2*(b+c)+a+d); //Rogers/Tanimoto
      break;
    case 22:
      d = NUMBIPART - (a+b+c);
      temp2 = (float)NUMBIPART*a;
      temp2 -= (float)((a+b)*(a+c));
      temp = (float)(NUMBIPART*NUMBIPART)*temp2;
      temp2 =(float)(a+b)*(a+c)*(b+d)*(c+d);
      value = (float)temp/temp2; //Eyraud
      break;
    case 23:
      temp = (float)(b+c);
      temp2 = (float)(2*a + b + c);
      value = (float)temp/temp2; //RF rate
      break;
    case 25:  
      value = (float)(b+c)/NUMBIPART; //Mean Manhattan
      break;
    case 26:
      d = NUMBIPART - (a+b+c);
      value = (float)(a+d)/NUMBIPART; //Simple Matching
      break;
    case 27:
      d = NUMBIPART - (a+b+c);
      value = (float)(a*d - b*c)/(sqrt((NUMBIPART*(a+b)*(a+c)))); //Dennis
      break;
    case 28:
      d = NUMBIPART - (a+b+c);
      sig1 = max(a,b)+max(c,d)+max(a,c)+max(b,d);
      sig2 = max(a+c,b+d) + max(a+b,c+d);     
      value = (float)(sig1-sig2)/(2*NUMBIPART - sig2); //Goodman/Kruskal
      break;
    case 29:
      value = (float)(NUMBIPART*pow((a-0.5),2))/((a+b)*(a+c)); //Fossum
      break;
    case 30:
      value = (float)(b+c)/(2*a + b + c); //Bray/Curtis
      break;
    case 31:
      d = NUMBIPART - (a+b+c);
      temp = (float)(fabs(a*(c+d)));
      temp2 = (float)(fabs(c*(a+b)));
      if (temp2 == 0)
	value = NAN;
      else
	value = (float)temp/temp2; //AMPLE
      break;
    case 32:
      d = NUMBIPART - (a+b+c);
      temp = (float)(sqrt(a*d) - sqrt(b*c));
      temp2 = (float)(sqrt(a*d) + sqrt(b*c));
      value = (float)temp/temp2; //Yule-W
      break;
    case 33:
      value = (float)(a - (b+c))/NUMBIPART; //FAITH
      break;
    case 34:
      value = (float)(a - (b+c))/(a+b+c); //FAITH-2
      break;
    default: 
      cerr << "Invalid option!" << endl;
      exit(0);
    }
}

void calculatew(float main, float aux1, float aux2, int distance_measure, float &value){

  float temp, temp2;
  value = 0.0;
  temp = 0.0;
  temp2 = 0.0;
  switch (distance_measure) {
  case 2: 
      value =(float)main/(aux1+ aux2 - main); //Jaccard/Tanimoto
      break;
  case 3:
    value = (float)(2*main)/(aux1+aux2); //Dice
    break;
  case 4:
    value = (float)main/NUMBIPART; //Russel/Rao
    break;
  case 5:
    value = (float)main/((2*aux1) + (2*aux2) -(3*main)); //Sokal/Sneath
    break;
  case 8:
    value = (float)main/(sqrt(aux1*aux2)); //Ochai/Cosine
    break;
  case 10: 
    value = (float)(NUMBIPART*main)/((aux1)*(aux2)); //Forbes
    break;
  case 12:
    value = (float)(pow(main,2))/(aux1*aux2); //Sorgenfrei
    break;
  case 14:
    temp = (float) (NUMBIPART*main);
    temp = temp - (aux1*aux2);
    temp2 =(float) (NUMBIPART*main + aux1*aux2); //Tarwid
    value = (float)temp/temp2;
    break;
  case 15:
    value = (float)main/aux1 + (float)main/aux2; //Johnson
    break;
  case 17:
    temp = (float)1/aux1;
    temp2 = (float)1/aux2;
    value = (float)(main/2)*(temp + temp2); //driver/kroeber
    break;
  case 18:
    value = (float)main/(sqrt(aux1*aux2)) - (float)max(aux1,aux2)/2; //fager/mcgowan 
    break;
  case 23:
    temp = main;
    temp2 = (float)aux1+(float)aux2;
    value = (float)temp/temp2; //weighted RF
    break;
  case 24:
    value = (float)main/2;
    break;
  default:
    cerr << "Invalid option!" << endl;
    exit(0);
  }
}
