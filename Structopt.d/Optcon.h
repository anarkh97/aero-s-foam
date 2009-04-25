#ifndef _OPTCON_H_
#define _OPTCON_H_

#ifdef STRUCTOPT

#include <Utils.d/resize_array.h>

class Optopr;

//------------------------------------------------------------------------------
//                      Equality and Inequality Constraints
//                                                      created  9/3/98 by Kurt
//------------------------------------------------------------------------------

class Optcon {

   public:

     int numcon;                          //Number of Constraints
     int numeqc;                          //Number of Equality Constraints
     int numieq;                          //Number of Inequality Constraints

     ResizeArray<double> valcon;          //Value of Constraint
     ResizeArray<double> gradcon;         //Gradient of Constraint
     ResizeArray<double> scalcon;         //Scaling Factor 
     ResizeArray<double> typcon;          //Type of Constraints 

     ResizeArray<double> savecon;         //Copy of Constraint Values

     ResizeArray<Optopr*>   opr;          //Functions of Constraint Critiria
          
   public:         
		      
     Optcon();

     void func();
     void grad();
     void print(FILE*);
     void printres(FILE*);

     double getScaledValue(int icon)    { return valcon[icon]*scalcon[icon]; }
     double getScaledGradient(int icon) { return gradcon[icon]*scalcon[icon]; }
     
     void removeConstraint(int);

     void saveValue();
     void restoreValue();
};

#endif

#endif
