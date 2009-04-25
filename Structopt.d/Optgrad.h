#ifndef _OPTGRAD_H_
#define _OPTGRAD_H_

#ifdef STRUCTOPT

#include <cstdio>
#include <Utils.d/resize_array.h>

class Optsol;
class Optpro;
class Structopt;
class Optfilter;

//------------------------------------------------------------------------------

class Optgrad {

     public:

        Optpro    *optpro;               //Current Optimization Problem
        Optsol    *optsol;               //Current Optimization Solver
        Structopt *structopt;            //Current Structural Opt.Prob.
        Optfilter *filterOp;             //Filteroperator for gradients

        int typ;                         //Type of Gradient Method
        int anatyp;                      //Dircet of Adjoint Method     
        int numtyp;                      //Forward or Central Difference Scheme
        int epstyp;                      //Relative or Absolute Increment
        int icount;                      //Counter for Gradient Evaluations
        int numFilCrit;                  //number of Criteria to be projected

        int* actCriteria;                //Flag for active criteria
        int* filterCrit;                 //List of Criteria to be projected
     
        double epsval;                   //Size of Distrubation
          
     public:
     
	void print(FILE*);	
	
	void setactpro ( Optpro *_optpro=0, Optsol *_optsol=0, 
	                 Structopt *_structopt=0 );
	
	void buildgrad ( graddata& );

        void grad         (int* active=0);

	void gradanalytic (int* active=0);
        void gradanalytic (int*,double**,int,int);

	void gradnumeric  ();
};

#endif

#endif
