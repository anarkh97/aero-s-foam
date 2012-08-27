#ifndef _OPTPRO_H_
#define _OPTPRO_H_

#ifdef STRUCTOPT

#include <Utils.d/resize_array.h>

class Optcrit;
class Optvar;
class Optopr;
class Optcon;
class Optobj;
class Optsol;
class Structopt;
class Relsol;

class criterion;
class Optopr;

struct critdata;
struct funcall;
struct absvardata;

//------------------------------------------------------------------------------
//                   Optimization Problem 
//                                          created  9/1/98 by Kurt
//------------------------------------------------------------------------------

class Optpro {

   public:

     int type;                         //Type of Problem: 0 - Design Optimization
                                       //                 1 - Reliability Analysis
   
     int numcrit;                      //Number of Optimization Criteria

     FILE* optunitout;                 //Output file

     int** OCCtable;                   //Objective/Constraint - Criteria table
   
     ResizeArray<Optcrit*>  opc;       //Optimization Criteria
		   
     Optobj *optobj;                   //Objective
          
     Optcon *optcon;                   //Constraints - Failure Criteria
     
     Optvar *optvar;                   //Optimization Variables

     Optsol *optsol;                   //Solution Strategy for Design Optimization
 
     Relsol *relsol;                   //Solution Strategy for Reliability Optimization

     int numElecStcVar;
     int numThermStcVar;

   public:         
		      
     Optpro(int);

     // build optimization problem

     void addCriteria   ( criterion & crit );

     void removeCriteria( int );

     void addObjective   ( double scal, funcall & func );     

     void addObjective   ( double scal, Optopr* );     

     void removeObjective( );

     void removeConstraint(int);

     void addConstraint ( int num, int typ, double scl, funcall & func );

     void addFailureCriteria  ( int num, funcall & func, int, double );

     void addDsgvar     ( absvardata & var );

     void addRndvar     ( absvardata & var );

     void addStcvar     ( int num, int typ, int loc1, int loc2,
                          funcall & func, int rfg=0 );

     void saveProblem     (double&,double*,int);

     void restoreProblem  (double&,double*,int);

     Optopr * addFunc   ( funcall & func );

     void inputcheck();

     // administration

     int checkAnalyticalSA();

     void setOutput(FILE* fp=0);
     void print(int hflag=1);

     // solution

     void solve(Structopt*);

     void buildOCC();
     void getActiveCriteria(int*,int*);


     // reliability specific functions

     int  checkFailcritOnly();

     void generateRandom();
     void computeRelDSA(double**,Optcrit**,int);

};

#endif

#endif
