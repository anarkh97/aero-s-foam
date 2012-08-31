#ifndef _OPTOBJ_H_
#define _OPTOBJ_H_

#ifdef STRUCTOPT

#include <Utils.d/resize_array.h>

class Optopr;

//------------------------------------------------------------------------------
//                             Objective
//                                                      created  9/1/98 by Kurt
//------------------------------------------------------------------------------

class Optobj {

   public:
		        
     int numobj;                            //Number of Components

     double valobj;                         //Value of Objective
     double gradobj;                        //Gradient of Objective
     double scalobj;                        //Scaling Factor 

     double saveobj;                        //Copy of Objective Value
     
     Optopr *multicrit;                     //Multicrit. Problem
          
   public:         
		      
     Optobj();

     void func();
     void grad();
     void print(FILE*);
     void printres(FILE*);

     double getScaledValue()    { return valobj*scalobj; }
     double getScaledGradient() { return gradobj*scalobj; }

     void   saveValue()         { saveobj = valobj;  }
     void   restoreValue()      { valobj  = saveobj; }
};

#endif

#endif

