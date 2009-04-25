#ifdef STRUCTOPT

#include <Structopt.d/Structopt_sd.h>

#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optobj.h>
#include <Structopt.d/Optopr.h>

//------------------------------------------------------------------------------
//                            Objective 
//                                          created  9/1/98 by Kurt
//------------------------------------------------------------------------------

Optobj::Optobj(){

      numobj=0;

      valobj=0.0;
      gradobj=0.0;
      scalobj=1.0;
            
      multicrit=0;      
}

//------------------------------------------------------------------------------

void Optobj::func() {
            
       valobj=multicrit->getval();
}

//------------------------------------------------------------------------------

void Optobj::grad() {
            
       gradobj=multicrit->getgrad();
}

//------------------------------------------------------------------------------

void Optobj::print(FILE * optunitout) {

     fprintf(optunitout,"\n\nObjective:\n");
     fprintf(optunitout,"==========\n\n");
     fprintf(optunitout,"Scaling Faktor: %5.2f\n\n",scalobj);
     
     multicrit->print(optunitout);
}

//------------------------------------------------------------------------------

void Optobj::printres(FILE * optunitout) {

     fprintf(optunitout,"\tObjective(s) (non-scaled): %12.5e\n\n",valobj);   
     
     multicrit->printres(optunitout);
}

#endif
