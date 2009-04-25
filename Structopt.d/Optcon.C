#ifdef STRUCTOPT

#include <Structopt.d/Structopt_sd.h>

#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optcon.h>
#include <Structopt.d/Optopr.h>

//------------------------------------------------------------------------------
//                            Objective 
//                                          created  9/1/98 by Kurt
//------------------------------------------------------------------------------

Optcon::Optcon():valcon(0),gradcon(0),scalcon(0),savecon(0),typcon(0),opr(0) 
{
     numcon  = 0;
     numeqc  = 0;                        
     numieq  = 0;                 
}

//------------------------------------------------------------------------------

void Optcon::func() {

     int i;
     
     for (i=0;i<numcon;i++) 
	valcon[i]=opr[i]->getval(); 
}

//------------------------------------------------------------------------------

void Optcon::grad() {

     int i;
     
     for (i=0;i<numcon;i++) {
       gradcon[i]=opr[i]->getgrad(); 
     }
}

//------------------------------------------------------------------------------

void Optcon::print(FILE * optunitout) {

     int i;

     fprintf(optunitout,"\n\nConstraints:\n");
     fprintf(optunitout,"============\n\n");

     fprintf(optunitout,"Number of Constraints           : %d\n\n",numcon);
     fprintf(optunitout,"Number of Equlity Constraints   : %d\n",numeqc);
     fprintf(optunitout,"Number of Inequaltiy Constraints: %d\n",numieq);

     if(numeqc > 0) {
      
       int j=0;
      
       fprintf(optunitout,"\nEquality Constraints\n\n");

       for (i=0;i<numcon;i++) {
             
          if ( ! typcon[i] ) { 

            j++;
            fprintf(optunitout,
	      "\t%d. Equality Constraint:  Scaling Factor: %5.2f\n\n",
	      j,scalcon[i]);
	  
	    opr[i]->print(optunitout); }                  
       }
     }

     if(numieq > 0) {
      
       int j=0;
      
       fprintf(optunitout,"\nInequality Constraints\n\n");

       for (i=0;i<numcon;i++) {
       
          if ( typcon[i] ) { 

            j++;
            fprintf(optunitout,
	      "\t%d. Inequality Constraint:  Scaling Factor: %5.2f\n\n",
	      j,scalcon[i]);
	  
            opr[i]->print(optunitout); }                  
       }
     }
}

//------------------------------------------------------------------------------

void Optcon::printres(FILE * optunitout) {

     int i;

     fprintf(optunitout,"\n\tConstraints:\n");

     for (i=0;i<numcon;i++) {
             
        fprintf(optunitout,"\n\t%d. Constraint (non-scaled): %12.5e\n\n",i+1,valcon[i]);
	  
	opr[i]->printres(optunitout); 
     }                  
}

//------------------------------------------------------------------------------

void Optcon::removeConstraint(int icon)
{
  if (icon > numcon-1) {
    fprintf(stderr,"ERROR: Constraint to be removed does not exist.\n");
    exit(-1);
  }

  if (icon != numcon-1) {
    fprintf(stderr,"ERROR: Only last constraint can be removed.\n");
    exit(-1);
  }

  // delete operator

  delete opr[icon];
  opr[icon]=0;

  // adjust number of constraints

  numcon--;

  if (typcon[icon])
    numieq--;
  else
    numeqc--;
}

//------------------------------------------------------------------------------

void Optcon::saveValue()
{
  int i;
  for (i=0;i<numcon;i++) savecon[i] = valcon[i];
}

//------------------------------------------------------------------------------

void Optcon::restoreValue()
{
  int i;
  for (i=0;i<numcon;i++)  valcon[i] = savecon[i];
}

#endif
