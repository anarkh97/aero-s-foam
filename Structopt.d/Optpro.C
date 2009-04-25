#ifdef STRUCTOPT

#include <stdio.h>
#include <cassert>

#include <Structopt.d/Optpro.h>

#include <Structopt.d/Optinp.h>
#include <Structopt.d/Optobj.h>
#include <Structopt.d/Optcon.h>
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Optcrit.h>
#include <Structopt.d/Optopr.h>
#include <Structopt.d/Optsol.h>
#include <Structopt.d/Relsol.h>

#include <Structopt.d/Structopt_sd.h>

#include <Driver.d/Domain.h>
#include <Structopt.d/Driver_opt.d/Domain_opt.h>


//------------------------------------------------------------------------------
//                   Optimization Problem 
//                                          created  9/1/98 by Kurt
//------------------------------------------------------------------------------
Optpro::Optpro(int _type):opc(0) {

      type     = _type;

      numcrit  = 0;
      OCCtable = 0;

      optobj = new Optobj;
      optcon = new Optcon;
      optvar = new Optvar;

      optsol = 0;
      relsol = 0;

      switch (type) 
      {
         case 0:
           optsol = new Optsol;
           break;
         case 1:
           relsol = new Relsol;
           break;
         default:
           fprintf(stderr,"ERROR: Wrong type of optimization problem\n");
      }

      numElecStcVar  = 0;
      numThermStcVar = 0;

}

//------------------------------------------------------------------------------

void Optpro::buildOCC() {

      // allocate and initialize OCCtable

      int numcon=optcon->numcon;

      OCCtable = new int*[numcon+1]; // +1 is objective
 
      int i,j;
      for (i=0;i<numcon+1;i++) {
         OCCtable[i] = new int[numcrit];
         for (j=0;j<numcrit;j++) OCCtable[i][j] = 0;
      }

      // build OCCtable
        
      optobj->multicrit->extractCriteria(OCCtable[0]);

      for (i=0;i<numcon;i++)  optcon->opr[i]->extractCriteria(OCCtable[i+1]);
}


//------------------------------------------------------------------------------

void Optpro::getActiveCriteria(int* actCon, int* actCrit) {

  int i,j;

  // initialize active Constraints

  if (actCon) {
     for (i=0;i<numcrit;i++) actCrit[i] = 0; }
  else {
     for (i=0;i<numcrit;i++) actCrit[i] = 1;
     return;
  }

  // set all criteria of objective active

  for (j=0;j<numcrit;j++) 
     if (OCCtable[0][j]) actCrit[j]=1;

  // set criteria of active constraints

  int numcon=optcon->numcon;

  for (i=0;i<numcon;i++) {
    if (actCon[i]) {
      for (j=0;j<numcrit;j++) 
         if (OCCtable[i+1][j]) actCrit[j]=1;
    }
  }
}

//------------------------------------------------------------------------------

int Optpro::checkAnalyticalSA() {
 
 int irc=0;

 switch (type)
 {
   case 0:  
     irc = optsol->getanagrad();
     break;
   case 1:
     irc = relsol->getanagrad();
     break;
 }
  
 return irc;
}

//------------------------------------------------------------------------------

void Optpro::setOutput(FILE* fp) {

  if (fp)
  {
    optunitout=fp;
    return;
  }

  //take global outputfile

  extern Domain* domain;
  Domain_opt* optdom = dynamic_cast<Domain_opt*>(domain);
  assert(optdom != 0);
 
  
  switch (type)
  {
     case 0:
       optunitout=optdom->getOptUnit();
       break;
     case 1:
       optunitout=optdom->getRelUnit();
       break;
  }
}


//------------------------------------------------------------------------------

void Optpro::print(int hflag) {

  int i;
  
  //Header    

  if (hflag)
  {
    fprintf(optunitout,"\n\n");
    fprintf(optunitout,"           Optimization-Reliability Analysis Module\n\n");
    fprintf(optunitout,"              Version: 1.01     Date: 5/12/02\n\n");
    fprintf(optunitout,"              Author: Kurt Maute\n\n");
  }

  switch (type)
  {
     case 0:
       fprintf(optunitout,"\nOptimization Criteria:\n");
       fprintf(optunitout,"======================\n\n");
       fprintf(optunitout,"Number of Optimization Criteria: %d\n\n",numcrit);
     
       for (i=0;i<numcrit;i++) opc[i]->print(optunitout);
  
       optobj->print(optunitout);
       optcon->print(optunitout);
       optvar->print(optunitout);
       break;

     case 1:
       fprintf(optunitout,"\nFailure Criteria:\n");
       fprintf(optunitout,"======================\n\n");
       fprintf(optunitout,"Number of Optimization Criteria: %d\n\n",numcrit);
     
       for (i=0;i<numcrit;i++) opc[i]->print(optunitout);
  
       optcon->print(optunitout);
       optvar->print(optunitout);
       break;

     default:
       fprintf(stderr,"ERROR: wrong problem type in optimization - reliability module\n");
       exit(-1);
  }

  fflush(optunitout);
}

//------------------------------------------------------------------------------

void Optpro::solve(Structopt* structopt) {
  
  switch (type)
  {
    case 0: //optimization
      if (optsol)
        { optsol->solve(this,structopt); }
      else 
	{
	  fprintf(stderr,"ERROR: no optimization solver specified\n");
	  exit(-1); 
	}
      break;

    case 1: //reliability
      if (relsol)
        relsol->solve(this,structopt);
      else {
        fprintf(stderr,"ERROR: no optimization solver specified\n");
        exit(-1); }
     break;

     default:
       fprintf(stderr,"ERROR: wrong problem type in optimization - reliability module\n");
       exit(-1);
  }
}

//------------------------------------------------------------------------------

void Optpro::generateRandom()
{
   optvar->generateRandom();
}

//------------------------------------------------------------------------------

int Optpro::checkFailcritOnly()
{
   int i;

   int irc=0;

   for (i=0;i<numcrit;i++) {
      if ( opc[i]->gettyp() != 15 ) irc=1;
   }

   // irc = 1  :  not only failure criteria  -> return 0
   // irc = 0  :  only failure criteria      -> return 1

   return ( irc == 1 ) ? 0 : 1;
}

//------------------------------------------------------------------------------

void Optpro::computeRelDSA(double ** gradmat,  Optcrit** lsfCrit, int numlsfCrit)
{
    // get gradient computation for current optimziation algorithm

    Optgrad * actgrad = optsol->getActOptgrad();

    // check if actgrad is analyticSA

    if ( actgrad->typ ) {
      fprintf(stderr,"WARNING: SA method for current algorithm is not analytic\n");
      fprintf(stderr,"         This is not compatible with SA of reliability");
      exit(-1);
    }       

    // add pseudo criterion

    int i;

    int* actList = new int[numlsfCrit]; 

    for (i=0;i<numlsfCrit;i++) {
       opc[numcrit] = lsfCrit[i];
       actList[i]   = numcrit;
       numcrit++;
    }

    // evaluate sensitivities of pseudo criteria

    actgrad->gradanalytic(actList,gradmat,numlsfCrit,-1);

    // remove pseudo criterion

    for (i=0;i<numlsfCrit;i++) {
       numcrit--;
       opc[numcrit] = 0;
    }

    delete [] actList;
}

//------------------------------------------------------------------------------

void Optpro::removeCriteria(int icrit)
{
   if (icrit != numcrit-1) {
     fprintf(stderr,"ERROR: only last criterion in list can be deleted\n");
     exit(-1);
   }

   delete opc[icrit];

   opc[icrit]=0;

   numcrit--;

}

//------------------------------------------------------------------------------

void Optpro::removeObjective()
{
   delete optobj->multicrit;
}

//------------------------------------------------------------------------------

void Optpro::removeConstraint(int icon)
{
   optcon->removeConstraint(icon);
}

//------------------------------------------------------------------------------

void Optpro::saveProblem(double &objSave, double* conSave, int numcon)
{
   int j;

   // save values used in solver

   objSave = optsol->obj;

   for (j=0;j<numcon;j++) { conSave[j] = optsol->con[j] ; }

   // save values stored in problem definiton

   optobj->saveValue();

   optcon->saveValue();
}

//------------------------------------------------------------------------------

void Optpro::restoreProblem(double &objSave, double* conSave, int numcon)
{ 
  int j;

  // restore values used in solver

  optsol->obj = objSave;

  for (j=0;j<numcon;j++) { optsol->con[j] = conSave[j]; }

  // restore values stored in problem definiton

  optobj->restoreValue();

  optcon->restoreValue();
}

#endif
