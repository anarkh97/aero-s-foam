#ifdef STRUCTOPT

#include <cstdlib>
#include <cstdio>

#include <Structopt.d/Optinp.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Relsol.h>
#include <Structopt.d/Optobj.h>
#include <Structopt.d/Optcon.h>
#include <Structopt.d/Optopr.h>

#include <Structopt.d/Structopt_sd.h>

//------------------------------------------------------------------------------
      
void Optpro::addFailureCriteria( int num, funcall & func, 
                                 int PMAflag, double targetRI )
{

  int inum = num-1;

  if ( optcon->opr[inum] ) {
    fprintf(stderr,
            "Error generating constraint: %d. constraint already defined\n",
            num);
    exit(-1);
  }
    
  if ( optcon->numcon < num ) optcon->numcon = num ;

  if (! PMAflag) {
    // RIA:  failure criteria is equality constraint (0), positve criteria: safe
    optcon->typcon[inum]  = 0;
    optcon->scalcon[inum] = 1.0;
  }
  else {
    // PMA:  failure criteria is new type (-1), positve criteria: safe
    optcon->typcon[inum]  = -1;
    optcon->scalcon[inum] = targetRI;
  }

  optcon->opr[inum] = addFunc(func);
  
  optcon->numeqc++;
}

//------------------------------------------------------------------------------

void Optpro::addRndvar( absvardata & var )
{
  int num=var.num-1;

  optvar->addAbsVar(num, Absvar::create(1,var)); 

  if ( var.igen ) {
  
    if ( ! optvar->getAbsVar(var.gena) ) { 
      fprintf(stderr,"Error generating rnd.var.: %d variable does not exists\n",
              var.gena+1);
      exit(-1);
    }

    if ( ! optvar->getAbsVar(var.genb) ) { 
      fprintf(stderr,"Error generating rnd.var.: %d variable does not exists\n",
              var.genb+1);
      exit(-1);
    }

    if ( optvar->getAbsVar(var.genb)->getDistribution()  !=  
         optvar->getAbsVar(var.gena)->getDistribution() )
    { 
      fprintf(stderr,"Error generating rnd.var.: variables %d and %d have different distributions\n",
              var.gena+1,var.genb+1);
      exit(-1);
    }
  
    double gnum = static_cast<double>(var.genb - var.gena);

    int    disa  = optvar->getAbsVar(var.gena)->getDistribution();   

    double meana  = optvar->getAbsVar(var.gena)->getMeanValue();
    double meanb  = optvar->getAbsVar(var.genb)->getMeanValue();
    double dmean  = (meanb - meana) / gnum;
    
    double sdeva  = optvar->getAbsVar(var.gena)->getStdDev();
    double sdevb  = optvar->getAbsVar(var.genb)->getStdDev();
    double dsdev  = (sdevb - sdeva) / gnum;

    double lowa  = optvar->getAbsVar(var.gena)->getlow();
    double lowb  = optvar->getAbsVar(var.genb)->getlow();
    double dlow  = (lowb - lowa) / gnum;

    double uppa  = optvar->getAbsVar(var.gena)->getupp();
    double uppb  = optvar->getAbsVar(var.genb)->getupp();
    double dupp  = (uppb - uppa) / gnum;
    
    int newnum;
    absvardata gvar;

    for (newnum=var.gena+1;newnum<var.genb;newnum+=var.gens) {

       num   = newnum+1;
       meana += dmean;
       sdeva += dsdev;
       lowa  += dlow;
       uppa  += dupp;
  
       if ( optvar->getAbsVar(newnum) ) {
         fprintf(stderr,
	         "Error generating abs.var.: %d. variable already defined\n",
                 num);
         exit(-1);
       }              

       gvar.num  = num;     
       gvar.mean = meana;
       gvar.sdev = sdeva;
       gvar.low  = lowa;
       gvar.upp  = uppa;

       optvar->addAbsVar(newnum, Absvar::create(1,gvar));
    }
  }
}


//------------------------------------------------------------------------------

void Relsol::addSolver( int & num, int & typ, nlpdata & param, graddata & grad )
{
    relalg[numalg]=Relalg::create(num,typ,param,grad);
  
    numalg++;                                  
}

#endif
