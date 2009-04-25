#ifdef STRUCTOPT

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <cassert>

#include <Structopt.d/Relsol.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optobj.h>
#include <Structopt.d/Optcon.h>
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Optinp.h>
#include <Structopt.d/Optopr.h>

#include <Structopt.d/Structopt_sd.h>

#include <Math.d/mathUtility.h>
#include <Timers.d/GetTime.h>

#include <Driver.d/Domain.h>
#include <Structopt.d/Driver_opt.d/Domain_opt.h>

//------------------------------------------------------------------------------

RelsolData::RelsolData() {

   lsfunction=0;
   relIndex=0;
   failureprob=0;
   pmaValue=0;


   gradlsfunction=0;
   gradrelIndex=0;
   gradfailureprob=0;
   gradpmaValue=0;

}

//------------------------------------------------------------------------------

void RelsolData::initialize(int _numlsf) 
{
   numlsf=_numlsf;

   if (!lsfunction)  lsfunction   = new double[numlsf];
   if (!failureprob) failureprob  = new double[numlsf];
}

//------------------------------------------------------------------------------

void RelsolData::initializeFORM() 
{
   if (!relIndex)    relIndex     = new double[numlsf];
   if (!pmaValue)    pmaValue     = new double[numlsf];
}

//------------------------------------------------------------------------------

void RelsolData::initializeGrad(int _numdsgvar) 
{
   numdsgvar=_numdsgvar;

   int i;

   if (!gradlsfunction) {
     gradlsfunction  = new double*[numlsf];
     for (i=0;i<numlsf;i++)
       gradlsfunction [i] = new double[numdsgvar];
   }

   if(!gradfailureprob) {
     gradfailureprob = new double*[numlsf];
     for (i=0;i<numlsf;i++)
       gradfailureprob[i] = new double[numdsgvar];
   }
}

//------------------------------------------------------------------------------

void RelsolData::initializeGradFORM() 
{
   int i;

   if (!gradrelIndex) {
     gradrelIndex = new double*[numlsf];
     for (i=0;i<numlsf;i++)
       gradrelIndex[i] = new double[numdsgvar];
   }

   if(!gradpmaValue) {
     gradpmaValue = new double*[numlsf];
     for (i=0;i<numlsf;i++)
       gradpmaValue[i] = new double[numdsgvar];
   }
}

//------------------------------------------------------------------------------

void RelsolData::cleanup() 
{
   delete [] lsfunction;

   if (gradlsfunction) {
     int i;
     for (i=0;i<numlsf;i++)
       delete [] gradlsfunction[i];
     delete [] gradlsfunction;
     gradlsfunction=0;
   }
   
   lsfunction=0;
   gradlsfunction=0;
}

//------------------------------------------------------------------------------

double RelsolData::getValue(int id, int typ)
{
   double val;

   if (id > numlsf-1) {
     fprintf(stderr,"Error: Failure criteria %d does not exist.\n",id+1);
     exit(-1);
   }

   switch(typ)
   {
     case 0:
       if (!failureprob) {
     	 fprintf(stderr,"Error: Failure probabilities data not available\n");
     	 exit(-1);
       }
       val=failureprob[id];
       break;
     case 1:
       if (!relIndex) {
     	 fprintf(stderr,"Error: Reliability index data not available\n");
     	 exit(-1);
       }
       val=relIndex[id];
       break;
     case 2:
       if (!pmaValue) {
     	 fprintf(stderr,"Error: Performance measure data not available\n");
     	 exit(-1);
       }
       val=pmaValue[id];
       break;
   }

   return val;
}

//------------------------------------------------------------------------------

double RelsolData::getGradient(int id, int ivar, int typ)
{
   double val;

   if (id > (numlsf-1) || ivar > (numdsgvar-1)) {
     fprintf(stderr,"Error: Gradient of Failure criteria %d does not exist for design variable %d.\n",id+1,ivar+1);
     exit(-1);
   }

   switch(typ)
   {
     case 0:
       if (!gradfailureprob) {
     	 fprintf(stderr,"Error: Gradient of Failure probabilities data not available\n");
     	 exit(-1);
       }
       val=gradfailureprob[id][ivar];
       break;
     case 1:
       if (!gradrelIndex) {
     	 fprintf(stderr,"Error: Gradient of Reliability index data not available\n");
     	 exit(-1);
       }
       val=gradrelIndex[id][ivar];
       break;
     case 2:
       if (!gradpmaValue) {
     	 fprintf(stderr,"Error: Gradient of Performance measure data not available\n");
     	 exit(-1);
       }
       val=gradpmaValue[id][ivar];
       break;
   }

   return val;
}

//------------------------------------------------------------------------------

Relsol::Relsol():relalg(0) {

   numalg=0;

   giter=0;
   
   numvar=0;
   numlsf=0;

   dsgopt=0;
   numdsgvar=0;

   solData = new RelsolData(); 
}

//------------------------------------------------------------------------------

int Relsol::getanagrad () {

   int i;
   int anagrad=0;

   // Watch: we get "1" form relalg if analytical SA

   for (i=0;i<numalg;i++) { if ( relalg[i]->getGradtyp() ) anagrad=1; }

   // return "1" if SA is needed
   
   return anagrad;
}  

//------------------------------------------------------------------------------

void Relsol::solve (Optpro *_relpro, Structopt *_structrel) {

   fprintf(stderr,"\n............................................\n");
   fprintf(stderr,  ".... Starting Reliability Analysis      ....\n");
   fprintf(stderr,  "............................................\n\n");

   double evalTime = getTime();

   extern Domain * domain;
   Domain_opt* optdom = dynamic_cast<Domain_opt*>(domain);
   assert(optdom != 0);
   
   relunitout=optdom->getRelUnit();

   relpro    = _relpro;
   structrel = dynamic_cast<Structopt_sd*>(_structrel);
   assert(structrel != 0);
 
   initialize();

   int i;
   for (i=0;i<numalg;i++) {
   
      liter=-1;
      relalg[i]->setOutput(relunitout);
      relalg[i]->print();
      relalg[i]->solve(this,solData);
   }
  
   // reset abstract random variables to mean value

   relpro->optvar->setMeanValue();

   cleanup();

   // Output CPU time for total function evaluation

   fprintf(stderr,"\n ... Time spent in reliability analysis: %e sec\n\n",
          (getTime()-evalTime)/1000.0);

   fprintf(stderr,"\n............................................\n");
   fprintf(stderr,  ".... Stopping Reliability Analysis      ....\n");
   fprintf(stderr,  "............................................\n\n");

}  

//------------------------------------------------------------------------------

void Relsol::initialize() {

   int i;

   numvar=relpro->optvar->nAbs();
   
   rndvar = new Absvar*[numvar];

   for (i=0;i<numvar;i++)
     rndvar[i] =  relpro->optvar->getAbsVar(i);

   numlsf=relpro->optcon->numcon;

   solData->initialize(numlsf);

   //check for driving optimization problem and analytical SA

   if ( /*structrel->structdom->*/optpro.get() != 0 ) {
 
     int agFlag = /*structrel->structdom->*/optpro->checkAnalyticalSA();

     // special treatment is only need if sensitivities in optimization problem
     // are computed analytically

     if (agFlag) {

       numdsgvar=/*structrel->structdom->*/optpro->optvar->nAbs();
       
       solData->initializeGrad(numdsgvar);
       
       dsgopt = 1;
     }
   }
}

//------------------------------------------------------------------------------

void Relsol::cleanup(){

   delete [] rndvar;

   solData->cleanup();
}

//------------------------------------------------------------------------------

void Relsol::scalfunc() {
   
   int i;
   for (i=0;i<relpro->optcon->numcon;i++)
      solData->lsfunction[i] = relpro->optcon->valcon[i] * relpro->optcon->scalcon[i];
}

//------------------------------------------------------------------------------

void Relsol::func(int iter) {

   relpro->optvar->updvar();                    //Update Structural Variables

   structrel->func();                           //Evaluation of System response

   relpro->optcon->func();                      //Evaluation of Limitstate function

   scalfunc();                                  //Scaling of Objective & Const.

   relpro->optvar->resetvar();                  //Reset Structural Variables

}

//------------------------------------------------------------------------------

void Relsol::generateRandom() 
{
   relpro->generateRandom();
}

//------------------------------------------------------------------------------

double Relsol::getFailprob(int id) { return solData->getValue(id,0); }

double Relsol::getRelindex(int id) { return solData->getValue(id,1); }

double Relsol::getPMAvalue(int id) { return solData->getValue(id,2); }


//------------------------------------------------------------------------------

double Relsol::getgradFailprob(int id, int iv) { return solData->getGradient(id,iv,0); }

double Relsol::getgradRelindex(int id, int iv) { return solData->getGradient(id,iv,1); }

double Relsol::getgradPMAvalue(int id, int iv) { return solData->getGradient(id,iv,2); }

//------------------------------------------------------------------------------

Relalg * Relalg::create( int inum, int ityp ,nlpdata& param , graddata& grad) 
{
    Relalg * newalg;
    
    int subalg;

    // consider FORM type algorithms

    if ( ityp > 99 ) {
      subalg  = ityp-100;
      ityp    = (ityp-subalg)/100;
    }

    switch (ityp) {

       case Relalg::MonteCarlo:
          newalg = new RelalgMonteCarlo(param);
	  break;
       case Relalg::Form:
          newalg = new RelalgFORM(subalg+10,param,grad);
	  break;
       default:
          fprintf(stderr," *** ERROR: Reliability analysis algorithm unkown\n");
          exit(-1);
    }
    
    newalg->num = inum;
    newalg->typ = ityp;

    return newalg;
}   

//------------------------------------------------------------------------------

RelalgMonteCarlo::RelalgMonteCarlo( nlpdata & param) {

    setDefault();

    if (param.rflag[0]) samples  = param.rval[0];
}

//------------------------------------------------------------------------------

void RelalgMonteCarlo::setDefault() {

   samples = 1.0;
}

//------------------------------------------------------------------------------

void RelalgMonteCarlo::solve(Relsol *_relsol, RelsolData* _solData) {

    relsol  = _relsol;
    solData = _solData; 

    int numlsf = relsol->numlsf;

    // initialize list of failure trails

    failure = new int[numlsf];

    int ilsf;
    for (ilsf=0;ilsf<numlsf;ilsf++) failure[ilsf]=0;

    // analytical SA in driving Dsg-Optimization not possible

    if (relsol->dsgopt)
      fprintf(stderr,"WARNING: analytical SA for MonteCarlo reliability analysis not possible\n");

    // main loop  

    char per = '%';

    double iter=0;

    int numvar = relsol->relpro->optvar->nAbs();

    srandom(unsigned(time(NULL)));

    while (iter<samples)
    {
       relsol->generateRandom();

       relsol->func(static_cast<int>(iter));

       // write and evaluate limit state functions

       fprintf(relunitout,"%c %d. Sample:\n",per,int(iter+1));

       for (ilsf=0;ilsf<numlsf;ilsf++) {

         fprintf(relunitout,"limitState(%d,%d) = %20.10e;\n",int(iter+1),ilsf+1,
                              solData->lsfunction[ilsf]);

         if (solData->lsfunction[ilsf] < 0) failure[ilsf]++;
       }

       int ivar;
       for (ivar=0;ivar<numvar;ivar++) {
         fprintf(relunitout,"randomvar(%d,%d ) = %20.10e;\n",int(iter+1),ivar+1,
		 relsol->relpro->optvar->getAbsVar(ivar)->getval());
       }

       fprintf(relunitout,"\n");

       fflush(relunitout);

       iter++;

    }

    // return failureprobability

    for (ilsf=0;ilsf<numlsf;ilsf++) 
      solData->failureprob[ilsf]=failure[ilsf]/samples;

    // print failure probabilities 

    printres();

    // delete temporary arrays
       
    delete [] failure;
}

//------------------------------------------------------------------------------

void RelalgMonteCarlo::func(int iter) {

  relsol->func(iter) ; 

}

//------------------------------------------------------------------------------

void RelalgMonteCarlo::grad() {

  fprintf(stderr,"WARNING: gradient evaluation for Monte Carlo Method is called\n"); 
  fprintf(stderr,"         nothing done\n"); 

}

//------------------------------------------------------------------------------

void RelalgMonteCarlo::print() {

   fprintf(relunitout,"\n\tReliability Analysis Strateg #%3d. : MonteCarlo\n",num);
   fprintf(relunitout,  "\t==============================================\n\n");

   fprintf(relunitout,"\tNumber of Samples ....................... : %e\n\n",samples);

}

//------------------------------------------------------------------------------

void RelalgMonteCarlo::printres() {

   fprintf(relunitout,"\n\tResult of %d. Reliability Analysis: MonteCarlo\n\n",num);
   fprintf(relunitout,  "\t====================================================\n\n");

   int numlsf = relsol->numlsf;

   int ilsf;

   for (ilsf=0;ilsf<numlsf;ilsf++) 
     fprintf(relunitout,"\tFailure Probability for Criteria %d :  %e\n",
     ilsf+1,solData->failureprob[ilsf]);

   fflush(relunitout);
}

//------------------------------------------------------------------------------

RelalgFORM::RelalgFORM(int opttyp, nlpdata& param, graddata& grad)
{
    // initialize variables

    uvals = 0;

    // create optimization problem

    optpro = new Optpro(0);

    // add solver

    int numsol=1;

    optpro->optsol->addSolver(numsol,opttyp,param,grad);
}


//------------------------------------------------------------------------------

void RelalgFORM::solve(Relsol *_relsol,RelsolData* _solData) 
{ 

  relsol  = _relsol;
  solData = _solData;

  // create arrays for solving FORM results

  solData->initializeFORM();

  if (relsol->dsgopt) solData->initializeGradFORM();

  // turn on FORM flag of random variables

  relsol->relpro->optvar->setFormFlag(1);

  // initialize FORM optimization problem

  initializeOptpro();

  // determine relibility index for each failure criteria

  int numlsf = relsol->numlsf;    

  int ilsf;
   
  for (ilsf=0;ilsf<numlsf;ilsf++)
  {
    // set RIA / PMA optimization problem

    setLSFactive(ilsf);

    // solve optimization problem

    optpro->solve(relsol->structrel);

    // store results

    saveResults(ilsf);

    // perform analytical SA for driving Dsg-Optimization 

    if (relsol->dsgopt) analyticDSA(ilsf);

    // reset RIA / PMA optimization problem

    resetLSFactive(ilsf);

  }

  // print reliability results

  printres();

  // clean up FORM optimization problem

  cleanOptpro();

  // turn on FORM flag of random variables

  relsol->relpro->optvar->setFormFlag(0);

}

//------------------------------------------------------------------------------

void RelalgFORM::initializeOptpro()
{

   // set output file in optimization solver

   optpro->optsol->optprotfile=relsol->relprotfile;
   optpro->optsol->fsize=relsol->fsize;

   // set optimization variables

   optpro->optvar=relsol->relpro->optvar; 

   // set optimization variables to mean value

   optpro->optvar->setMeanValue();

   // add reliability index (beta) as LAST optimization criterion

   int numcrit=relsol->relpro->numcrit;

   numcrit++;

   criterion RIcriterion;
   RIcriterion.initialize();

   RIcriterion.num=numcrit;
   RIcriterion.data.typ=14;

   relsol->relpro->addCriteria(RIcriterion);

   // copy LSF as optimization criteria

   numcrit=relsol->relpro->numcrit;

   optpro->numcrit=numcrit;

   int icrit;

   for (icrit=0;icrit<numcrit;icrit++)
     optpro->opc[icrit]=relsol->relpro->opc[icrit]; 

   // allocate memory for storing results

   if (!uvals) {

     int numlsf = relsol->numlsf;    

     uvals = new double*[numlsf];

     int numvar = relsol->relpro->optvar->nAbs();

     int ilsf,ivar;
     for (ilsf=0;ilsf<numlsf;ilsf++) {
       uvals[ilsf] = new double[numvar]; 
       for (ivar=0;ivar<numvar;ivar++) 
         uvals[ilsf][ivar]=relsol->relpro->optvar->getAbsVar(ivar)->getMeanValue();
     }
   }
}

//------------------------------------------------------------------------------

void RelalgFORM::setLSFactive(int ilsf)
{
   // check for FORM - typ: RIA / PMA 

  int FORMtyp = static_cast<int>(relsol->relpro->optcon->typcon[ilsf]);

   // construct optimization problem

   switch ( FORMtyp ) {

     case 0:                    // RIA
        setRIAproblem(ilsf);
        break;
     case -1:                   // PMA
        setPMAproblem(ilsf);
        break;
     default:
        fprintf(stderr,"ERROR: FORM type not implemented\n"); 
        exit(-1);
   }
}

//------------------------------------------------------------------------------

void RelalgFORM::resetLSFactive(int ilsf)
{
   // check for FORM - typ: RIA / PMA 

  int FORMtyp = static_cast<int>(relsol->relpro->optcon->typcon[ilsf]);

   // construct optimization problem

   switch ( FORMtyp ) {

     case 0:                    // RIA
        resetRIAproblem(ilsf);
        break;
     case -1:                   // PMA
        resetPMAproblem(ilsf);
        break;
     default:
        fprintf(stderr,"ERROR: FORM type not implemented\n"); 
        exit(-1);
   }
}

//------------------------------------------------------------------------------

void RelalgFORM::setRIAproblem(int ilsf)
{
   // reset optimization variables to mean value

   optpro->optvar->setMeanValue(uvals[ilsf]);

   // number of optimization criteria 
  
   int numcrit=optpro->numcrit;

   // build optjective:  1/initial_value * SUM { CRIT[] }

   funcall   objfunc;
   funcdata* objdata=buildFuncdata();

   objfunc.typ=0;
   objfunc.fdata=objdata;

   objdata->numopr=1;
   objdata->numgen=0;
   objdata->oprtyp[0]=0; 
   objdata->oprnum[0]=numcrit-1;
   objdata->a	  [0]=1.0; 
   objdata->p	  [0]=1.0;
   objdata->b	  [0]=0.0; 

   // for augmented formulation
   
   fprintf(stderr," ... WATCH: augmented formulation for RIA used\n");
   
   objdata->numopr=2;
   objdata->oprtyp [1]=1; 
   objdata->oprnum [1]=0;
   objdata->a	   [1]=1.0; 
   objdata->p	   [1]=2.0;
   objdata->b	   [1]=0.0; 
   objdata->subfunc[0]=relsol->relpro->optcon->opr[ilsf]->extractFuncall();
   
   double l2norm2=optpro->optvar->getFormL2norm2();
   
   //if (l2norm2 == 0 )
   //  l2norm2 = 1.0;
   //else
   //  l2norm2=abs(1.0/l2norm2);
 
   fprintf(stderr," ... WATCH: scaling of objective in MPP (RIA) search skipped\n");

   l2norm2 = 1.0;

   optpro->addObjective(l2norm2,objfunc);

   // optimization problem has only one equality constraint

   optpro->optcon->numcon=1;
   optpro->optcon->numeqc=1;

   // set optimization constraint

   optpro->optcon->typcon [0] = 0;
   optpro->optcon->scalcon[0] = 1.0;
   optpro->optcon->opr    [0] = relsol->relpro->optcon->opr[ilsf];

   // print optimization problem ignoring constraints
   
   optpro->print(0);
}

//------------------------------------------------------------------------------

void RelalgFORM::resetRIAproblem(int ilsf)
{
   // remove objective in optpro

   optpro->removeObjective();
}

//------------------------------------------------------------------------------

void RelalgFORM::setPMAproblem(int ilsf)
{
   // reset optimization variables to mean value

   optpro->optvar->setMeanValue(uvals[ilsf]);

   // number of optimization criteria 
  
   int numcrit=optpro->numcrit;

   // build optjective:    

   optpro->addObjective(1.0,relsol->relpro->optcon->opr[ilsf]);

   // reset number of constraints

   optpro->optcon->numcon=0;
   optpro->optcon->numeqc=0;

   // build constraint on reliability index:  
   // SUM { CRIT[last criteria] / target reliability - 1.0 }

   double targetRI = relsol->relpro->optcon->scalcon[ilsf];

   funcall   constfunc;
   funcdata* constdata=buildFuncdata();

   constfunc.typ=0;
   constfunc.fdata=constdata;

   constdata->numopr=1;
   constdata->numgen=0;
   constdata->oprtyp[0]=0; 
   constdata->oprnum[0]=numcrit-1;
   constdata->a	  [0]=1.0/targetRI/targetRI; 
   constdata->p	  [0]=1.0;
   constdata->b	  [0]=-1.0; 
  
   optpro->addConstraint(1,0,1.0,constfunc);

   // print optimization problem ignoring constraints
   
   optpro->print(0);
}

//------------------------------------------------------------------------------

void RelalgFORM::resetPMAproblem(int ilsf)
{
   // remove constraint 

   optpro->removeConstraint(0);
}

//------------------------------------------------------------------------------

void RelalgFORM::saveResults(int ilsf)
{
   // check for FORM - typ: RIA / PMA 

  int FORMtyp = static_cast<int>(relsol->relpro->optcon->typcon[ilsf]);

   // construct optimization problem

   switch ( FORMtyp ) {

     case 0:                    // RIA
        saveRIAresults(ilsf);
        break;
     case -1:                   // PMA
        savePMAresults(ilsf);
        break;
     default:
        fprintf(stderr,"ERROR: FORM type not implemented\n"); 
        exit(-1);
   }

   // save optimization variables 

   int numvar = relsol->relpro->optvar->nAbs();

   int ivar;

   int zeroFlag = 0;

   if (zeroFlag) {
     fprintf(stderr," ... WATCH: FORM variables are reset to zero\n");
     for (ivar=0;ivar<numvar;ivar++) 
       uvals[ilsf][ivar]=0.0;
   }
   else {
     fprintf(stderr," ... WATCH: FORM variables are set to previous results\n");
     for (ivar=0;ivar<numvar;ivar++) 
       uvals[ilsf][ivar]=optpro->optvar->getAbsVar(ivar)->getval();
   }
}

//------------------------------------------------------------------------------

void RelalgFORM::saveRIAresults(int ilsf)
{
    // evaluate reliability index

    double b2 = optpro->optobj->valobj;    
    double sb = optpro->optsol->relSwitch;

    double sbeta = sb*sqrt(b2);

    solData->relIndex[ilsf] = sbeta;

    // evaluate failure probabilities

    solData->failureprob[ilsf]=normCdf(-sbeta); 

    // evaluate performance value

    solData->pmaValue[ilsf]=optpro->optcon->valcon[0]; 
}

//------------------------------------------------------------------------------

void RelalgFORM::savePMAresults(int ilsf)
{
    // evaluate performance value

    solData->pmaValue[ilsf]=optpro->optobj->valobj; 

    // evaluate reliability index

    double targetRI = relsol->relpro->optcon->scalcon[ilsf];

    double con  = optpro->optcon->valcon[0];

    double b2   = (1.0+con)*targetRI*targetRI; 

    double beta = sqrt(b2);

    solData->relIndex[ilsf] = beta;

    // evaluate failure probabilities

    solData->failureprob[ilsf]=normCdf(-beta); 
}

//------------------------------------------------------------------------------

void RelalgFORM::analyticDSA(int ilsf)
{ 
   // set most probable point

   relsol->relpro->optvar->updvar();     

   // get norm of limit state function gradients

   double gnorm = optpro->optsol->getGnorm();

   // extract criteria of current failure criterion
   // which is assumed to be
   // for RIA : the first and only constraint
   // for PMA : the objective

   int* occtab; 

   int FORMtyp = static_cast<int>(relsol->relpro->optcon->typcon[ilsf]);

   switch ( FORMtyp ) {
     case 0:                    
        occtab = optpro->OCCtable[1];   // case: RIA
        break;
     case -1:                   
        occtab = optpro->OCCtable[0];   // case: PMA
        break;
   }

   int numcrit = optpro->numcrit;

   int numActCrit=0;

   int i;
   for (i=0;i<numcrit;i++) 
      if (occtab[i] > 0) numActCrit++;

   Optcrit** critList = new Optcrit*[numActCrit];
   
   numActCrit=0;
   for (i=0;i<numcrit;i++) {
      if (occtab[i] > 0) {
         critList[numActCrit] = optpro->opc[i];
         numActCrit++;
      }
   }

   // temporary array for sensitivities

   int numdsgvar = relsol->numdsgvar;

   double *  pgradmat = new double [numActCrit*numdsgvar];
   double **  gradmat = new double*[numActCrit];

   for (i=0;i<numActCrit;i++) gradmat[i]=pgradmat+i*numdsgvar;  

   // call optimization problem to determine SA of individual criteria

   /*relsol->structrel->structdom->*/optpro->computeRelDSA(gradmat,critList,numActCrit);

   // compute sensitivities of reliability index (RIA) and performance measure (PMA)

   // Watch: gradients of reliability index   are only correct in the context of RIA approach 
   // Watch: gradients of performance measure are only correct in the context of PMA approach 

   int j;

   for (i=0;i<relsol->numdsgvar;i++) {
      
      for (j=0;j<numActCrit;j++) critList[j]->grad = gradmat[j][i];
   
      double dgds;

      switch ( FORMtyp ) {
        case 0:                    
          dgds = optpro->optcon->opr[0]->getgrad();  // case: RIA
        break;
        case -1:                   
          dgds = optpro->optobj->multicrit->getgrad();  // case: PMA
        break;
      }
  
      fprintf(stderr,"\n\n ... for dsg.var. %d:  dgds = %e\n",i+1,dgds);

      solData->gradrelIndex[ilsf][i] = dgds/gnorm;

      fprintf(stderr,"\n ... for dsg.var. %d:  gradBeta = %e\n",i+1,solData->gradrelIndex[ilsf][i]);

      solData->gradpmaValue[ilsf][i] = dgds;
   }

   // compute sensitivities of failure probability 

   // Watch: gradients of failure probability are only correct in the context RIA approach 

   double sqotpi = 1.0/sqrt(8.0*atan(1.0));

   double sbeta  = solData->relIndex[ilsf];

   fprintf(stderr,"\n ... Beta = %e\n",sbeta);
 
   double eb2o2  = exp(-sbeta*sbeta/2.0);

   for (i=0;i<relsol->numdsgvar;i++) {

      double gp = - sqotpi * eb2o2 * solData->gradrelIndex[ilsf][i]; 

      fprintf(stderr,"\n ... for dsg.var. %d:  gradProb = %e\n",i+1,gp);

      solData->gradfailureprob[ilsf][i]= gp;
   }

   // reset random structural variables

    relsol->relpro->optvar->resetvar();                  

   // delete temporary arrays

   delete [] critList;
   delete [] pgradmat;
   delete [] gradmat;
}

//------------------------------------------------------------------------------

void RelalgFORM::cleanOptpro() 
{
   int numcrit=relsol->relpro->numcrit;

   // remove last criteria in relpro

   relsol->relpro->removeCriteria(numcrit-1);
}

//------------------------------------------------------------------------------

void RelalgFORM::func(int iter) 
{
  fprintf(stderr,"in FORM::func - do nothing\n");

  relsol->func(iter) ; 
}

//------------------------------------------------------------------------------

void RelalgFORM::grad() 
{
  fprintf(stderr,"in FORM::grad - do nothing\n");
}

//------------------------------------------------------------------------------

void RelalgFORM::setOutput(FILE* op)
{ 
   relunitout=op; 
   optpro->setOutput(op); 
}

//------------------------------------------------------------------------------

void RelalgFORM::print() {   


   fprintf(relunitout,"\nReliability Analysis Strateg #%3d. : FORM\n",num);
   fprintf(relunitout,  "==============================================\n\n");

}

//------------------------------------------------------------------------------

void RelalgFORM::printres() 
{

   fprintf(relunitout,"\nResult of %d. Reliability Analysis: FORM\n\n",num);
   fprintf(relunitout,  "===============================================\n\n");

   int numlsf = relsol->numlsf;

   int ilsf;

   for (ilsf=0;ilsf<numlsf;ilsf++) {

     int FORMtyp = static_cast<int>(relsol->relpro->optcon->typcon[ilsf]);

     switch ( FORMtyp ) {
       case 0:                   
         fprintf(relunitout,"\tCriteria %d: Reliability Index Approach\n\n",ilsf+1);
         break;
       case -1:
         fprintf(relunitout,"\tCriteria %d: Performance Measure Approach\n\n",ilsf+1);
         break;
       default:
         fprintf(stderr,"ERROR: FORM type not implemented\n"); 
         exit(-1);
     }

     fprintf(relunitout,"\tReliability Index   :  %e\n",
     solData->relIndex[ilsf]);

     fprintf(relunitout,"\tFailure Probability :  %e\n",
     solData->failureprob[ilsf]);

     fprintf(relunitout,"\tPreformance Index   :  %e\n",
     solData->pmaValue[ilsf]);

     fprintf(relunitout,"\n\n\n");
   }

   fflush(relunitout);

}

//------------------------------------------------------------------------------

int RelalgFORM::getGradtyp () {

   return optpro->checkAnalyticalSA();
}  

//------------------------------------------------------------------------------

double RelalgFORM::normCdf(double x)
{
 double sum,val,d,did;
 long i,n;
 int sign;

 double pi=4.0*atan(1.0);

 sum=0.0;
 val=0.0;
 sign=1;

 if (x<0) { 
   x=-x;
   sign=-1;
 }

 if(x<.000001) return(.5);

 if(sign==1&&x>6) return(1.0);

 if(sign==-1&&x>6) return(0.0);

 n=static_cast<long>(floor(100000*x));
 d=x/n;

 for(i=1;i<n;i++) {
   did = i*d;
   sum+=2.0*exp((-.5)*did*did);
 }

 sum+=1+exp((-.5)*x*x);

 val=.5+(1/sqrt(2*pi))*(d/2)*sum;

 if(sign==1) 
   return(val);
 else 
   return(1-val);
}

#endif
