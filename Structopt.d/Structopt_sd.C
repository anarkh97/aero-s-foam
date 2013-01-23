#ifdef STRUCTOPT

#include <cassert>
#include <unistd.h>

#include <Utils.d/MyComplex.h>

#include <Structopt.d/Optinp.h>
#include <Structopt.d/Structopt_sd.h>
#include <Structopt.d/Relsol.h>

#include <Math.d/mathUtility.h>
#include <Element.d/Element.h>
#include <Timers.d/GetTime.h>

//------------------------------------------------------------------------------

Structopt_sd::Structopt_sd(int tp, Domain_opt* d, Optpro* o):elecattrList(0),thermattrList(0)
{

  type=tp;

  structdom = d;
  optpro    = o;

  analysisData = 0;  
    
  staticpros    = 0;
  staticsolvers = 0;
  c_staticpros    = 0;
  c_staticsolvers = 0;
  numStaticProb = 0;

  /*
    nonlinpros    = 0;
    nonlinsolvers = 0;  
    numNlnstcProb = 0;
  
    eigenpros     = 0;
    eigensolvers  = 0;
    numEigenProb  = 0;

    dynpros       = 0;
    dynsolvers    = 0;
    numDynamProb  = 0;
    
    NLdynpros     = 0;
    NLdynsolvers  = 0;
    numNLdynProb  = 0;

    numFVibrProb  = 0;
    fvibrpros     = 0;
    fvibrsolvers  = 0;
  */

  numvar  = 0 ;
  numFunc = 0 ; 

  anagrdType = 0;
  optInfset  = 0;

  reliabilityFlag = 0;
  reliabilityProb = 0;
  reliabilityStrc = 0;

  failcritOnly = 0;

  designoptFlag = 0;
  designoptStrc = 0;

  maxstep       = 0;
  sendInitDisp  = 1;

  numCurElecAttrVar = 0;

  elecattr     = 0; 
  gradelecattr = 0;

  numElecAttrVar     = 0;
  numGradElecAttrVar = 0;

  numCurThermAttrVar  = 0;
  numThermAttrVar     = 0;
  numGradThermAttrVar = 0;  
   
  thermattr     = 0; 
  gradthermattr = 0;
  
}
 
//------------------------------------------------------------------------------

void Structopt_sd::build(Domain_opt *domain, Optpro *_optpro) 
{    
  //initialize optimization procedure  
  structdom =  domain;
  optpro    = _optpro;  
  //create structural problem types
  initStructProbSol();  
  //initialize reliability analysis
  if ( structdom->getReliabilityFlag() ) initReliabilityAnalysis();

  //set criteria and structural variables
  optpro->optvar->setvalptr(this);
  optpro->optvar->updvar();
  //initialize arrays for analytical sensitivity analysis
  //if ( optpro->checkAnalyticalSA() ) 
    { initAnalyticalSA(); }
  
#ifdef AEROELASTIC
  //set aeroelastic flag and get fluid parameters from fluid
  aeroact = 1; 
  structdom->rcvExtDomainpar();
     
  //handle non-structural based criteria
  structdom->informExtCrit(optpro);  
#endif

  //open output-files  
  geoSource->openOutputFiles();               
  return;
}

//------------------------------------------------------------------------------

void Structopt_sd::initStructProbSol(Structopt_sd* stcopt) // default argument is 0
{
  //implemented problem types: single domain static  analysis
  //                           single domain modal   analysis
  //                           single domain dynamic analysis
  
  //watch ! last modal analysis since this can be subproblem of
  //        static or dynamic problems 
  
  numAnalysis   = structdom->getNumAnalysis();  
  numStaticProb = structdom->getNumStaticAna();

  /*
  numNlnstcProb = structdom->getNumNlnstcAna();
  numEigenProb  = structdom->getNumEigenAna();
  numDynamProb  = structdom->getNumDynamAna();
  numNLdynProb  = structdom->getNumNLdynAna();
  numFVibrProb  = structdom->getNumFVibrAna();
  */
  
  // initialize storage arrays for problem types
  if (numStaticProb) 
    {
      if(structdom->isComplex())
	{
	  c_staticpros    = new SingleDomainStatic_opt<DComplex,ComplexVector,ComplexSolver>*[numStaticProb];
	  c_staticsolvers = new StaticSolver_opt<DComplex, 
	    ComplexSolver, ComplexVector,
	    SingleDomainPostProcessor_opt<DComplex, ComplexVector, ComplexSolver>,
	    SingleDomainStatic_opt<DComplex,ComplexVector,ComplexSolver>,
	    ComplexVector >
	    *[numStaticProb];
	}
      else
	{
	  staticpros    = new SingleDomainStatic_opt<double,Vector,Solver>*[numStaticProb];
	  staticsolvers = new StaticSolver_opt<double, 
	    Solver, Vector,
	    SingleDomainPostProcessor_opt< double, Vector, Solver>,
	    SingleDomainStatic_opt<double,Vector,Solver>,
	    ComplexVector >
	    *[numStaticProb];
	}
    }
  
  /*
  if (numNlnstcProb) {
    nonlinpros    = new NonLinStatic*[numNlnstcProb];
    nonlinsolvers = new NLStaticSolver<Solver,Vector,
      SingleDomainPostProcessor<double,Vector>,
      NonLinStatic, GeomState  >*[numNlnstcProb];
  }
  
  if (numEigenProb) {
    eigenpros    = new SingleDomainEigen*[numEigenProb];
    eigensolvers = new EigenSolver<DynamMat,Vector,VectorSet,
      SDEigenPostProcessor,
      SingleDomainEigen>*[numEigenProb];
  }
  
  if (numDynamProb) {
    dynpros      = new SingleDomainDynamic*[numDynamProb];
    dynsolvers   = new DynamicSolver<DynamMat,Vector,SDDynamPostProcessor,
      SingleDomainDynamic>*[numDynamProb];
  }
  
  if (numNLdynProb) {
    NLdynpros    = new NonLinDynamic*[numNLdynProb];
    NLdynsolvers = new NLDynamSolver <Solver, Vector,
      SDDynamPostProcessor, NonLinDynamic,
      GeomState>*[numNLdynProb];
  }
  
  if (numFVibrProb)
    {
      fvibrpros     = new SingleDomainFVibr*[numFVibrProb];
      fvibrsolvers  = new FVibrSolver*[numFVibrProb];
    }
  */
  
  // allocate and initialize analysis maps
  
  analysisMap = new int[numAnalysis];
  
  int countStaticProb = 0;
  int countNlnstcProb = 0;
  int countEigenProb  = 0;
  int countDynamProb  = 0;
  int countNLdynProb  = 0;
  int countFVibrProb  = 0;
  
  int ina;     
  for (ina=0;ina<numAnalysis;ina++) {
    
    structdom->activateAnalysis(ina,0);  
    
    switch (structdom->probType()) {
      
    case SolverInfo::Static:
      analysisMap[ina] = countStaticProb;
      countStaticProb++;
      break;
      /*
    case SolverInfo::NonLinStatic:
      analysisMap[ina] = countNlnstcProb;
      countNlnstcProb++;
      break;
    case SolverInfo::Modal:
      analysisMap[ina] = countEigenProb;
      countEigenProb++;
      break;
    case SolverInfo::Dynamic:
      analysisMap[ina] = countDynamProb;
      countDynamProb++;
      break;
    case SolverInfo::NonLinDynam:
      analysisMap[ina] = countNLdynProb;
      countNLdynProb++;
      break;
    case SolverInfo::FVibr:
      analysisMap[ina] = countFVibrProb;
      countFVibrProb++;
      break;
      */
    default:
      fprintf(stderr, "Problem Type %d is not implemented in Optimization Module\n", structdom->probType());
      exit(-1);
    }
  }
  
  // initialize individual analysis problems
  
  for (ina=0;ina<numAnalysis;ina++) {
    
    structdom->activateAnalysis(ina,0);  
    
    int apb = analysisMap[ina];
    
    switch (structdom->probType()) {
      
    case SolverInfo::Static:
      if (stcopt) {
	staticpros[apb]    =  stcopt->staticpros[apb];
	staticsolvers[apb] =  stcopt->staticsolvers[apb];
      }
      else {
	if(structdom->isComplex())
	  {
	    c_staticpros[apb] = new SingleDomainStatic_opt<DComplex,ComplexVector,
	      ComplexSolver>(structdom);
	    
	    c_staticsolvers[apb] = new StaticSolver_opt<DComplex, ComplexSolver, ComplexVector,
	      SingleDomainPostProcessor_opt<DComplex, ComplexVector, ComplexSolver>,
	      SingleDomainStatic_opt<DComplex,ComplexVector,ComplexSolver>,
	      ComplexVector >
	      (c_staticpros[apb]);
	
	    c_staticsolvers[apb]->optPreSolve();
	  }
	else
	  {
	    staticpros[apb] = new SingleDomainStatic_opt<double,Vector,
	      Solver>(structdom);
	    
	    staticsolvers[apb] = new StaticSolver_opt<double, Solver, Vector,
	      SingleDomainPostProcessor_opt<double, Vector, Solver>,
	      SingleDomainStatic_opt<double,Vector,Solver>,
	      ComplexVector >
	      (staticpros[apb]);
	
	    staticsolvers[apb]->optPreSolve();
	  }
      }
      
      break;
      /*
    case SolverInfo::NonLinStatic:	  
      if (stcopt) {
	nonlinpros[apb]    =  stcopt->nonlinpros[apb];
	nonlinsolvers[apb] =  stcopt->nonlinsolvers[apb];
      }
      else {
	nonlinpros[apb]    = new NonLinStatic(structdom);	 
	nonlinsolvers[apb] = new NLStaticSolver <Solver,Vector,
	  SingleDomainPostProcessor<double,Vector>,
	  NonLinStatic, GeomState> (nonlinpros[apb]);
	
	nonlinsolvers[apb]->NLoptPreSolve();
      }
      
      break;	      
    case SolverInfo::Modal:
      if (stcopt) {
	eigenpros[apb]    = stcopt->eigenpros[apb];
	eigensolvers[apb] = stcopt->eigensolvers[apb];
      }
      else {
	eigenpros[apb]    = new SingleDomainEigen(structdom);
	
	eigensolvers[apb] =  EigenSolver<DynamMat, Vector, VectorSet,
	  SDEigenPostProcessor,SingleDomainEigen>
	  ::buildEigenSolver(eigenpros[apb]);
	
	eigensolvers[apb]->optPreSolve();
      }	
      break;
    case SolverInfo::Dynamic:
      if (stcopt) {
	dynpros[apb]    = stcopt->dynpros[apb];
	dynsolvers[apb] = stcopt->dynsolvers[apb];
      }
      else {
	dynpros[apb]    = new SingleDomainDynamic(structdom);

	dynsolvers[apb] = new DynamicSolver <DynamMat,Vector,
	  SDDynamPostProcessor,
	  SingleDomainDynamic>
	  (dynpros[apb]);
 
	dynsolvers[apb]->optPreSolve(this);
      }
      break;
    case SolverInfo::NonLinDynam:
      if (stcopt) {
	NLdynpros[apb]    = stcopt->NLdynpros[apb];
	NLdynsolvers[apb] = stcopt->NLdynsolvers[apb];
      }
      else {
	NLdynpros[apb] = new NonLinDynamic(structdom);

	NLdynsolvers[apb] = new NLDynamSolver <Solver, Vector,
	  SDDynamPostProcessor, NonLinDynamic,
	  GeomState>(NLdynpros[apb]);
 
	NLdynsolvers[apb]->optPreSolve(this);
      }  
      break;
    case SolverInfo::FVibr:
      if (stcopt) {
	fvibrpros[apb]    =  stcopt->fvibrpros[apb];
	fvibrsolvers[apb] =  stcopt->fvibrsolvers[apb];
      }
      else {
	fvibrpros[apb] = new SingleDomainFVibr(structdom);
	fvibrsolvers[apb] = new FVibrSolver(fvibrpros[apb]);
	fvibrsolvers[apb]->optPreSolve();
      }
      break;
      */
    default:
      fprintf(stderr, "Problem Type %d is not implemented in Optimization Module\n", structdom->probType());    	  
      exit(-1);
    }
  }
}

//------------------------------------------------------------------------------

void Structopt_sd::initReliabilityAnalysis()
{
  // check if calling module is already reliability analysis 

  if (type) return;

  structdom->reloptInitialize(this);

  reliabilityFlag = 1;
  reliabilityProb = /*structdom->*/relpro.get();
  reliabilityStrc = /*structdom->*/dynamic_cast<Structopt_sd*>(structrel.get());

  failcritOnly = optpro->checkFailcritOnly();
}

//------------------------------------------------------------------------------
void Structopt_sd::initAnalyticalSA()
{

  int i,j,ina;

  // build gradient arrays
  structdom->buildOptGrad();
  /*
  // assign gradient arrays to nonlinear problems
  for (ina=0;ina<numNlnstcProb;ina++) nonlinpros[ina]->buildOptGradNL();
  for (ina=0;ina<numNLdynProb;ina++)  NLdynpros[ina]->buildOptGradNL();
  */

  //initialize pointers for derivatives
  optpro->optvar->setgradptr(this);
  //initialize all abstract and physical variables
  optpro->optvar->initgrad();
  
  //build optimization variable - elmement influence table
  //if not already read parsing the optimization input file

  int numabs = optpro->optvar->nAbs();

  if ( ! optInfset ) {
  
    optInf = new OptActInfo*[numabs];

    for(i=0;i<numabs;i++) {
      optpro->optvar->updgrad(i,1.0);
      optInf[i]=structdom->buildOptInf();
    }
  } 

  // write optimization variable - elmement influence table
 
  FILE * foptinf = fopen("optinf.data","w");
 
  for (i=0;i<numabs;i++) {
 
    int  size = optInf[i]->size();
    char* typ = optInf[i]->typname();
   
    fprintf(foptinf,"# Abstract variable: %d\n",i+1); 
   
    for (j=0;j<size;j++) {
      fprintf(foptinf,"%d   %d   %s\n",i+1,(optInf[i]->getEntry(j))+1,typ); 
    }

  }  
 
  fclose(foptinf); 

}

//------------------------------------------------------------------------------
void Structopt_sd::buildInopt(Domain_opt* domain, Optpro* relpro, Structopt_sd* stcopt) 
{

  //initialize optimization procedure

  structdom =  domain;
  optpro    =  relpro;
 
  //set flag and save outer optimization interface

  designoptFlag = 1;
  designoptStrc = stcopt;

  //create structural problem types

  initStructProbSol(stcopt);
 
  //set criteria and structural variables
  numvar  = optpro->optvar->nStc();
  optpro->optvar->setvalptr(this);

  //initialize arrays for analytical sensitivity analysis
  if ( optpro->checkAnalyticalSA() ) initAnalyticalSA() ;  

  //set aeroelastic flag, fluid parameters obtained only by master problem


#ifdef AEROELASTIC

  aeroact = 1; 

#endif

}

//------------------------------------------------------------------------------

void Structopt_sd::cleanup() {
  // so far only dummy routine

  int ifail=1;

  int ina;
  for (ina=0;ina<numAnalysis;ina++) {
 
    structdom->activateAnalysis(ina);  

    switch (structdom->probType()) {

    case SolverInfo::Static:
      ifail=0;
      break;
      /*
    case SolverInfo::NonLinStatic:
      ifail=0;
      break;
    case SolverInfo::Modal:
      ifail=0;
      break;
    case SolverInfo::Dynamic:      
      ifail=0;
      break;
    case SolverInfo::NonLinDynam:
      ifail=0;
      break;
    case SolverInfo::FVibr:
      ifail=0; 
      break;
      */
    default:
      assert(0);
    }
  }      
  
  if (ifail) {
    fprintf(stderr,"Problem Type not implemented in Optimization Module");
    exit(-1);
  }
}

//------------------------------------------------------------------------------

void Structopt_sd::initOptInf(int numvar)
{
  optInf = new OptActInfo*[numvar];
 
  for(int ivar=0; ivar<numvar; ++ivar)
    optInf[ivar] = new OptActInfo; 
}

//------------------------------------------------------------------------------

void Structopt_sd::setoptInf(int ivar)
{
  structdom->setoptInf(optInf[ivar]);
}


//------------------------------------------------------------------------------

double * Structopt_sd::getptrnodattr( int& loc1, int& loc2 ) {

  double *p=0;
	  
  return p;
}

//------------------------------------------------------------------------------

double * Structopt_sd::getptrnodcord( int& loc1, int& loc2 ) {

  double *p=0;

  CoordSet &nodes = structdom->getNodes();
	   
  switch (loc2) {
	   	   
  case 0: 
    p=&(nodes[loc1]->x); 
    break;
  case 1: 
    p=&(nodes[loc1]->y); 
    break;
  case 2: 
    p=&(nodes[loc1]->z); 
    break;
  default: 
    p=0; 
  }
      	  
  return p;
}


//------------------------------------------------------------------------------

double * Structopt_sd::getptrnodalforce( int& loc1, int& loc2 ) {

  // variable force only works for one analysis so far
  // quick solution: analysis Id and load case Id are set here
      
  int anaId = 0;
  int anaLC = 0;

  if (numAnalysis > 1){
    fprintf(stderr," *** ERROR: variable forces only possible\n");
    fprintf(stderr,"            for 1 analysis\n");
    exit(-1);
  }

  structdom->activateAnalysis(anaId);  

  int numLCs = structdom->getNumLC(); 

  if (numLCs > 1){
    fprintf(stderr," *** ERROR: variable forces only possible\n");
    fprintf(stderr,"            for 1 loadcase\n");
    exit(-1);
  }

  // get pointers to variable force
      
  int numNBC = structdom->nNeumann(anaId,anaLC);  
  BCond* nbc = structdom->getNBC(anaId,anaLC);  

  double *p=0;
 
  int ibc;
  for(ibc=0;ibc<numNBC;ibc++)
    if (nbc[ibc].nnum == loc1 && nbc[ibc].dofnum == loc2)
      {
	p = &(nbc[ibc].val);
	break;
      }
	
  if ( !p )	
    {	    	  
      fprintf(stderr,"ERROR:  Optimization variable: Force in Node %d\n",loc1+1);
      fprintf(stderr,"        for DOF %d could not be set. Stop.\n\n",loc2+1);
      exit(-1);
    }
      
  return p;	     
}


//------------------------------------------------------------------------------

double * Structopt_sd::getptrgradnodattr( int& loc1, int& loc2 ) {

  double *p=0;
	  
  return p;
}

//------------------------------------------------------------------------------

double * Structopt_sd::getptrgradnodcord( int& loc1, int& loc2 ) {

  double *p=0;

  CoordSet &gradnodes = structdom->getGradNodes();
	   
  switch (loc2) {
	   	   
  case 0: 
    p=&(gradnodes[loc1]->x); 
    break;
  case 1: 
    p=&(gradnodes[loc1]->y); 
    break;
  case 2: 
    p=&(gradnodes[loc1]->z); 
    break;
  default: 
    p=0; 
  }
      	  
  return p;
}

//------------------------------------------------------------------------------

double * Structopt_sd::getptrgradnodalforce( int& loc1, int& loc2 ) {

  int numNBC     = structdom->nNeumann();  
  BCond* gradnbc = structdom->getGradNBC();  

  double *p=0;
 
  int ibc;
  for(ibc=0;ibc<numNBC;ibc++)
    if (gradnbc[ibc].nnum == loc1 && gradnbc[ibc].dofnum == loc2 )
      {
	p = &(gradnbc[ibc].val);
	break;
      }
	
  if ( !p )	
    {	    	  
      fprintf(stderr,"ERROR:  Optimization variable: Force in Node %d\n",loc1+1);
      fprintf(stderr,"        for DOF %d could not be set. Stop.\n\n",loc2+1);
      exit(-1);
    }
      
  return p;	     
}


//------------------------------------------------------------------------------

double * Structopt_sd::getptrfluidattr( int& loc1) {

  double *p=0;

#ifdef AEROELASTIC	  
  double *fluidattr = structdom->getFluidAttr();
	
  if (!fluidattr) return p;   	      
   
  switch (loc1) {
  case 0: 
    p=&(fluidattr[0]); 
    break;
  case 1: 
    p=&(fluidattr[1]); 
    break;
  case 2: 
    p=&(fluidattr[2]); 
    break;
  case 3: 
    p=&(fluidattr[3]); 
    break;
  case 4: 
    p=&(fluidattr[4]); 
    break;
  case 5: 
    p=&(fluidattr[5]); 
    break;
  case 6: 
    p=&(fluidattr[6]); 
    break;
  case 7: 
    p=&(fluidattr[7]); 
    break;
  case 8: 
    p=&(fluidattr[8]); 
    break;
  default: 
    fprintf(stderr,"Fluid attribute not implemented\n");
    exit(-1);
  }
#endif

  return p;
}

//------------------------------------------------------------------------------

double * Structopt_sd::getptrgradfluidattr( int& loc1) {

  double *p=0;
	  
#ifdef AEROELASTIC	  
  double *gradfluidattr = structdom->getgradFluidAttr();
	   	      
  switch (loc1) {
  case 0: 
    p=&(gradfluidattr[0]); 
    break;
  case 1: 
    p=&(gradfluidattr[1]); 
    break;
  case 2: 
    p=&(gradfluidattr[2]); 
    break;
  case 3: 
    p=&(gradfluidattr[3]); 
    break;
  case 4: 
    p=&(gradfluidattr[4]); 
    break;
  case 5: 
    p=&(gradfluidattr[5]); 
    break;
  case 6: 
    p=&(gradfluidattr[6]); 
    break;
  case 7: 
    p=&(gradfluidattr[7]); 
    break;
  case 8: 
    p=&(gradfluidattr[8]); 
    break;
  default: 
    fprintf(stderr,"Fluid attribute not implemented\n");
    exit(-1);
  }
#endif
  return p;
}

//------------------------------------------------------------------------------

double *Structopt_sd::getptrelecattr(int &num,int &loc1,int &loc2) 
{

  int pos = 3*(numElecAttrVar);
  double *p = 0;

  // initialize array to store electrostatic attributes and gradients

  if (!elecattr) elecattr = new double[3*optpro->numElecStcVar];

  // initialize positions in resize arrays before obtaining pointer to them

  elecattr[pos]   = 0.0;
  elecattr[pos+1] = double (loc1);
  elecattr[pos+2] = double (loc2);

  elecattrList[num-1] = numElecAttrVar;

#ifdef AEROELASTIC	  
   
  switch (loc2) {
  case 1:   // elastic modulus (used for weighting of efield)
  case 3:   // density
  case 5:   // permittivity
  case 100: // voltage
    p = &(elecattr[pos]); 
    break;
  default: 
    fprintf(stderr,"Electrostatic attribute %d not implemented\n\n",loc2);
    exit(-1);
  }

#endif // AEROELASTIC

  numElecAttrVar++;
 
  return p;

}

//------------------------------------------------------------------------------

double *Structopt_sd::getptrgradelecattr(int &num,int &loc2) 
{

  double *p = 0;

  // initialize array to store electrostatic attributes and gradients

  if (!gradelecattr) gradelecattr = new double[optpro->numElecStcVar];

  gradelecattr[numGradElecAttrVar] = 0.0;

#ifdef AEROELASTIC	  
   
  switch (loc2) {
  case 1:   // elastic modulus (used for weighting of efield)
  case 3:   // density
  case 5:   // permittivity
  case 100: // voltage
    p = &(gradelecattr[numGradElecAttrVar]); 
    break;
  default: 
    fprintf(stderr,"Electrostatic attribute %d not implemented\n\n",loc2);
    exit(-1);
  }

#endif // AEROELASTIC

  numGradElecAttrVar++;
  return p;

}

//------------------------------------------------------------------------------

double *Structopt_sd::getptrthermattr(int &num,int &loc1,int &loc2) 
{

  int pos = 3*(numThermAttrVar);
  double *p = 0;

  // initialize array to store thermal attributes and gradients

  if (!thermattr) thermattr = new double[3*optpro->numThermStcVar];

  // initialize positions in resize arrays before obtaining pointer to them

  thermattr[pos]   = 0.0;
  thermattr[pos+1] = double (loc1);
  thermattr[pos+2] = double (loc2);

  thermattrList[num-1] = numThermAttrVar;

#ifdef AEROELASTIC	  
   
  switch (loc2) {
  case 2: //radiation
    p = &(thermattr[pos]);
    break;
  case 4: // convectivity
    p = &(thermattr[pos]); 
    break;
  case 5: // conductivity
    p = &(thermattr[pos]); 
    break;
  default: 
    fprintf(stderr,"Thermal attribute %d not implemented\n\n",loc2);
    exit(-1);
  }

#endif // AEROELASTIC

  numThermAttrVar++;
 
  return p;

}

//------------------------------------------------------------------------------

double *Structopt_sd::getptrgradthermattr(int &num,int &loc2) 
{

  double *p = 0;

  // initialize array to store thermal attributes and gradients

  if (!gradthermattr) gradthermattr = new double[optpro->numThermStcVar];

  gradthermattr[numGradThermAttrVar] = 0.0;

#ifdef AEROELASTIC	  
   
  switch (loc2) {
  case 2: // radiation
    p = &(gradthermattr[numGradThermAttrVar]); 
    break;
  case 4: // convectivity
    p = &(gradthermattr[numGradThermAttrVar]); 
    break;
  case 5: // conductivity
    p = &(gradthermattr[numGradThermAttrVar]); 
    break;
  default: 
    fprintf(stderr,"Thermal attribute %d not implemented\n\n",loc2);
    exit(-1);
  }

#endif // AEROELASTIC

  numGradThermAttrVar++;
  return p;

}

//------------------------------------------------------------------------------

void Structopt_sd::func() {

  int ifail=1;

  double funcTime = getTime();

  // first solve reliability analysis problem

  if (reliabilityFlag) {

    // need to set failcritOnly=0 otherwise analytic SA in reliability module
    // is not working (skipping to evaluate SA)

    double raTime = getTime();

    int fcoSave  = failcritOnly;
    failcritOnly = 0;

    reliabilityProb->solve(reliabilityStrc);

    failcritOnly = fcoSave;

    fprintf(stderr,"\n ... Time spent in reliability analysis: %e sec\n\n",
            (getTime()-raTime)/1000.0);
  }

  // now solve remaining analysis problems

  if ( ! failcritOnly ) {

#ifdef AEROELASTIC
    sndOptpar(-1,-1);
#endif
    
    int ina;
    for (ina=0;ina<numAnalysis;ina++) {
 
      structdom->activateAnalysis(ina);

      int apb = analysisMap[ina];

      switch (structdom->probType()) {
      case SolverInfo::Static:
	if(structdom->isComplex())
	  {
	    c_staticsolvers[apb]->optSolve();	//Static Analysis
	    ifail=0;
	  }
	else
	  {
	    staticsolvers[apb]->optSolve();	//Static Analysis
	    ifail=0;
	  }
	break;
	/*
      case SolverInfo::NonLinStatic:
	nonlinsolvers[apb]->optSolve();       //Nonlinear Analysis
	ifail=0;
	break;
      case SolverInfo::Modal:
	eigensolvers[apb]->optSolve();	//Eigenvalue Analysis
	ifail=0;
	break;
      case SolverInfo::Dynamic:
	dynsolvers[apb]->optSolve(this);      //Dynamic Analysis
	ifail=0;
	break;
      case SolverInfo::NonLinDynam:           //Nonlinear Dynamic Analysis
	NLdynsolvers[apb]->optSolve(this);
	ifail=0;
	break;
      case SolverInfo::FVibr:           // Forced Vibrations
	fvibrsolvers[apb]->optSolve();
	ifail=0;
	break;
	*/
      default:
	assert(0);
      }
    }
 
    if (ifail) {
      fprintf(stderr,"Problem Type not implemented in Optimization Module\n");
      exit(-1);
    }
  }

  // Evaluation of Opt.Crit

  evaluate();                           

  // Count function calls   

  numFunc++;                             

  // Output CPU time for total function evaluation

  fprintf(stderr,"\n ... Time spent in function evaluation: %e sec\n\n",
          (getTime()-funcTime)/1000.0);
}

//------------------------------------------------------------------------------

void Structopt_sd::evaluate() {

  // get optimization results from fluid

#ifdef AEROELASTIC

  // save aerodynamic forces of fluid solution
  
  if ( ! failcritOnly ) structdom->saveExtDomainRes();

#endif

  int i;
  for (i=0;i<optpro->numcrit;i++) { 

    optpro->opc[i]->evaluate(this);
  }
}

//------------------------------------------------------------------------------

void Structopt_sd::postProcessing(int giter ) {

  double time=giter;

  int ifail=1;

  int ina;
  for (ina=0;ina<numAnalysis;ina++) {
 
    structdom->activateAnalysis(ina);  
    
    int apb = analysisMap[ina];

    switch (structdom->probType()) {
    case SolverInfo::Static:
      if(structdom->isComplex())
	{
	  c_staticsolvers[apb]->optPostProcessor(time);
	  ifail=0;
	}
      else
	{
	  staticsolvers[apb]->optPostProcessor(time);
	  ifail=0;
	}	
      break;
      /*
    case SolverInfo::NonLinStatic:
      nonlinsolvers[apb]->optPostProcessor(time);
      ifail=0;
      break;
    case SolverInfo::Modal:
      eigensolvers[apb]->optPostProcessor(time);
      ifail=0;
      break;
    case SolverInfo::Dynamic:
      //dynsolvers[apb]->optPostProcessor(time);
      ifail=0;
      break;
    case SolverInfo::NonLinDynam:
      NLdynsolvers[apb]->optPostProcessor(time);
      ifail=0;
      break;
    case SolverInfo::FVibr:
      fvibrsolvers[apb]->optPostProcessor(time);
      ifail=0;
      break;
      */
    default:
      assert(0);
    }
  }

  if (ifail) {
    fprintf(stderr,"Problem Type not implemented in Optimization Module");
    exit(-1);
  }
}

//------------------------------------------------------------------------------

int Structopt_sd::getAnalysisData(int anaId, double& time)
{
  // initialize 

  bcx = 0;
  
  sol = 0;
  vel = 0;
  acc = 0;

  gs = 0;
  aC = 0;
  
  lam      = 1.0;
  lam_anal = 1.0;

  gStiff   = 0;
  kelArray = 0;

  gStiff = 0;
  gMass  = 0;

  kelArray = 0;

  cbcx = 0;
  csol = 0;

  // activate problem type and get pointers

  structdom->activateAnalysis(anaId);  

  int apb = analysisMap[anaId];
  
  int lc;
  
  switch (structdom->probType()) {

  case SolverInfo::Static:
    if(domain->isComplex())
      {
	lc        = (time < 0) ? 0 : static_cast<int>(time);
	csol      = c_staticsolvers[apb]->getpsol(lc);
	cbcx      = c_staticpros[apb]->getbc();
	kelArray  = c_staticpros[apb]->getkelArray();
      }
    else
      {
	lc       = (time < 0) ? 0 : static_cast<int>(time);
	sol      = staticsolvers[apb]->getpsol(lc);
	bcx      = staticpros[apb]->getbc();
	kelArray = staticpros[apb]->getkelArray();
      }
    break;
	/*
  case SolverInfo::NonLinStatic:
    lc       = (time < 0) ? 0 : (int)time; 
    aC       = nonlinpros[apb]->getCorotators();
    sol      = nonlinsolvers[apb]->getSol(lc);
    gs       = nonlinsolvers[apb]->getGeomState(lc);
    bcx      = nonlinpros[apb]->getbc();
    lam      = nonlinsolvers[apb]->getLambdaCrit(lc);
    if (lam == 0.0) lam      = nonlinsolvers[apb]->getLambdaLC(lc);    
    break;
  case SolverInfo::Dynamic:
    sol      = dynsolvers[apb]->getpDis();
    vel      = dynsolvers[apb]->getpVel();
    acc      = dynsolvers[apb]->getpAcc();
    gStiff   = dynpros[apb]->getpK(dynsolvers[apb]->getpOps());
    gMass    = dynpros[apb]->getpM(dynsolvers[apb]->getpOps());
    gDamp    = dynpros[apb]->getpC(dynsolvers[apb]->getpOps());
    bcx      = dynpros[apb]->boundaryValue();
    break;
  case SolverInfo::NonLinDynam:
    gs       = NLdynsolvers[apb]->getGeomState();
    aC       = NLdynpros[apb]->getCorotators();
    sol      = NLdynsolvers[apb]->getSol(); 
    bcx      = NLdynpros[apb]->boundaryValue();
    lam      = NLdynpros[apb]->getconvLambda();
    break;
  case SolverInfo::FVibr:
    lc       = (time < 0) ? 0 : static_cast<int>(time);
    csol     = fvibrsolvers[apb]->getpsol(lc);
    cbcx     = fvibrpros[apb]->getbc();
    kelArray = 0;
    break;
    */
  default:
    fprintf(stderr,"Error: Solver Type not implemented in Optimization\n");
    exit(-1);
  }

  return structdom->probType();
}

//------------------------------------------------------------------------------

int Structopt_sd::getAnalysisGradData(int anaId, double& time)
{
  // initialize 

  bcx = 0;
  
  sol  = 0;
  vel  = 0;
  grad = 0;
  adj  = 0;

  gs = 0;
  aC = 0;

  gStiff   = 0;
  kelArray = 0;

  gStiff = 0;
  gMass  = 0;

  kelArray = 0;

  cbcx = 0;
  csol = 0;
  cgrad = 0;
  cadj = 0;

  // activate problem type and get pointers

  structdom->activateAnalysis(anaId);  

  int apb = analysisMap[anaId];
  
  int lc;
  
  switch (structdom->probType()) {

  case SolverInfo::Static:
    if(structdom->isComplex())
      {
	lc       = (time < 0) ? 0 : static_cast<int>(time);
	csol	 = c_staticsolvers[apb]->getpsol(lc);
	cgrad    = c_staticsolvers[apb]->getpgrad(lc);
	cbcx	 = c_staticpros[apb]->getbc();
	kelArray = c_staticpros[apb]->getkelArray();
      }
    else
      {
	lc       = (time < 0) ? 0 : static_cast<int>(time);
	sol	 = staticsolvers[apb]->getpsol(lc);
	grad     = staticsolvers[apb]->getpgrad(lc);
	bcx	 = staticpros[apb]->getbc();
	kelArray = staticpros[apb]->getkelArray();
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
    lc	 = (time < 0) ? 0 : (int)time;
    gs	 = nonlinsolvers[apb]->getGeomState(lc);
    aC	 = nonlinpros[apb]->getCorotators();
    sol      = nonlinsolvers[apb]->getSol(lc);
    grad     = nonlinsolvers[apb]->getGradSol(lc);
    bcx      = nonlinpros[apb]->getbc();
    break;
  case SolverInfo::Dynamic:
    sol      = dynsolvers[apb]->getpDis();
    grad	 = dynsolvers[apb]->getpgDis();
    bcx      = dynpros[apb]->boundaryValue();
    break;
  case SolverInfo::NonLinDynam:
    gs       = NLdynsolvers[apb]->getGeomState();
    aC       = NLdynpros[apb]->getCorotators();
    sol      = NLdynsolvers[apb]->getSol(); 
    grad	 = NLdynsolvers[apb]->getpgDis();
    bcx      = NLdynpros[apb]->boundaryValue();
    break;
  case SolverInfo::FVibr:
    lc       = (time < 0) ? 0 : static_cast<int>(time);
    csol     = fvibrsolvers[apb]->getpsol(lc);
    cgrad    = fvibrsolvers[apb]->getpgrad(lc);
    cbcx     = fvibrpros[apb]->getbc();
    kelArray = 0;
    break;
    */
  default:
    fprintf(stderr,"Error: Solver Type not implemented in Optimization\n");
    exit(-1);
        
  }

  return structdom->probType();
}

//------------------------------------------------------------------------------

int Structopt_sd::getAnalysisAdjData(int anaId, double& time, int zeroAll)
{
  // initialize 
  bcx = 0;
  
  sol  = 0;
  vel  = 0;
  grad = 0;
  adj  = 0;

  gs = 0;
  aC = 0;

  gStiff   = 0;
  kelArray = 0;

  gStiff = 0;
  gMass  = 0;

  kelArray = 0;
  
  dqdlam = 0.0;
  
  csol  = 0;
  cbcx  = 0;
  cgrad = 0;
  cadj  = 0;

  // zero all adjoint vectors
  
  if (zeroAll) {
  
    int ina;
    for (ina=0;ina<numAnalysis;ina++) {

      adj = 0;
      cadj = 0;

      structdom->activateAnalysis(ina);  
    
      int apb = analysisMap[anaId];

      switch (structdom->probType()) {
      case SolverInfo::Static:
	if(structdom->isComplex())
	  {
	    cadj   = c_staticsolvers[apb]->getpadj();
	  }
	else
	  {
	    adj   = staticsolvers[apb]->getpadj();
	  }
	break;
	/*
      case SolverInfo::NonLinStatic:
	adj   = nonlinsolvers[apb]->getAdj();
	break; 	
      case SolverInfo::Dynamic:
	adj   = dynsolvers[apb]->getpadj();
	break;
      case SolverInfo::NonLinDynam:
	adj   = NLdynsolvers[apb]->getpadj();
	break;
      case SolverInfo::FVibr:
	cadj   = fvibrsolvers[apb]->getpadj();
	break;
	*/
      default:
	assert(0);
      }
    
      if (adj) adj->zero();
      if (cadj) cadj->zero();
    }
  }

  // activate problem type and get pointers

  structdom->activateAnalysis(anaId);  

  int apb = analysisMap[anaId];
  
  int lc;
  
  switch (structdom->probType()) {

  case SolverInfo::Static:
    if(structdom->isComplex())
      {
	lc       = (time < 0) ? 0 : static_cast<int>(time);
	csol	 = c_staticsolvers[apb]->getpsol(lc);
	cadj   	 = c_staticsolvers[apb]->getpadj();
	cbcx	 = c_staticpros[apb]->getbc();
	kelArray = c_staticpros[apb]->getkelArray();
      }
    else
      {
	lc       = (time < 0) ? 0 : static_cast<int>(time);
	sol	 = staticsolvers[apb]->getpsol(lc);
	adj   	 = staticsolvers[apb]->getpadj();
	bcx	 = staticpros[apb]->getbc();
	kelArray = staticpros[apb]->getkelArray();
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
    lc	 = (time < 0) ? 0 : (int)time;
    gs	 = nonlinsolvers[apb]->getGeomState(lc,1);
    aC	 = nonlinpros[apb]->getCorotators();
    sol      = nonlinsolvers[apb]->getSol(lc);
    adj      = nonlinsolvers[apb]->getAdj();
    bcx      = nonlinpros[apb]->getbc();
    kelArray = nonlinpros[apb]->getkelArray();
    break;
  case SolverInfo::Dynamic:
    sol      = dynsolvers[apb]->getpDis();
    adj 	 = dynsolvers[apb]->getpadj();
    bcx      = dynpros[apb]->boundaryValue();
    break;
  case SolverInfo::NonLinDynam:
    gs       = NLdynsolvers[apb]->getGeomState();
    aC       = NLdynpros[apb]->getCorotators();
    sol      = NLdynsolvers[apb]->getSol(); 
    adj	 = NLdynsolvers[apb]->getpadj();
    bcx      = NLdynpros[apb]->boundaryValue();
    kelArray = NLdynpros[apb]->getkelArray();
    break;
  case SolverInfo::FVibr:
    lc       = (time < 0) ? 0 : static_cast<int>(time);
    csol	 = fvibrsolvers[apb]->getpsol(lc);
    cadj   	 = fvibrsolvers[apb]->getpadj();
    cbcx	 = fvibrpros[apb]->getbc();
    kelArray = 0;
    break;
    */
  default:
    fprintf(stderr,"Error: Solver Type not implemented in Optimization\n");
    exit(-1);
  }

  
  return structdom->probType();
}

//------------------------------------------------------------------------------

double Structopt_sd::getstrainenergy(int* eleList, int listSize, double& time,
				     int anaId) 
{
  
  double val;  

  int probType = getAnalysisData(anaId,time);
  
  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val = structdom->getStrainEnergy(*sol,bcx,gStiff,kelArray,eleList,
					 listSize);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    val = structdom->getNLStrainEnergy(gs,aC,eleList,listSize);
    break;
    */
  default:
    fprintf(stderr,"Error: strain energy can only be evaluated for\n");
    fprintf(stderr,"       static/dynamic and linear/nonlinear problems\n");
    exit(-1);
  }

  return val;

}

//------------------------------------------------------------------------------
double Structopt_sd::getgradstrainenergy(int* eleList, int listSize, double& time,
                                      int anaId) 
{
  double val;
  int probType = getAnalysisGradData(anaId, time);
  switch (probType) 
    {      
    case SolverInfo::Static:
    case SolverInfo::Dynamic:     
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val = structdom->getGradStrainEnergy(*sol,*grad,bcx,kelArray,eleList,
					     listSize);
      }
      break;
      /*
    case SolverInfo::NonLinStatic:
    case SolverInfo::NonLinDynam:
      val = structdom->getGradNLStrainEnergy(gs,*grad, aC,eleList,listSize);
      break;
      */
    default:
      fprintf(stderr,"Error: derivatives of strain energy can only be\n");
      fprintf(stderr,"       evaluated for (NL)static and dynamic problems\n");
      fprintf(stderr,"...... stop\n");
      exit(-1);
    }  
  return val;
} 

//------------------------------------------------------------------------------
double Structopt_sd::getgradpartstrainenergy(int* eleList, int listSize, 
                                          double& time, int anaId) {

  double val;

  int probType = getAnalysisData(anaId,time);

  switch (structdom->probType()) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val=structdom->getGradPartStrainEnergy(*sol,bcx,kelArray,eleList,
					       listSize);
      }
    break;
    /*
  case SolverInfo::NonLinDynam:
  case SolverInfo::NonLinStatic:
    val=structdom->getGradNLPartStrainEnergy(gs,aC,kelArray,eleList,listSize);
    break;
    */
  default:
    fprintf(stderr,"Error: derivatives of strain energy can only be \n");
    fprintf(stderr,"	     evaluated for static and dynamic problems \n");
    exit(-1);
  }

  return val;

} 

//------------------------------------------------------------------------------

void Structopt_sd::getgraddustrainenergy(int* eleList, int listSize, double& time,
                                      int anaId) 
{
     
  int probType = getAnalysisAdjData(anaId, time);

  switch (structdom->probType()) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {    
	adj->zero();  
	structdom->getGradduStrainEnergy(*sol,*adj,bcx,kelArray,eleList,
					 listSize);
      }
    break;
    /*
  case SolverInfo::NonLinDynam:
  case SolverInfo::NonLinStatic:
    adj->zero();  
    structdom->getGradNLduStrainEnergy(gs,*adj,aC,kelArray,eleList,
				       listSize);
    break;
    */
  default:
    fprintf(stderr,"Error: derivatives of strain energy can only be \n");
    fprintf(stderr,"	     evaluated for static and dynamic problems \n");
    exit(-1);
  }

}

//------------------------------------------------------------------------------

double Structopt_sd::getmass(int* eleList, int listSize, int anaId) 
{
  
  double val = structdom->getStructureMass(eleList,listSize);
  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradmass(int* eleList, int listSize, int anaId) 
{
  double val = structdom->getGradStructureMass(eleList,listSize);          
  return val;
}

//------------------------------------------------------------------------------

void Structopt_sd::getgraddumass(int* eleList, int listSize, int anaId) 
{
  double time = 0.0;
  int probType = getAnalysisAdjData(anaId,time);
  if(structdom->isComplex())
    {
      cadj->zero();
    }
  else
    {
      adj->zero();
    }
}

//------------------------------------------------------------------------------

double Structopt_sd::getMomOfInertia(int* eleList, int listSize, int momFlag,
                                  int anaId)
{
     
  if (momFlag < 3 || momFlag > 5) {
    fprintf(stderr," Error: wrong reference axis for moment of inertia\n");
    fprintf(stderr,"       stop\n");
    exit(-1);
  }
  double val = structdom->getMomOfInertia(eleList,listSize,momFlag);          
  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradMomOfInertia(int* eleList, int listSize, int momFlag,
                                      int anaId)
{
  double val = structdom->getGradMomOfInertia(eleList,listSize,momFlag);          
  return val;
}

//------------------------------------------------------------------------------

void Structopt_sd::getgradduMomOfInertia(int* eleList, int listSize, int momFlag, 
                                      int anaId) 
{
  double time = 0.0; 

  int probType = getAnalysisAdjData(anaId,time);

  adj->zero();
}

//------------------------------------------------------------------------------

double Structopt_sd::getfrequency(int ieig, int anaId) 
{
          
  Vector * eigval;

  int apb = analysisMap[anaId];

  structdom->activateAnalysis(anaId);  
  /*
  if (structdom->probType() == SolverInfo::Modal) {
    eigval = eigensolvers[apb]->getpeigval();
  }
  else {
  */
    fprintf(stderr,"Error: eigenfrequencies can only be evaluated\n");
    fprintf(stderr,"	    for eigen problems ...... stop\n");
    exit(-1);
    /*}*/

  double val = sqrt( (*eigval)[ieig] ) / (2.0*3.14159265);     
  
  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradfrequency(int ieig, int anaId ) 
{

  Vector * eigval;
  Vector * eigrad;

  int apb = analysisMap[anaId];

  structdom->activateAnalysis(anaId);  
  
  /*
  if (structdom->probType() == SolverInfo::Modal) {
    eigval = eigensolvers[apb]->getpeigval();
    eigrad = eigensolvers[apb]->getpeigrad();
  }
  else {
  */
    fprintf(stderr,"Error: derivatives of eigenfrequencies can only be\n");
    fprintf(stderr,"       evaluated for eigen problems ...... stop\n");
    exit(-1);
    /*}*/

  double grad = (*eigrad)[ieig] / ( 4.0*3.14159265*sqrt( (*eigval)[ieig] ) );
    
  return grad;

}
    
//------------------------------------------------------------------------------

double Structopt_sd::getcontrolcost(int anaId) 
{
          
  double val;

  int apb = analysisMap[anaId];

  structdom->activateAnalysis(anaId);  
  
  /*
  if (structdom->probType() == SolverInfo::Modal) {
    val = eigensolvers[apb]->getControlCost();
  }
  else {  
  */
    fprintf(stderr,"Error: Control cost can only be evaluated\n");
    fprintf(stderr,"       for eigen problems ...... stop\n");
    exit(-1);
    /*}*/

  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getkineticenergy(double& time, int anaId) 
{

  double val;

  int probType = getAnalysisData(anaId,time);
     
  switch (probType) {
    /*
  case SolverInfo::Dynamic:
    val=structdom->getKineticEnergy(*vel, gMass);
    break;
    */
  default:
    fprintf(stderr,
	    "Error: kinetic energy make sens for dynamic problems only. stop\n");
    exit(-1);
  }

  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradkineticenergy(double& time, int anaId) 
{

  double val=0.0;

  fprintf(stderr,"Error: analytical derivatives of kinetic energy are not\n");
  fprintf(stderr,"	 implemented yet .... stop\n");
  exit(-1);
     
  return val;

} 

//------------------------------------------------------------------------------

double Structopt_sd::getdampingenergy(double& time, int anaId) 
{

  double val;

  int probType = getAnalysisData(anaId,time);
     
  switch (probType) {
    
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val=structdom->getDampingEnergy(*vel, gDamp);
      }
    break;
  default:
    fprintf(stderr,
	    "Error: kinetic energy make sens for dynamic problems only. stop\n");
    exit(-1);
  }

  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgraddampingenergy(double& time, int anaId) {

  double val=0.0;

  fprintf(stderr,"Error: analytical derivatives of damping energy are not\n");
  fprintf(stderr,"	 implemented yet .... stop\n");
  exit(-1);
     
  return val;

} 

//----------------------------------------------------------------------

double Structopt_sd::getExtVal(int crtNum, int ind, int anaId, double& time)
{ 
 
  double val = 0.0;

  int probType = getAnalysisData(anaId,time);

  switch (probType) {
    /*
  case SolverInfo::Dynamic:
  case SolverInfo::NonLinDynam:
    val=0;
    break;
    */
  default:
    fprintf(stderr,"Error: Aeroforces only for dynamic problems in aeroelasticity.\n"); 
    fprintf(stderr,"Stop in getaeroforce - 1 !\n");
    exit(-1);
  }

#ifdef AEROELASTIC
  val=structdom->getExtVal(crtNum,ind);
#endif
 
  return val;

}


//------------------------------------------------------------------------------

double Structopt_sd::getaeroforce(int idir, double& time, int anaId)
{

  double val;

  int probType = getAnalysisData(anaId,time);

  switch (probType) {
    /*
  case SolverInfo::Dynamic:
  case SolverInfo::NonLinDynam:
    val=0;
    break;
    */
  default:
    fprintf(stderr,"Error: Aeroforces only for dynamic problems in aeroelasticity.\n"); 
    fprintf(stderr,"Stop in getaeroforce - 1 !\n");
    exit(-1);
  }

#ifdef AEROELASTIC
  val=structdom->getAeroForce(idir);
#else
  fprintf(stderr,"Error: Aeroforces only for dynamic problems in aeroelasticity.\n"); 
  fprintf(stderr,"Stop in getaeroforce - 2 !\n");
  exit(-1);
#endif      

  return val;

} 

//------------------------------------------------------------------------------

double Structopt_sd::getgradaeroforce(int idir, double& time, int anaId) 
{

  double val;

  int probType = getAnalysisData(anaId,time);

  switch (probType) {
    /*
  case SolverInfo::Dynamic:
  case SolverInfo::NonLinDynam:
    val=0;
    break;
    */
  default:
    fprintf(stderr,"Error: Aeroforces only for dynamic problems in aeroelasticity.\n"); 
    fprintf(stderr,"Stop in getaeroforce - 1 !\n");
    exit(-1);
  }

#ifdef AEROELASTIC
  if (anagrdType == 1)
    val=structdom->getgradAeroForce(idir);
#else
  fprintf(stderr,"Error: Derivative of Aeroforces only for dynamic problems in aeroelasticity.\n"); 
  fprintf(stderr," Stop in getgradaeroforce - 2 !\n");
  exit(-1);
#endif      
     
  return val;

} 

//------------------------------------------------------------------------------

void Structopt_sd::getgradduaeroforce(int idir, double& time, int anaId)
{

  // purpose: zero adj

  int probType = getAnalysisAdjData(anaId,time);

  switch (probType) {
    /*
  case SolverInfo::Dynamic:
  case SolverInfo::NonLinDynam:
    adj->zero();  
    break;
    */
  default:
    fprintf(stderr,"Error: Aeroforces only for dynamic problems in aeroelasticity.\n"); 
    fprintf(stderr," Stop in getgradduaeroforce - 1 !\n");
    exit(-1);
  }

} 

//-------------------------------------------------------------------------

void Structopt_sd::getgradduExtVal(int crtNum, int anaId, double& time)
{

  // purpose: zero adj

  int probType = getAnalysisAdjData(anaId,time);

  switch (probType) {
    /*
  case SolverInfo::Dynamic:
  case SolverInfo::NonLinDynam:
    adj->zero();  
    break;
    */
  default:
    fprintf(stderr,"Error: External Criteria only for dynamic problems.\n"); 
    fprintf(stderr," Stop in getgradduaeroforce - 1 !\n");
    exit(-1);
  }

} 

//------------------------------------------------------------------------------

double Structopt_sd::getestat(int idir,int *eleList,int listSize,double &time,
                           int anaId ) 
{

  double val=0.0;

  int probType = getAnalysisData(anaId,time);

  switch (probType) {

  case SolverInfo::Static:

    switch (idir) {
      /*
    case 0: // sum of eforces
      structdom->getElectricField(*sol,bcx,-1,3,0);
      structdom->getElecForceOnStruc(-1,5,0);
      val = structdom->getSumElecForce();
      break;
      */
    case 1:
      if(structdom->isComplex())
	{
	  assert(0);
	}
      else
	{
	  val = structdom->getStrainEnergy(*sol,bcx,gStiff,kelArray,eleList,listSize);
	}
      break;
    default:
      fprintf(stderr,"\nElectrostatic criteria %d not implemented. Exiting!\n",idir);
      exit(-1);
      break;
    }
    break;
    /*
  case SolverInfo::Dynamic:

    switch (idir) {
#ifdef AEROELASTIC
    case 0: // sum of eforces
      break;
    case 1:
      break;
#endif // AEROELASTIC
   
    default:
      fprintf(stderr,"\nElectrostatic criteria %d not implemented. Exiting!\n",idir);
      exit(-1);
      break;
    }
    break; 
    */
  default:
    assert(0);
  }

  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradestat(int idir,int *eleList,int listSize,double &time,
                               int anaId) 
{

  double val=0.0;

  int probType = getAnalysisGradData(anaId, time);

  switch (probType) {

  case SolverInfo::Static:

    switch (idir) {
      /*
    case 0:
      // need gradients of electric field for electric force gradients
      structdom->getGradElectricField(*sol,*grad,bcx);
      val = structdom->getGradElectricForce();
      break;
      */
    case 1:
      if(structdom->isComplex())
	{
	  assert(0);
	}
      else
	{
	  val = structdom->getGradStrainEnergy(*sol,*grad,bcx,kelArray,eleList,listSize);
	}
    break;
    default:
      fprintf(stderr," \n Total grad of electrostatic criteria %d not implemented. Exiting!\n\n",idir);
      exit(-1);
    }
    break;
    /*
  case SolverInfo::Dynamic:
    switch (idir) {
#ifdef AEROELASTIC
    case 0: // sum of eforces
      break;
    case 1:
      break;
#endif // AEROELASTIC
    default:
      fprintf(stderr,"\n Total grad of electrostatic criteria %d not implemented. Exiting!\n",idir);
      exit(-1);
      break;
    }
    break; 
    */
  default:
    assert(0);
  }

  return val;
} 

//------------------------------------------------------------------------------

double Structopt_sd::getgradpartestat(int idir,int *eleList,int listSize,
                                   double &time,int anaId) 
{

  double val=0.0;

  switch (structdom->probType()) {

  case SolverInfo::Static:

    switch (idir) {
    case 1:
      if(structdom->isComplex())
      {
	assert(0);
      }
      else
	{
	  val = structdom->getGradPartStrainEnergy(*sol,bcx,kelArray,
						   eleList,listSize);
	}
      break;
    default:
      fprintf(stderr," \n Partial grad of electrostatic criteria %d not implemented. Exiting!\n\n",idir);
      exit(-1);
    }
    break;
    /*
  case SolverInfo::Dynamic:
    switch (idir) {
#ifdef AEROELASTIC
    case 0: // sum of eforces
      break;
    case 1:
      break;
#endif // AEROELASTIC
    */
    default:
      fprintf(stderr,"\n Partial grad of electrostatic criteria %d not implemented. Exiting!\n",idir);
      exit(-1);
      break;
  }

  return val;
} 

//------------------------------------------------------------------------------

double Structopt_sd::getDCmass(double& powFac,int* eleList, int listSize, double& time, int anaId) 
{
     
  double val;

  int probType = getAnalysisData(anaId, time);
    
  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val = structdom->getDCmass(*sol, bcx, powFac, eleList, listSize);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    */
  default:
    fprintf(stderr,"Error: mass in deformed configuration can only \n");
    fprintf(stderr,"       be evaluated for linear static problems \n");
    fprintf(stderr,"       ..... stop\n");
    exit(-1);
  }
 
  return val;
}

//------------------------------------------------------------------------------

double Structopt_sd::getgradDCmass(double& powFac, int* eleList, int listSize, double& time, int anaId) 
{

  double val;

  int probType = getAnalysisGradData(anaId, time);

  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val = structdom->getGradDCmass(*sol, *grad, bcx, powFac, eleList, listSize);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    */
  default: 
    fprintf(stderr,"Error: derivative of mass in deformed configuration \n");
    fprintf(stderr,"       can only be evaluated for linear static problems \n");
    fprintf(stderr,"       ..... stop\n");
    exit(-1);
  }

  return val;
}

//------------------------------------------------------------------------------

double Structopt_sd::getgradpartDCmass(double& powFac, int* eleList, int listSize, double& time, int anaId) 
{
  double val;

  int probType = getAnalysisData(anaId, time);
    
  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val = structdom->getGradPartDCmass(*sol, bcx, powFac, eleList, listSize);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    */
  default:
    fprintf(stderr,"Error: mass in deformed configuration can only \n");
    fprintf(stderr,"       be evaluated for linear static problems \n");
    fprintf(stderr,"       ..... stop\n");
    exit(-1);
  }
 
  return val;
}

//------------------------------------------------------------------------------

void Structopt_sd::getgradduDCmass(double& powFac, int* eleList, int listSize, double& time, int anaId) 
{
  int probType = getAnalysisAdjData(anaId, time);

  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	adj->zero();
	structdom->getGradduDCmass(*sol, *adj, bcx, powFac, eleList, listSize);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    */
  default:
    fprintf(stderr,"Error: derivatives of nodal stress can only be\n");
    fprintf(stderr,"       evaluated for static and dynamic problems\n");
    fprintf(stderr,"       ... stop\n");
    exit(-1);
  }
}

//------------------------------------------------------------------------------

double Structopt_sd::getnodstr(int inode, int typ, double& time, int anaId) 
{

  int surface,stressIndex;
  getStressInfo(typ,surface,stressIndex);  

  double val;

  int probType = getAnalysisData(anaId, time);
    
  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val   = structdom->getNodalStressStrain(*sol, bcx, inode,
						stressIndex, surface);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    val  = structdom->getNLNodalStressStrain(gs, aC, inode,
					     stressIndex, surface);
    break;
  case SolverInfo::FVibr:
    val   = structdom->getNodalStressStrain(*csol, cbcx, inode,
					    stressIndex, surface);
    break;
    */
  default:
    fprintf(stderr,"Error: nodal stress can only be evaluated for\n");
    fprintf(stderr,"       static and dynamic problems ..... stop\n");
    exit(-1);
  }
 
  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradnodstr(int inode, int typ, double& time, int anaId) 
{

  int surface,stressIndex;
  getStressInfo(typ,surface,stressIndex);  
  
  double val;

  int probType = getAnalysisGradData(anaId, time);

  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val   = structdom->getGradNodalStressStrain(*sol, *grad, bcx, inode,
						    stressIndex, surface);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    val   = structdom->getNLGradNodalStressStrain(gs,aC,*grad,
						  inode,stressIndex,surface);
    break;
  case SolverInfo::FVibr:
    val   = structdom->getGradNodalStressStrain(*csol, *cgrad, cbcx, inode,
						stressIndex, surface);
    break;
    */
  default: 
    fprintf(stderr,"Error: derivatives of nodal stress can only be\n");
    fprintf(stderr,"       evaluated for static and dynamic problems\n");
    fprintf(stderr,"       .... stop\n");
    exit(-1);
  }

  return val;

}

//------------------------------------------------------------------------------
double Structopt_sd::getgradpartnodstr(int inode, int typ, double& time, int anaId) {

  int surface,stressIndex;
  getStressInfo(typ,surface,stressIndex);  

  double val;

  int probType = getAnalysisData(anaId,time);
  
  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val   = structdom->getGradPartNodalStressStrain(*sol, bcx, inode,
							stressIndex, surface);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    val=structdom->getNLGradPartNodalStressStrain(gs, aC, inode,
						  stressIndex, surface);
    break;
  case SolverInfo::FVibr:
    val   = structdom->getGradPartNodalStressStrain(*csol, cbcx, inode,
						    stressIndex, surface);
    break;
    */
  default:
    fprintf(stderr,"Error: derivatives of nodal stress can only be\n");
    fprintf(stderr,"       evaluated for static and dynamic problems\n");
    fprintf(stderr,"       .... stop\n");
    exit(-1);
  }

  return val;

}

//------------------------------------------------------------------------------
void Structopt_sd::getgraddunodstr(int inode, int typ, double& time, int anaId) 
{

  int surface,stressIndex;
  getStressInfo(typ,surface,stressIndex);  
    
  int probType = getAnalysisAdjData(anaId, time);

  switch (probType) {
  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	adj->zero();
	structdom->getGradduNodalStressStrain(*sol, *adj, bcx,
					      inode, stressIndex, surface);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    adj->zero();
    structdom->getNLGradduNodalStressStrain(gs, *adj, aC, inode,
					    stressIndex, surface);
    break;
  case SolverInfo::FVibr:
    cadj->zero();
    structdom->getGradduNodalStressStrain(*csol, *cadj, cbcx,
					  inode, stressIndex, surface);
    break;
    */
  default:
    fprintf(stderr,"Error: derivatives of nodal stress can only be\n");
    fprintf(stderr,"       evaluated for static and dynamic problems\n");
    fprintf(stderr,"       ... stop\n");
    exit(-1);
  }
}

//------------------------------------------------------------------------------

double Structopt_sd::getdisp( int node, int typ, int dva, double& time,
                           int anaId) 
{
    
  double val;

  int probType = getAnalysisData(anaId,time);

  int apb = analysisMap[anaId];
   
  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::NonLinStatic:
    if ( dva ) {
      fprintf(stderr,
	      "Error: no velocities or accerlations in statics .... stop !\n");
      exit(-1);
    }
    if(structdom->isComplex())
      {
	val=ScalarTypes::Real(structdom->getNodalDisp(*csol, cbcx, node, typ, true));
      }
    else
      {
	val=structdom->getNodalDisp(*sol, bcx, node, typ, true);
      }
    break;
    /*
  case SolverInfo::Dynamic:
    switch (dva) {
    case 0:
      val=structdom->getNodalDisp(*sol, bcx, node, typ, lam);
      break;
    case 1:
      val=structdom->getNodalDisp(*vel, bcx, node, typ, lam);
      break;
    case 2:
      val=structdom->getNodalDisp(*acc, bcx, node, typ, lam);
      sol	= dynsolvers[apb]->getpAcc();
    } 	 
    break;
  case SolverInfo::NonLinDynam:
    switch (dva) {
    case 0:
      val=structdom->getNodalDisp(*sol, bcx, node, typ, lam);
      break;
    case 1:
      fprintf(stderr,"Derivative of velocity not applied in NLdynamics");
      exit(-1);
      break;
    case 2:
      fprintf(stderr,"Derivative of acceleration not applied in NLdynamics");
      exit(-1);
    }
    break;
  case SolverInfo::FVibr:
    if ( dva ) 
      {
	fprintf(stderr,
		"Error: no velocities or accerlations in Forced Vibrations .... stop !\n");
	exit(-1);
      }
    val = structdom->getNodalDisp(*csol, cbcx, node, typ, lam);
    break;
    */
  default:
    fprintf(stderr,"Error: nodal displacements|velocities|accelerations can only\n");
    fprintf(stderr,"       be evaluated for static and dynamic problems ... stop\n\n");
    exit(-1);
  }

  return val;

}				       
				       
//------------------------------------------------------------------------------

double Structopt_sd::getgraddisp( int node, int typ, int dva, double& time,
                               int anaId) 
{

  if ( dva ) {
    fprintf(stderr,
	    "Error: no velocities or accerlations derivatives implemented. stop !\n");
    exit(-1);
  }	   

  double val;

  int probType = getAnalysisGradData(anaId,time);

  switch (probType) {
  case SolverInfo::Static:
  case SolverInfo::NonLinStatic:
  case SolverInfo::Dynamic:
  case SolverInfo::NonLinDynam:
    if(structdom->isComplex())
      {
	val=ScalarTypes::Real(structdom->getGradNodalDisp(*csol, *cgrad, cbcx, node, typ, true));
      }
    else
      {
	val=structdom->getGradNodalDisp(*sol, *grad, bcx, node, typ, true);
      }
    break;
    /*
  case SolverInfo::FVibr:
    val = structdom->getGradNodalDisp(*csol, *cgrad, cbcx, node, typ);
    break;
    */
  default:
    fprintf(stderr,"Error: derivatives of nodal displacements can only be\n");
    fprintf(stderr,"       evaluated for static and dynamic problems\n");
    fprintf(stderr,"...... stop\n");
    exit(-1);
  }


  return val;

}				       
				       
//------------------------------------------------------------------------------

void Structopt_sd::getgraddudisp( int node, int typ, int dva, double& time,
                               int anaId) 
{

  if ( dva ) {
    fprintf(stderr,
	    "Error: no velocities or accerlations derivatives implemented. stop !\n");
    exit(-1);
  }	   

  int probType = getAnalysisAdjData(anaId,time,1);

  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::NonLinStatic:
  case SolverInfo::Dynamic:
  case SolverInfo::NonLinDynam:
    if(structdom->isComplex())
      {
	cadj->zero();
	structdom->getGradduNodalDisp(*csol,*cadj,cbcx,node,typ, true);
      }
    else
      {
	adj->zero();
	structdom->getGradduNodalDisp(*sol,*adj,bcx,node,typ, true);
      }
    break;
    /*
  case SolverInfo::FVibr:
    structdom->getGradduNodalDisp(*csol, *cadj, cbcx, node, typ);
    break;
    */
  default:
    fprintf(stderr,"Error: derivatives of nodal displacements can only be\n");
    fprintf(stderr,"       evaluated for static and dynamic problems\n");
    fprintf(stderr,"...... stop\n");
    exit(-1);
  }
}				       

//--------------------------------------------------------------------------

double Structopt_sd::getLambda(double& time,int anaId)
{

  double val;

  int probType = getAnalysisData(anaId,time);

  switch (probType) {
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    val = lam;
    break;
    */
  default:
    fprintf(stderr,"Error: load control parameter lambda only defined \n");
    fprintf(stderr,"       for nonlinear problems.\n");
    exit(-1);
  }

  return val;
  

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradLambda(double& time, int anaId)
{

  double val;

  int probType = getAnalysisGradData(anaId,time);
  int apb = analysisMap[anaId];
  
  switch (probType) {
    /*
  case SolverInfo::NonLinStatic:
    val = nonlinsolvers[apb]->getgradLambda();
    break;
  case SolverInfo::NonLinDynam:
    val = NLdynsolvers[apb]->getgradLambda();
    break;
    */
  default:
    fprintf(stderr,"Error: load control parameter lambda only defined \n");
    fprintf(stderr,"       for nonlinear problems.\n");
    exit(-1);
  }

  return val;
  

}  
  
//------------------------------------------------------------------------------

void Structopt_sd::getgraddulambda(double& time, int anaId) 
{

  double val;

  int probType = getAnalysisAdjData(anaId,time);

  switch (probType) {
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    adj->zero();
    dqdlam = 1.0;
    break;
    */
  default:
    fprintf(stderr,"Error: load control parameter lambda only defined \n");
    fprintf(stderr,"       for nonlinear problems.\n");
    exit(-1);
  }



}  
  
//------------------------------------------------------------------------------

double Structopt_sd::getinternalforce( int inode, int typ, double& time, int anaId)
{
    
  double val;

  int probType = getAnalysisData(anaId,time);

  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val = structdom->getNodalInternalForce(*sol, bcx, inode, typ);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    val = structdom->getNodalInternalForceNL(*gs, aC, inode, typ);
    break;
    */
  default:
    fprintf(stderr,"Error: nodal internal forces can only be evaluated \n");
    fprintf(stderr,"       for static and dynamic problems.\n");
    exit(-1);
  }

  return val;

}				       

//------------------------------------------------------------------------------

double Structopt_sd::getnodalpos( int inode, int typ, double& time, int anaId)
{    
  double val = structdom->getVariationNodalPos(inode,typ);  
  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradnodalpos( int inode, int typ, double& time, int anaId)
{
  double grad = structdom->getGradVariationNodalPos(inode,typ); 
  return grad;
}

//------------------------------------------------------------------------------

double Structopt_sd::getgradinternalforce(int inode,int typ,double& time,int anaId)
{
  double val;
  int probType = getAnalysisGradData(anaId, time);
  switch (probType) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val = structdom->getGradNodalInternalForce(*sol, *grad, bcx, inode, typ);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    val = structdom->getGradNodalInternalForceNL(*gs, aC, *grad, inode, typ);
    break;
    */
  default:
    fprintf(stderr,"Error: derivatives of nodal internal forces can only be\n");
    fprintf(stderr,"       evaluated for static and dynamic problems\n");
    exit(-1);
  }

  return val;

}				       

//------------------------------------------------------------------------------

double Structopt_sd::getgradpartinternalforce(int inode,int typ,double& time,
                                           int anaId)
{

  double val;

  int probType = getAnalysisData(anaId,time);

  switch (structdom->probType()) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	val   = structdom->getGradPartNodalInternalForce(*sol, bcx, inode, typ);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    val   = structdom->getGradPartNodalInternalForceNL(*gs, aC, inode, typ);
    break;
    */
  default:
    fprintf(stderr,"Error: derivatives of nodal internal force can only be\n");
    fprintf(stderr,"       evaluated for linear static and dynamic problems.\n");
    exit(-1);
  }

  return val;

}

//------------------------------------------------------------------------------

void Structopt_sd::getgradduinternalforce(int inode,int typ,double& time,int anaId)
{

  int probType = getAnalysisAdjData(anaId, time);
    
  switch (structdom->probType()) {

  case SolverInfo::Static:
  case SolverInfo::Dynamic:
    if(structdom->isComplex())
      {
	assert(0);
      }
    else
      {
	adj->zero();
	structdom->getGradduNodalInternalForce(*sol, *adj, bcx, inode, typ);
      }
    break;
    /*
  case SolverInfo::NonLinStatic:
  case SolverInfo::NonLinDynam:
    adj->zero();
    structdom->getGradduNodalInternalForceNL(*gs, aC, *adj, inode, typ);
    break;
    */
  default:    
    fprintf(stderr,"Error: derivatives of nodal internal forces can only be\n");
    fprintf(stderr,"       evaluated for linear (quasi) static problems.\n");
    exit(-1);
  }

}				       
  
//------------------------------------------------------------------------------

double Structopt_sd::getstressint(int typ, double& refStress, double& powFac, 
                               int areaFlag, int nodeleTyp, 
                               int* list, int listsize, double& time, int anaId)
{

  int surface,stressIndex;
  getStressInfo(typ,surface,stressIndex);  

  int probType = getAnalysisData(anaId,time);

  double retVal;

  if (! nodeleTyp ) {

    switch (probType) {
    
    case SolverInfo::Static:
    case SolverInfo::Dynamic:
      if(structdom->isComplex())
	{
	  assert(0);
	}
      else
	{
	  retVal=structdom->getStressIntegral(*sol, bcx, refStress, 
					      powFac,list,listsize, 
					      stressIndex, areaFlag);
	}
      break;
      /*
    case SolverInfo::NonLinStatic:
    case SolverInfo::NonLinDynam:
      retVal=structdom->getNLStressIntegral(gs, aC, refStress, 
					    powFac, list, listsize, 
					    stressIndex, areaFlag);
      break;  	   
      */
    default:
      fprintf(stderr,"Error: Stress-Integral can only be evaluated for\n");
      fprintf(stderr,"       static and dynamic problems.\n");
      exit(-1);		      
    }
  }
  else {

    double nodStress;
    double sumVal=0.0;
    double maxval=0.0;

    int numnodes = structdom->numNode();

    int inode;
    for (inode=0;inode<numnodes;inode++){

      switch (probType) {
      case SolverInfo::Static:
      case SolverInfo::Dynamic:
	if(structdom->isComplex())
	  {
	    assert(0);
	  }
	else
	  {
	    nodStress=structdom->getNodalStressStrain(*sol, bcx, inode,
						      stressIndex, surface);
	  }
	break;
	/*
      case SolverInfo::NonLinStatic:
      case SolverInfo::NonLinDynam:
	nodStress=structdom->getNLNodalStressStrain(gs, aC, inode, 
						    stressIndex,surface);
	break;
	*/
      default:
	fprintf(stderr,"Error: Stress-Integral can only be evaluated for\n");
	fprintf(stderr,"       static and dynamic problems.\n");
	exit(-1);	               
      }
      sumVal+= pow(nodStress/refStress,powFac);
      maxval=max(maxval,nodStress);
    }

    double quoVal = areaFlag ? numnodes : 1.0;

    retVal=pow(sumVal/quoVal,1.0/powFac);
  }

  return retVal;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradstressint(int typ, double& refStress, double& powFac, 
                                   int areaFlag, int nodeleTyp, int* list,
                                   int listsize,double& time,int anaId)
{

  int surface,stressIndex;
  getStressInfo(typ,surface,stressIndex);  

  int probType = getAnalysisGradData(anaId, time);
 
  double retVal;

  if (! nodeleTyp ) {

    switch (probType) {    
    case SolverInfo::Static:
    case SolverInfo::Dynamic:
      if(structdom->isComplex())
	{
	  assert(0);
	}
      else
	{
	  retVal=structdom->getGradStressIntegral(*sol, *grad, bcx, refStress, 
						  powFac,list,listsize,
						  stressIndex, areaFlag);
	}
      break;
      /*
    case SolverInfo::NonLinStatic:
    case SolverInfo::NonLinDynam:
      retVal=structdom->getNLGradStressIntegral(gs, aC, *grad, refStress, 
						powFac, list, listsize, 
						stressIndex, areaFlag);
      break;
      */
    default:
      fprintf(stderr,"Error: Stress-Integral derivative can only be evaluated \n");
      fprintf(stderr,"       for static and dynamic problems.\n");
      exit(-1);		      
    }
  }
  else {

    double sumVal=0.0;
    double gradVal=0.0;
    double nodStress, gradStress;

    int numnodes = structdom->numNode();

    int inode;
    for (inode=0;inode<numnodes;inode++){

      switch (probType) {
      
      case SolverInfo::Static:
      case SolverInfo::Dynamic:
	if(structdom->isComplex())
	  {
	    assert(0);
	  }
	else
	  {
	    nodStress  = structdom->getNodalStressStrain(*sol, bcx, inode,
							 stressIndex, surface);
	    
	    gradStress = structdom->getGradNodalStressStrain(*sol,*grad,bcx,
							     inode,stressIndex,
							     surface);
	  }
	break;
	/*
      case SolverInfo::NonLinStatic:
      case SolverInfo::NonLinDynam:
	nodStress  = structdom->getNLNodalStressStrain(gs, aC, inode, 
						       stressIndex, surface);
    	  						 
	gradStress = structdom->getNLGradNodalStressStrain(gs, aC, *grad, 
							   inode, stressIndex,surface);
	break;
	*/
      default:
	fprintf(stderr,"Error: Stress-Integral derivative can only be evaluated \n");
	fprintf(stderr,"	 for static and dynamic problems.\n");
	exit(-1);
      }	  			
      sumVal  += pow(nodStress/refStress,powFac);
      gradVal += gradStress * pow(nodStress/refStress,powFac-1) /refStress;
    }

    double quoVal = areaFlag ? numnodes : 1.0;

    retVal=pow(sumVal/quoVal,1.0/powFac-1.0)*gradVal/quoVal;
  
  }

  return retVal;

}

//------------------------------------------------------------------------------

void Structopt_sd::getgraddustressint(int typ, double& refStress, double& powFac, 
                                   int areaFlag, int nodeleTyp,int* list,
                                   int listsize,double& time,int anaId)
{    
  int surface,stressIndex;
  getStressInfo(typ,surface,stressIndex);  

  int probType = getAnalysisAdjData(anaId, time);
  
  adj->zero();

  if (! nodeleTyp ) {
 
    switch (probType) {
    
    case SolverInfo::Static:
    case SolverInfo::Dynamic:
      if(structdom->isComplex())
	{
	  assert(0);
	}
      else
	{
	  structdom->getGradduStressIntegral(*sol, *adj, bcx, refStress, powFac, 
					     list,listsize,stressIndex, areaFlag);
	}
      break;
      /*
    case SolverInfo::NonLinStatic:
    case SolverInfo::NonLinDynam:
      structdom->getNLGradduStressIntegral(gs, aC, *adj, refStress, powFac, list,
					   listsize, stressIndex, areaFlag);
      break;
      */
    default:
      fprintf(stderr,"Error: Stress-Integral derivative can only be evaluated \n");
      fprintf(stderr,"       for static and dynamic problems.\n");
      exit(-1);		      
    }
  }
  else {

    double sumVal=0.0;
    double nodStress;

    int numnodes = structdom->numNode();

    Vector tmpAdj(adj->size());

    int inode;
    for (inode=0;inode<numnodes;inode++){

      switch (probType) {
      
      case SolverInfo::Static:
      case SolverInfo::Dynamic:
	if(structdom->isComplex())
	  {
	    assert(0);
	  }
	else
	  {
	    nodStress = structdom->getNodalStressStrain(*sol, bcx, inode,
							stressIndex, surface);
	  }
	break;
	/*
      case SolverInfo::NonLinStatic:
      case SolverInfo::NonLinDynam:
	nodStress = structdom->getNLNodalStressStrain(gs,aC, inode, 
						      stressIndex, surface);
	*/
      default:
	fprintf(stderr,"Error: Stress-Integral derivative can only be evaluated \n");
	fprintf(stderr,"	 for static and dynamic problems.\n");
	exit(-1);
      }	  			
      sumVal  += pow(nodStress/refStress,powFac);
    }

    double quoVal = areaFlag ? numnodes : 1.0;

    sumVal = pow(sumVal/quoVal,1.0/powFac-1.0)/quoVal;

    for (inode=0;inode<numnodes;inode++){

      tmpAdj.zero();
      
      switch (probType) {
      
      case SolverInfo::Static:
      case SolverInfo::Dynamic:
	if(structdom->isComplex())
	  {
	    assert(0);
	  }
	else
	  {
	    nodStress  = structdom->getNodalStressStrain(*sol, bcx, inode,
							 stressIndex, surface);
	  }
	break;
	/*
      case SolverInfo::NonLinStatic:
      case SolverInfo::NonLinDynam:
	nodStress  = structdom->getNLNodalStressStrain(gs,aC, inode, 
						       stressIndex, surface);  
	*/
      default:
	fprintf(stderr,"Error: Stress-Integral derivative can only be evaluated \n");
	fprintf(stderr,"	 for static and dynamic problems.\n");
	exit(-1);
      }	  			

      double gradVal    = sumVal * pow(nodStress/refStress,powFac-1) /refStress;
      
      switch (probType) {
      
      case SolverInfo::Static:
      case SolverInfo::Dynamic:
	if(structdom->isComplex())
	  {
	    assert(0);
	  }
	else
	  {
	    structdom->getGradduNodalStressStrain(*sol, tmpAdj, bcx, inode, 
						  stressIndex, surface);
	  }
	break;
	/*
      case SolverInfo::NonLinStatic:
      case SolverInfo::NonLinDynam:
	structdom->getNLGradduNodalStressStrain(gs, tmpAdj, aC, inode,
						stressIndex, surface);
	break;
	*/
      default:
	fprintf(stderr,"Error: Stress-Integral derivative can only be evaluated \n");
	fprintf(stderr,"	 for static and dynamic problems.\n");
	exit(-1);
      }	  			
      (*adj) +=  gradVal * tmpAdj;
    }
  }

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradpartstressint(int typ,double& refStress,double& powFac,
                                       int areaFlag,int nodeleTyp,int* list,
                                       int listsize,double& time,int anaId)
{

  int surface,stressIndex;
  getStressInfo(typ,surface,stressIndex);  

  int probType = getAnalysisData(anaId,time);
 
  double retVal;

  if (! nodeleTyp ) {

    switch (probType) {
    
    case SolverInfo::Static:
    case SolverInfo::Dynamic:
      if(structdom->isComplex())
	{
	  assert(0);
	}
      else
	{
	  retVal=structdom->getGradPartStressIntegral(*sol, *grad, bcx, refStress,
						      powFac, list, listsize, 
						      stressIndex, areaFlag);
	}
      break;
      /*
    case SolverInfo::NonLinStatic:
    case SolverInfo::NonLinDynam:
      retVal=structdom->getNLGradPartStressIntegral(gs, aC,*grad, refStress, 
						    powFac, list, listsize,
						    stressIndex, areaFlag);
      break;
      */
    default:
      fprintf(stderr,"Error: Stress-Integral derivative can only be evaluated \n");
      fprintf(stderr,"       for static and dynamic problems.\n");
      exit(-1);
    }
  }
  else {

    double sumVal=0.0;
    double gradVal=0.0;
    double nodStress, gradStress;

    int numnodes = structdom->numNode();

    int inode;
    for (inode=0;inode<numnodes;inode++){

      switch (probType) {
      
      case SolverInfo::Static:
      case SolverInfo::Dynamic:
	if(structdom->isComplex())
	  {
	    assert(0);
	  }
	else
	  {
	    nodStress  = structdom->getNodalStressStrain(*sol, bcx, inode,
							 stressIndex, surface);
	    gradStress = structdom->getGradPartNodalStressStrain(*sol, bcx, inode,
								 stressIndex,
								 surface);
	  }
	break;
	/*
      case SolverInfo::NonLinStatic:
      case SolverInfo::NonLinDynam:
	nodStress  = structdom->getNLNodalStressStrain(gs,aC, inode, 
						       stressIndex, surface);
	gradStress = structdom->getNLGradPartNodalStressStrain(gs, aC, inode,
							       stressIndex, surface);  					
	*/
      default:
	fprintf(stderr,"Error: Stress-Integral derivative can only be evaluated \n");
	fprintf(stderr,"	 for static and dynamic problems.\n");
	exit(-1);
      }	  			
      sumVal  += pow(nodStress/refStress,powFac);
      gradVal += gradStress * pow(nodStress/refStress,powFac-1) /refStress;
    }

    double quoVal = areaFlag ? numnodes : 1.0;

    retVal=pow(sumVal/quoVal,1.0/powFac-1.0)*gradVal/quoVal;
  
  }

  return retVal;

}

//------------------------------------------------------------------------------

double Structopt_sd::getdisplevel(int aFlag, int quoFlag, int size, 
				  int * nodeid, int*distyp, double* refVal, 
				  double* difVal, double& powFac, double& time, 
				  int anaId)
{
  double val;
  int probType = getAnalysisGradData(anaId,time);
  switch (probType) 
    {
    case SolverInfo::Static:
      if(structdom->isComplex())
	{
	  val=structdom->getDisplevel(*csol, cbcx, aFlag, quoFlag, size,
				      nodeid, distyp, refVal, difVal, powFac);
	}
      else
	{
	  val=structdom->getDisplevel(*sol, bcx, aFlag, quoFlag, size,
				      nodeid, distyp, refVal, difVal, powFac);
	}
      break;
    default:
      assert(0);
    }
  return val;
}				       
				       
//------------------------------------------------------------------------------

double Structopt_sd::getgraddisplevel(int aFlag, int quoFlag, int size, 
                                   int * nodeid, int*distyp, double* refVal, 
                                   double* difVal, double& powFac, double& time,
				   int anaId)
{
  double val;
  int probType = getAnalysisGradData(anaId,time);
  int apb = analysisMap[anaId];   
  switch (probType) 
    {
    case SolverInfo::Static:
      if(structdom->isComplex())
	{
	  val = structdom->getGradDisplevel(*csol,  *cgrad, cbcx, aFlag, quoFlag, size, nodeid, distyp, 
					    refVal, difVal, powFac);
	}
      else
	{
	  val = structdom->getGradDisplevel(*sol,  *grad, bcx, aFlag, quoFlag, size, nodeid, distyp, 
					    refVal, difVal, powFac);
	}
      break;
    default:
      assert(0);
    }
  return val;
}				       


//------------------------------------------------------------------------------
void Structopt_sd::getgraddudisplevel(int aFlag, int quoFlag, int size, 
				      int * nodeid, int*distyp, double* refVal, 
				      double* difVal, double& powFac, double& time,
				      int anaId )

{
  int probType = getAnalysisAdjData(anaId,time);
  int apb = analysisMap[anaId];   
  switch (probType) 
    {
    case SolverInfo::Static:
      if(structdom->isComplex())
	{
	  cadj->zero();
	  structdom->getGradduDisplevel(*csol,  *cadj, cbcx, aFlag, quoFlag, size, nodeid, distyp, 
					refVal, difVal, powFac);
	}
      else
	{
	  adj->zero();
	  structdom->getGradduDisplevel(*sol,  *adj, bcx, aFlag, quoFlag, size, nodeid, distyp, 
					refVal, difVal, powFac);
	}
      break;
    default:
      assert(0);
    }
  return;
} 

//------------------------------------------------------------------------------

double Structopt_sd::getstresslevel(int aFlag, int quoFlag, int eFlag, int size, 
                                 int * nodeid, int*stresstyp, double* refVal, 
                                 double* difVal, double& powFac, double& time,
				 int anaId)
{

  // check input data

  int numele = structdom->numElements();
  int numnod = structdom->numNodes();
  
  int ierror=0;

  int i;
  for (i=0;i<size;i++) {
    if (!eFlag && nodeid[i] >= numele) ierror=1;
    if ( eFlag && nodeid[i] >= numnod) ierror=1;
  }
  
  if (ierror) {
    fprintf(stderr,"Error: node/element id > numnodes/numelem\n");
    fprintf(stderr,"	  in stress leveling........... stop\n");
    exit(-1);
  }

  // evaluate stress level criterion

  int probType = getAnalysisData(anaId, time);

  int node,typ;

  int iele=-1;

  double stressVal,locDif;

  double aveVal = 0.0;
  double sumVal = 0.0;

  double * tmpVal = static_cast<double*>(dbg_alloca(sizeof(double)*size));

  for (i=0;i<size;i++) 
    {   
      node   = nodeid[i];
      typ    = stresstyp[i];

      if (!eFlag) iele = node;

      int surface,stressIndex;
      getStressInfo(typ,surface,stressIndex);  

      switch (structdom->probType()) {
      case SolverInfo::Static:
      case SolverInfo::Dynamic:
	if(structdom->isComplex())
	  {
	    assert(0);
	  }
	else
	  {
	    tmpVal[i] = structdom->getNodalStressStrain(*sol, bcx, node,
							stressIndex, surface,iele);
	  }
        break;
	/*
      case SolverInfo::NonLinStatic:
      case SolverInfo::NonLinDynam:
        tmpVal[i] = structdom->getNLNodalStressStrain(gs, aC, node,
          					      stressIndex, surface);
        break;
      case SolverInfo::FVibr:
 	tmpVal[i] = structdom->getNodalStressStrain(*csol, cbcx, node,
						    stressIndex, surface,iele);
        break;
	*/
      default:
        fprintf(stderr,"ERROR: wrong analysis type in stress leveling\n");
	exit(-1);
      }
    }

  if (aFlag) 
    {
      for (i=0;i<size;i++) aveVal += tmpVal[i];
      aveVal = aveVal / size;
    } 
  
  for (i=0;i<size;i++)
    {
      if (aFlag) 
	locDif = aveVal;
      else
	locDif = difVal[i];

      stressVal = (tmpVal[i]-locDif)/refVal[i];

      sumVal+= pow(stressVal,powFac);
    }  

  double quoVal = quoFlag ? size : 1.0;

  return pow(sumVal/quoVal,1.0/powFac);

}				       
				       
//------------------------------------------------------------------------------

double Structopt_sd::getgradstresslevel(int aFlag,int quoFlag,int eFlag,int size, 
                                     int *nodeid,int*stresstyp,double* refVal, 
                                     double* difVal,double& powFac,double& time,
				     int anaId)
{
  int probType = getAnalysisGradData(anaId, time);

  int i,node,typ;

  int iele=-1;

  double stressVal,gradstressVal,locDif,gradlocDif;

  double aveVal     = 0.0;
  double gradaveVal = 0.0;
  double sumVal     = 0.0;
  double gradsumVal = 0.0;

  double * tmpVal = static_cast<double*>(dbg_alloca(sizeof(double)*size));
  double * tmpGrd = static_cast<double*>(dbg_alloca(sizeof(double)*size));

  for (i=0;i<size;i++) 
    {   
      node	= nodeid[i];
      typ	= stresstyp[i];
 
      if (!eFlag) iele = node;

      int surface,stressIndex;
      getStressInfo(typ,surface,stressIndex);  

      switch (structdom->probType()) {
      case SolverInfo::Static:
      case SolverInfo::Dynamic: 
	if(structdom->isComplex())
	  {
	    assert(0);
	  }
	else
	  {
	    tmpVal[i] = structdom->getNodalStressStrain(*sol, bcx, node,
							stressIndex, surface,iele);
	    tmpGrd[i] = structdom->getGradNodalStressStrain(*sol, *grad, bcx, node,
							    stressIndex,surface,iele);
	  }
        break;
	/*
      case SolverInfo::NonLinStatic:
      case SolverInfo::NonLinDynam:
        tmpVal[i] = structdom->getNLNodalStressStrain(gs, aC, node,
         	    				      stressIndex, surface);
        tmpGrd[i] = structdom->getNLGradNodalStressStrain(gs,aC,*grad,node,
                                                          stressIndex,surface);
        break;
      case SolverInfo::FVibr:
	tmpVal[i] = structdom->getNodalStressStrain(*csol, cbcx, node,
						    stressIndex, surface,iele);
	tmpGrd[i] = structdom->getGradNodalStressStrain(*csol, *cgrad, cbcx, node,
	 					        stressIndex,surface,iele);
        break;
	*/
      default:
        fprintf(stderr,"ERROR: wrong analysis type in stress leveling\n");
	exit(-1);
      }
    }

  if (aFlag) 
    {
      for (i=0;i<size;i++) {
	aveVal     += tmpVal[i];
	gradaveVal += tmpGrd[i];
      }
      aveVal     = aveVal / size;
      gradaveVal = gradaveVal / size;
    } 

  for (i=0;i<size;i++)
    {
      if (aFlag) {
	locDif     = aveVal;
	gradlocDif = gradaveVal;
      }
      else {
	locDif     = difVal[i];
	gradlocDif = 0.0;
      }

      gradstressVal = (tmpGrd[i] - gradlocDif)/refVal[i];

      stressVal     = (tmpVal[i] - locDif)/refVal[i];

      sumVal    += pow(stressVal,powFac);
    
      gradsumVal+= gradstressVal * pow(stressVal,powFac-1.0);

    }  

  double quoVal = quoFlag ? size : 1.0;

  return  pow(sumVal/quoVal,1.0/powFac-1.0)*gradsumVal/quoVal;

}				       
				       
//------------------------------------------------------------------------------

double Structopt_sd::getgradpartstresslevel(int aFlag,int quoFlag,int eFlag,
                                         int size,int * nodeid,int*stresstyp,
                                         double* refVal,double* difVal, 
                                         double& powFac, double& time,int anaId)
{

  int probType = getAnalysisData(anaId, time);

  int i,node,typ;

  int iele=-1;

  double stressVal,locDif,gradstressVal,gradlocDif;

  double aveVal     = 0.0;
  double gradaveVal = 0.0;
  double sumVal     = 0.0;
  double gradsumVal = 0.0;

  double * tmpVal = static_cast<double*>(dbg_alloca(sizeof(double)*size));
  double * tmpGrd = static_cast<double*>(dbg_alloca(sizeof(double)*size));

  for (i=0;i<size;i++) 
    {   
      node	= nodeid[i];
      typ	= stresstyp[i];

      if (!eFlag) iele = node;

      int surface,stressIndex;
      getStressInfo(typ,surface,stressIndex);  

      switch (structdom->probType()) {
      case SolverInfo::Static:
      case SolverInfo::Dynamic:
	if(structdom->isComplex())
	  {
	    assert(0);
	  }
	else
	  {
	    tmpVal[i] = structdom->getNodalStressStrain(*sol, bcx, node,
							stressIndex, surface,iele);
	    tmpGrd[i] = structdom->getGradPartNodalStressStrain(*sol, bcx, node,
								stressIndex,surface,iele);
	  }
        break;
	/*
      case SolverInfo::NonLinStatic:
      case SolverInfo::NonLinDynam:
        tmpVal[i] = structdom->getNLNodalStressStrain(gs, aC, node,
  		        			      stressIndex, surface);
        tmpGrd[i] = structdom->getNLGradPartNodalStressStrain(gs, aC, node,
							      stressIndex,surface);
        break;
      case SolverInfo::FVibr:
	tmpVal[i] = structdom->getNodalStressStrain(*csol, cbcx, node,
						    stressIndex, surface,iele);
	tmpGrd[i] = structdom->getGradPartNodalStressStrain(*csol, cbcx, node,
							    stressIndex,surface,iele);
        break;
	*/
      default:
        fprintf(stderr,"ERROR: wrong analysis type in stress leveling\n");
	exit(-1);
      }
    }


  if (aFlag) 
    {
      for (i=0;i<size;i++) {
	aveVal     += tmpVal[i];
	gradaveVal += tmpGrd[i];
      }
      aveVal     = aveVal / size;
      gradaveVal = gradaveVal / size;
    } 

  for (i=0;i<size;i++)
    {
      if (aFlag) {
	locDif     = aveVal;
	gradlocDif = gradaveVal;
      }
      else {
	locDif     = difVal[i];
	gradlocDif = 0.0;
      }

      gradstressVal = (tmpGrd[i] - gradlocDif)/refVal[i];

      stressVal     = (tmpVal[i] - locDif)/refVal[i];

      sumVal       += pow(stressVal,powFac);
    
      gradsumVal   += gradstressVal * pow(stressVal,powFac-1.0);

    }  

  double quoVal = quoFlag ? size : 1.0;

  return  pow(sumVal/quoVal,1.0/powFac-1.0)*gradsumVal/quoVal;

}				       

//------------------------------------------------------------------------------

void Structopt_sd::getgraddustresslevel(int aFlag,int quoFlag,int eFlag,int size, 
                                     int * nodeid,int*stresstyp,double* refVal, 
                                     double* difVal,double& powFac,double& time,
				     int anaId)

{

  int probType = getAnalysisAdjData(anaId, time, 1);
  
  int i,j,node,typ;

  int iele=-1;

  double stressVal,locDif;

  double aveVal     = 0.0;
  double sumVal     = 0.0;

  double * tmpVal = static_cast<double*>(dbg_alloca(sizeof(double)*size));

  Vector*        tmpVec  = 0;  
  ComplexVector* ctmpVec = 0;

  if (adj)   tmpVec = new Vector(*adj);
  if (cadj) ctmpVec = new ComplexVector(*cadj);
  
  for (i=0;i<size;i++) 
    {   
      node   = nodeid[i];
      typ    = stresstyp[i];

      if (!eFlag) iele = node;

      int surface,stressIndex;
      getStressInfo(typ,surface,stressIndex);  

      switch (structdom->probType()) {
      case SolverInfo::Static:
      case SolverInfo::Dynamic: 
	if(structdom->isComplex())
	  {
	    assert(0);
	  }
	else
	  {
	    tmpVal[i] = structdom->getNodalStressStrain(*sol, bcx, node,
							stressIndex,surface,iele);
	  }
        break;
	/*
      case SolverInfo::NonLinStatic:
      case SolverInfo::NonLinDynam:
        tmpVal[i] = structdom->getNLNodalStressStrain(gs, aC, node,
  						      stressIndex,surface);
        break;
      case SolverInfo::FVibr:
	tmpVal[i] = structdom->getNodalStressStrain(*csol, cbcx, node,
	        				    stressIndex,surface,iele);
        break;
	*/
      default:
        fprintf(stderr,"ERROR: wrong analysis type in stress leveling\n");
	exit(-1);
      }
    }

  if (aFlag) 
    {
      for (i=0;i<size;i++) aveVal += tmpVal[i];
      aveVal     = aveVal / size;
    } 

  for (i=0;i<size;i++)
    {
      if (aFlag) 
	locDif = aveVal;
      else
	locDif = difVal[i];
 
      stressVal = (tmpVal[i] - locDif)/refVal[i];

      sumVal   += pow(stressVal,powFac);
    }  

  double quoVal = quoFlag ? size : 1.0;

  double allFac = pow(sumVal/quoVal,1.0/powFac-1.0) / quoVal; 

  for (i=0;i<size;i++)
    {
      if (aFlag) 
	locDif = aveVal;
      else
	locDif = difVal[i];

      stressVal = (tmpVal[i]-locDif)/refVal[i];
      sumVal    = pow(stressVal,powFac-1.0)/refVal[i];

      node   = nodeid[i];
      typ    = stresstyp[i];

      if (!eFlag) iele = node;

      int surface,stressIndex;
      getStressInfo(typ,surface,stressIndex);  

      switch (structdom->probType()) {
      case SolverInfo::Static:
      case SolverInfo::Dynamic: 
	if(structdom->isComplex())
	  {
	    assert(0);
	  }
	else
	  {
	    tmpVec->zero();
	    structdom->getGradduNodalStressStrain(*sol, *tmpVec, bcx,
						  node, stressIndex, surface,iele);
	    (*adj)  += allFac * sumVal * (*tmpVec);
	  }
        break;
	/*
      case SolverInfo::NonLinStatic:
      case SolverInfo::NonLinDynam:
        tmpVec->zero();
        structdom->getNLGradduNodalStressStrain(gs, *tmpVec, aC, node,
	                                        stressIndex, surface);
        (*adj)  += allFac * sumVal * (*tmpVec);
        break;
      case SolverInfo::FVibr:
        ctmpVec->zero();
        structdom->getGradduNodalStressStrain(*csol, *ctmpVec, cbcx,
					      node, stressIndex, surface,iele);

	(*cadj)  += allFac * sumVal * (*ctmpVec);
        break;
	*/
      default:
        fprintf(stderr,"ERROR: wrong analysis type in stress leveling\n");
	exit(-1);
      }

      if (aFlag) {

	for (j=0;j<size;j++) 
	  {

	    stressVal = (tmpVal[j]-locDif)/refVal[j];

	    sumVal = pow(stressVal,powFac-1.0) / refVal[j] / quoVal;

	    switch (structdom->probType()) {
	    case SolverInfo::Static:
	    case SolverInfo::Dynamic: 
	    case SolverInfo::NonLinStatic:
	    case SolverInfo::NonLinDynam:
	      (*adj) -= allFac * sumVal * (*tmpVec);
	      break;
	      /*
	    case SolverInfo::FVibr:
	      (*cadj)  -= allFac * sumVal * (*ctmpVec);
	      */
	    default:
	      fprintf(stderr,"ERROR: wrong analysis type in stress leveling\n");
	      exit(-1);
	    }
	  }
      }
    }  

  if (tmpVec)  delete tmpVec;
  if (ctmpVec) delete ctmpVec;
} 

//------------------------------------------------------------------------------

double Structopt_sd::getFailprob(int typ, int icrit)
{

  double val;
  
  switch(typ) {
  case 0:      // failure probability
    val=reliabilityProb->relsol->getFailprob(icrit);
    break;
  case 1:      // reliability index
    val=reliabilityProb->relsol->getRelindex(icrit);
    break;
  case 2:      // reliability index
    val=reliabilityProb->relsol->getPMAvalue(icrit);
    break;
  default:
    fprintf(stderr,
	    "Error: Reliability criterion not implemented\n");
  }
        
  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradFailprob(int typ, int icrit, int ivar)
{
  double val;

  // get id of abstract variable 
  // special treatment needed for adjoint approach

  int actv = ivar;

  if (ivar < 0) actv= optpro->optvar->getIOldVar();

  // get gradients stored in reliability solver class
  
  switch(typ) {
  case 0:      // failure probability
    val=reliabilityProb->relsol->getgradFailprob(icrit,actv);
    break;
  case 1:      // reliability index
    val=reliabilityProb->relsol->getgradRelindex(icrit,actv);
    break;
  case 2:      // reliability index
    val=reliabilityProb->relsol->getgradPMAvalue(icrit,actv);
    break;
  default:
    fprintf(stderr,
	    "Error: Reliability criterion not implemented\n");
  }

  return val;
}

//------------------------------------------------------------------------------

double Structopt_sd::getPullIn(int anaId) 
{
  
  double val;  

#ifdef AEROELASTIC

  // send opt parameters so that electrostatic knows what to do
  sndOptpar(-1,-1);

  double time = 0.0;

  int probType = getAnalysisData(anaId, time);

  int apb = analysisMap[anaId];

  switch (probType) {
    /*
  case SolverInfo::Dynamic:
    val    = dynsolvers[apb]->strcEmEigSolve(1);
    break;
    */
  default:
    fprintf(stderr,"Error: pull-in can only be evaluated for\n");
    fprintf(stderr,"       linear electro-mechanical problems.\n");
    exit(-1);
  }

#else

  fprintf(stderr,"Error: pull-in can only be evaluated for\n");
  fprintf(stderr,"       linear electro-mechanical problems...... stop\n");
  exit(-1);

#endif

  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradPullIn(int actvar,int anaId) 
{
  
  double val;  

#ifdef AEROELASTIC

  // send opt parameters so that electrostatic knows what to do
  sndOptpar(actvar,-1);

  double time = 0.0;

  int probType = getAnalysisData(anaId, time);

  int apb = analysisMap[anaId];

  switch (probType) {
    
  case SolverInfo::Dynamic:
    val    = dynsolvers[apb]->gradStrcEmEigSolve();
    break;
  default:
    fprintf(stderr,"Error: gradients of pull-in can only be evaluated for\n");
    fprintf(stderr,"       linear electro-mechanical problems.\n");
    exit(-1);
  }
#else

  fprintf(stderr,"Error: gradients of pull-in can only be evaluated for\n");
  fprintf(stderr,"       linear electro-mechanical problems...... stop\n");
  exit(-1);

#endif

  return val;

}

//------------------------------------------------------------------------------

double Structopt_sd::getgradpartPullIn(int anaId) 
{
  
  double val;  

#ifdef AEROELASTIC

  double time = 0.0;

  int probType = getAnalysisData(anaId, time);

  int apb = analysisMap[anaId];

  switch (probType) {

  case SolverInfo::Dynamic:
    val    = dynsolvers[apb]->gradStrcEmEigSolve();
    break;
  default:
    fprintf(stderr,"Error: gradients of pull-in can only be evaluated for\n");
    fprintf(stderr,"       linear electro-mechanical problems...... stop\n");
    exit(-1);
  }
#else

  fprintf(stderr,"Error: gradients of pull-in can only be evaluated for\n");
  fprintf(stderr,"	  linear electro-mechanical problems...... stop\n");
  exit(-1);

#endif

  return val;

}

//------------------------------------------------------------------------------

void Structopt_sd::graddirect(int ivar, double** gc, int numAC, int* listAC) 
{

  anagrdType = 1;                       //Set type of SA

  // evaluate gradients of deterministic criteria 

  if ( ! failcritOnly ) {

#ifdef AEROELASTIC
    sndOptpar(ivar,-1);
#endif

    int ifail = 1;
 
    int ina;
    for (ina=0;ina<numAnalysis;ina++) {
 
      structdom->activateAnalysis(ina);
 
      int apb = analysisMap[ina];
 
      switch (structdom->probType()) {
      case SolverInfo::Static:
	if(structdom->isComplex())
	  {
	    c_staticsolvers[apb]->optSolveDuds();
	    ifail = 0;
	  }
	else
	  {
	    staticsolvers[apb]->optSolveDuds();
	    ifail = 0;
	  }
	break;
	/*
      case SolverInfo::NonLinStatic:
	if (ivar == 0)  nonlinsolvers[apb]->optStablePreSolve();
	nonlinsolvers[apb]->optSolveDuds();            
	ifail=0;
	break;
      case SolverInfo::Modal:
	eigensolvers[apb]->optSolveDwds();
	ifail = 0;
	break;
      case SolverInfo::Dynamic:
	dynsolvers[apb]->optDstateds(this,anagrdType);
	ifail = 0;
	break;
      case SolverInfo::NonLinDynam:
	//need the stabilty presolve here
	NLdynsolvers[apb]->optDstateds(this,anagrdType);
	ifail = 0;
	break;
      case SolverInfo::FVibr:
	fvibrsolvers[apb]->optSolveDuds();
	ifail = 0;
	break;
	*/
      default:
	assert(0);
      }
    }

    if (ifail) {
      fprintf(stderr,"Error: direct analytical sensitivity analysis is implemented\n");
      fprintf(stderr,"       only for static, eigen and dynamic single domain problems\n");
      exit(-1);
    }

#ifdef AEROELASTIC
    structdom->savegradExtDomainRes();
#endif
  }
   
  gradcrit(ivar,gc,numAC,listAC);

}

//------------------------------------------------------------------------------
void Structopt_sd::gradadjoint(int icrit, double** gc) 
{

  anagrdType = 2;                       //Set type of SA

  int probType = optpro->opc[icrit]->getAnalysisType();

  switch ( probType ) {

  case Optcrit::StressAnalysis:
    gradadjointStress(icrit, gc);
    break;

  case Optcrit::ModalAnalysis:
    gradadjointModal(icrit, gc);
    break;

  case Optcrit::GeometricAnalysis:
    gradadjointGeom(icrit, gc);
    break;
    
  case Optcrit::ExtCritAnalysis:
    gradadjointStress(icrit, gc);
    break; 

  default:
    fprintf(stderr,"Error: adjoint analytical sensitivity analysis \n");
    fprintf(stderr,"       Analysis Type is not implemented        \n");
    exit(-1);

  }
}

//------------------------------------------------------------------------------
void Structopt_sd::gradadjointModal(int icrit, double** gc) 
{

  // get analysis Id

  int anaId = optpro->opc[icrit]->getAnalysisId();

  int apb   = analysisMap[anaId];

  // check whether modal problem is defined

  structdom->activateAnalysis(anaId);  

  if (structdom->probType() != SolverInfo::Modal) {
    fprintf(stderr,"Error: adjoint sensitivity analysis for modal criteria require\n"); 
    fprintf(stderr,"       modal analysis  ... stop\n");
    exit(-1);
  }
  /*
  // get eigenvalue id of current criterion

  int ieig = optpro->opc[icrit]->getEigvalId();

  // get derivatives of eigenvalue wrt. all optimization variables

  Vector * eigval;
  Vector * eigrad;

  int i;
  int numabs = optpro->optvar->numabs; 

  for (i=0;i<numabs;i++) {                   //Loop over all opt.var.

    optpro->optvar->updgrad(i,1.0);         //Set Variation of Abs.&Struct. Var

    structdom->setoptInf(optInf[i]);        //Set Opt.Var.- Ele. Influence

    eigensolvers[apb]->optSolveDwds();

    eigval = eigensolvers[apb]->getpeigval();
    eigrad = eigensolvers[apb]->getpeigrad();

    gc[icrit][i] = (*eigrad)[ieig] / ( 4.0*3.14159265*sqrt( (*eigval)[ieig] ));
  } 
  */
}

//------------------------------------------------------------------------------

void Structopt_sd::gradadjointStress(int icrit, double** gc)
{

  // need the analysis ID in sndOptpar, to send to other domains

  int anaId = optpro->opc[icrit]->getAnalysisId();

  activeAnalysis = anaId;

#ifdef AEROELASTIC
  sndOptpar(-1,icrit);  
#endif

  // Set design velocity = 0

  zerograd();

  // Build dq:du

  optpro->opc[icrit]->graddu(this);

  //Solve adjoint problem: K(-1) * dq:du

  int apb   = analysisMap[anaId];

  int numabs = optpro->optvar->nAbs(); 

  double time = optpro->opc[icrit]->getTime();  //Get load case of criterion
  int lc = (time < 0) ? 0 : static_cast<int>(time);
 
  switch (structdom->probType()) {	 
  case SolverInfo::Static:  	  
    if(structdom->isComplex())
      {
	c_staticsolvers[apb]->optAdjSolve();
      }
    else
      {
	staticsolvers[apb]->optAdjSolve();
      }
    //grad = staticsolvers[apb]->getpgrad();
    break;
    /*
  case SolverInfo::NonLinStatic:
    if (icrit == 0)  nonlinsolvers[apb]->optStablePreSolve();    
    nonlinsolvers[apb]->optAdjSolve(dqdlam);
    //grad = nonlinsolvers[apb]->getGradSol();    
    break;
  case SolverInfo::Dynamic: 	      //Determine K(-1) * dq:du
    //need the preprocess Solve Here
    dynsolvers[apb]->optDstateds(this,anagrdType);
    //grad  = dynsolvers[apb]->getpgDis();
    break;
  case SolverInfo::NonLinDynam:
    NLdynsolvers[apb]->optDstateds(this,anagrdType,dqdlam);
    //grad = NLdynsolvers[apb]->getpgDis();
    break;
  case SolverInfo::FVibr:
    fvibrsolvers[apb]->optAdjSolve(lc);
    break;
    */
  default:
    fprintf(stderr,"Error: adjoint analytical sensitivity analysis is \n");
    fprintf(stderr,"       implemented only for static and dynamic problems.\n");
  }

  // Post-multiply with pseudo-load vector
   
  for (int i=0;i<numabs;i++) {              //Loop over all opt.var.

    //fprintf(stderr,"\r ... Post-multiplying optimization variable %d",i+1);

    optpro->optvar->updgrad(i,1.0);        //Set Variation of Abs.&Struct. Var

    structdom->setoptInf(optInf[i]);       //Set Opt.Var.- Ele. Influence
     
    optpro->opc[icrit]->gradpart(this);    //Determine partial(q)/partial(s)
       
    double grad1 = optpro->opc[icrit]->grad;

    int stcvarTyp = 0;
    
#ifdef AEROELASTIC

    stcvarTyp = structdom->getOptvarTyp();

    switch (stcvarTyp) {
    case Optvar::attribute:
      sndOptpar(i,-1);
      break;
    case Optvar::composite:
      sndOptpar(i,-1);
      break;
    case Optvar::elecattr:
      sndOptpar(i,-2);
      break;
    case Optvar::thermattr:
      sndOptpar(i,-3);
      break;
    case Optvar::mpattr:
      sndOptpar(i,-4);
      break;
    default:
      sndOptpar(i,icrit);  
    }

#endif

    double grad2 = 0; 
    //Determine (adjoint(T) * dK:ds * u)
    
    switch (structdom->probType()) 
      {   
      case SolverInfo::Static:
	if(structdom->isComplex())
	  {
	    grad2 = c_staticsolvers[apb]->optGetAdjointPseudo(lc);
	  }
	else
	  {
	    grad2 = staticsolvers[apb]->optGetAdjointPseudo(lc);	    
	  }
	break;
	/*
      case SolverInfo::NonLinStatic:
	grad2 = nonlinsolvers[apb]->optGetAdjointPseudo();
	break;
      case SolverInfo::Dynamic:
	grad2 = dynsolvers[apb]->optGetAdjointPseudo(stcvarTyp);
	break;
      case SolverInfo::NonLinDynam:
	grad2 = NLdynsolvers[apb]->optGetAdjointPseudo(stcvarTyp);
	break;
      case SolverInfo::FVibr:
	grad2 = fvibrsolvers[apb]->optGetAdjointPseudo(lc).real();
	break;
	*/
      default:
	fprintf(stderr,"Error: adjoint analytical sensitivity analysis is \n");
	fprintf(stderr,"       implemented only for static and dynamic problems.\n");	 
	exit(-1);
      }

#ifdef AEROELASTIC
    if (stcvarTyp != Optvar::attribute && stcvarTyp != Optvar::composite) {

      // so far the gradients of the criteria are always saved on dir=0
      structdom->savegradExtDomainRes(1);  

      // gradients of shape variations only come from the fluid and meshmotion

      grad2 += structdom->getgradAeroForce(0);
      grad2 += structdom->getgradMeshRes(0);

      // are these checks on the variable types necessary

      if (stcvarTyp == Optvar::elecattr || stcvarTyp == Optvar::mpattr) {
        grad2 += structdom->getgradElecRes(0);
      }
      if (stcvarTyp == Optvar::thermattr || stcvarTyp == Optvar::mpattr) {
        grad2 += structdom->getgradThermRes(0);
      }

    } 
#endif

    gc[icrit][i] = grad1 + grad2;      //Store dq:ds
  }

}

//------------------------------------------------------------------------------
void Structopt_sd::gradadjointGeom(int icrit, double** gc) 
{
  zerograd();                                  //Set design velocity = 0

  int i;
  int numabs = optpro->optvar->nAbs(); 

  double time = optpro->opc[icrit]->getTime();      //Get load case of criterion
  int lc = (time < 0) ? 0 : static_cast<int>(time);
  
  for (i=0;i<numabs;i++) {                          //Loop over all opt.var.

    optpro->optvar->updgrad(i,1.0);                //Set Variation of Abs.&Struct. Var

    structdom->setoptInf(optInf[i]);               //Set Opt.Var.- Ele. Influence
     
    optpro->opc[icrit]->gradpart(this);            //Determine partial(q)/partial(s)
       
    gc[icrit][i] = optpro->opc[icrit]->grad;       //Store dq:ds
  }
}

//------------------------------------------------------------------------------

void Structopt_sd::zerograd() 
{
  structdom->zeroGrad();
}

//------------------------------------------------------------------------------

int Structopt_sd::initProcess( double & stateTime, double & dt)
{
  maxstep = 0;
 
  return static_cast<int>(stateTime/dt);

}

//------------------------------------------------------------------------------

void Structopt_sd::procEvaluate(double time, double step ) {

  double tmax = time + 0.5 * step;
  double tmin = time - 0.5 * step;
  
  int i;
  for (i=0;i<optpro->numcrit;i++) { 
    optpro->opc[i]->evaluate(this, tmax, tmin);
  }
}

//------------------------------------------------------------------------------

void Structopt_sd::procSavevar()
{
  int numabs = optpro->optvar->nAbs();

  int iabs;
  for (iabs=0;iabs<numabs;iabs++) 
    optpro->optvar->getAbsVar(iabs)->saveval();

  maxstep=0;
}

//------------------------------------------------------------------------------

int Structopt_sd::procComputeMaxSteps(double& stepid)
{
  int mxstp = 0;

  int numabs  = optpro->optvar->nAbs();

  int minstep = analysisData->transMin;

  double limitdsgvel = analysisData->designVel;

  double machAcc  = 1.0e-12;

  double locid = max(1.0,stepid);

  int iabs;

  for (iabs=0;iabs<numabs;iabs++) {  

    optpro->optvar->getAbsVar(iabs)->setnewval();
  
    double newvar   = optpro->optvar->getAbsVar(iabs)->getnewval();
    double actvar   = optpro->optvar->getAbsVar(iabs)->getval();
    double delvar   = abs(newvar - actvar);

    if ( delvar < machAcc ) continue;

    double maxdsgvel = optpro->optvar->maxDesignVelocity(iabs);

    // design change by non-shape variation

    double one = 1.0;
    maxdsgvel = max(one,maxdsgvel);

    int maxdsg = static_cast<int>(( delvar * maxdsgvel ) / ( limitdsgvel * locid ));

    mxstp=max(mxstp,maxdsg);
    mxstp=max(mxstp,minstep);
  }

  return mxstp;
}

//------------------------------------------------------------------------------

void Structopt_sd::procComputeStepSize(int mxstp)
{
  double twopi    = 6.2831853;
  double pihalf   = 1.5707963;
  double machAcc  = 1.0e-12;

  int numabs = optpro->optvar->nAbs();

  int iabs;

  for (iabs=0;iabs<numabs;iabs++) {  
  
    double newvar   = optpro->optvar->getAbsVar(iabs)->getnewval();
    double actvar   = optpro->optvar->getAbsVar(iabs)->getval();
    double delvar   = abs(newvar - actvar);

    if ( delvar < machAcc ) continue;
  
    double sum    = 0.0;

    int is;
    for (is=0;is<mxstp;is++) 
      sum += 1.0 + sin(twopi*static_cast<double>(is+1)/mxstp - pihalf);
      
    optpro->optvar->getAbsVar(iabs)->setdynscl((newvar - actvar) / sum);
  }
}

//------------------------------------------------------------------------------

int Structopt_sd::procTransition(int nstp, int mxstp)
{

  double twopi    = 6.2831853;
  double pihalf   = 1.5707963;
  double machAcc  = 1.0e-12;

  // initialize return value

  int retval=2;

  // compute intermediate values for variables

  int numabs = optpro->optvar->nAbs();
  int iabs;

  for (iabs=0;iabs<numabs;iabs++) {
  
    double newvar = optpro->optvar->getAbsVar(iabs)->getnewval();
    double actvar = optpro->optvar->getAbsVar(iabs)->getval();
    double dynscl = optpro->optvar->getAbsVar(iabs)->getdynscl();
    double delvar   = abs(newvar - actvar);

    if ( delvar < machAcc ) continue;
        
    retval = 1;

    actvar += dynscl *(1.0 + sin(twopi*static_cast<double>(nstp)/mxstp - pihalf));
    
    double reach = (newvar-actvar)*dynscl;
   
    if ( reach > 0.0 && nstp < mxstp) 
      retval=0;
    else  
      actvar=newvar;

    optpro->optvar->getAbsVar(iabs)->setval(actvar);
  }     
  
  // update structural variables

  optpro->optvar->updvar();

  return retval;
}


//------------------------------------------------------------------------------

int Structopt_sd::procSetvar(double time, double step ) {

  // return value:  2 - optimization variables not changed; stop transition
  //                1 - optimization variables changed; stop transition
  //                0 - optimization variables changed; continue transition

  if ( ! analysisData->transFlag ) return 2;

  if ( time == 0.0 && analysisData->initBlend == 0 ) return 2;
  
  int steplimit = analysisData->transMax;

  if ( steplimit <= 0 ) return 2;
  
  int numstep = static_cast<int>(time/step);

  // initialize and save old state for which equilibrium has been computed

  // numstep = 0  initialize 
  //              we assume that if reliability probelm is involved
  //              a transition process is always needed
  // numstep < 0  save old state

  //if (numstep == 0 ) {
  //  maxstep=0;
  //  if (!reliabilityFlag) return 0;
  //}

  if (numstep < 0) {  
    procSavevar();
    if (reliabilityFlag) reliabilityStrc->procSavevar();
    if (designoptFlag)   designoptStrc->procSavevar();
    return 0;
  }

  // determine maximum number of transition steps

  if ( ! maxstep ) {
  
    maxstep = procComputeMaxSteps(step);

    if (reliabilityFlag) maxstep = max(maxstep,
				       reliabilityStrc->procComputeMaxSteps(step));
    if (designoptFlag)   maxstep = max(maxstep,
				       designoptStrc->procComputeMaxSteps(step));

    if (maxstep == 0) return 2;
    
    maxstep=min(steplimit,maxstep);

    fprintf(stderr,"\n\n ... Transition in %d steps\n\n",maxstep);
    
    // determine step size for each abstract optimization variable

    procComputeStepSize(maxstep);

    if (designoptFlag)   designoptStrc->procComputeStepSize(maxstep);
    if (reliabilityFlag) reliabilityStrc->procComputeStepSize(maxstep);

  }

  // reset reliability variables
  // WATCH: both optimization and random variables need to be reset
  //        before they can be updated
  
  optpro->optvar->resetvar();

  if (designoptFlag)   designoptStrc->optpro->optvar->resetvar();
  if (reliabilityFlag) reliabilityStrc->optpro->optvar->resetvar();

  // compute new optimization varialbe in transition process
  // WATCH: first update optimization, then random variables

  fprintf(stderr,"\n\n ... Transition Step: %d\n\n",numstep);

  int retval=2;

  if (designoptFlag)   {
    retval = min( retval , designoptStrc->procTransition(numstep,maxstep) );
    retval = min( retval , procTransition(numstep,maxstep) ); 
  }
  else if (reliabilityFlag) {
    retval = min( retval , procTransition(numstep,maxstep) );
    retval = min( retval , reliabilityStrc->procTransition(numstep,maxstep) );
  }
  else
    retval = procTransition(numstep,maxstep);
  
  return retval;
}  

//------------------------------------------------------------------------------

void
Structopt_sd::getMidPointMass(double* totmas, double** midPoint)
{
  int numabs = optpro->optvar->nAbs(); 

  int i;
  for (i=0;i<numabs;i++) {                     //Loop over all opt.var.

    structdom->setoptInf(optInf[i]);          //Set Opt.Var.- Ele. Influence

    //Compute midpoint and mass of 
    //influence area   
    structdom->getMidPointMass(totmas[i], midPoint[i]);
  }
}

//------------------------------------------------------------------------------

double 
Structopt_sd::getforml2norm2()
{
  return optpro->optvar->getFormL2norm2();
}

//------------------------------------------------------------------------------

double 
Structopt_sd::getgradforml2norm2()
{
  return optpro->optvar->getgradFormL2norm2();
}

//------------------------------------------------------------------------------

void 
Structopt_sd::getStressInfo(int typ, int& surface, int& stressIndex)
{
  surface = 1;
  if( typ%10 ==  0 ) surface = 3;   // lower surface
  if( typ%10 ==  1 ) surface = 2;   // middle surface
  if( typ%10 ==  2 ) surface = 1;   // upper surface
    
  stressIndex = -1;
  if( (typ-typ%10) ==  0 ) stressIndex=6;    // von Mises
  if( (typ-typ%10) == 10 ) stressIndex=-1;   // 1. principal
  if( (typ-typ%10) == 20 ) stressIndex=-1;   // 2. principal
  if( (typ-typ%10) == 30 ) stressIndex=0;    // sxx
  if( (typ-typ%10) == 40 ) stressIndex=1;    // syy
  if( (typ-typ%10) == 50 ) stressIndex=2;    // szz

  // temporary 

  if (stressIndex == -1) {
    fprintf(stderr,"ERROR: stress typ not implemented yet\n");
    exit(-1);
  }
}

//------------------------------------------------------------------------------
#ifdef AEROELASTIC	  

void Structopt_sd::sndOptpar(int ivar, int icrit)
{
  // optParam:
  // 0      : number of variables
  // 1      : number of criteria
  // 2      : type of sensitivity analysis
  // 3      : current variable or criteria
  // 4      : analysis id - only need for adjoint
  // 5 - 13 : fluid attributes (numFldAttr)
  // 14- 22 : gradients of fluid attributes (numFldAttr)
  // 23- 29 : flags for fluid criteria: (numFldRes)
  //          23: Fx; 24: Fy; 25: Fz; 
  //          26: Mx; 27: My; 28: Mz; 
  //          29: Boom
  // 30- ?1 : elec info 
  // ?1- ?2 : thermal info
   
  // WATCH:  when merging with electrostatic, make sure that optParam is o.k. 

  int length  = numBaseOptpar + 2*numFldAttr + numFldRes;

  if (ivar < 0 & icrit < 0) anagrdType = 0;

  // add different thermal and elec lengths for analysis and sensitivities

  switch (anagrdType) {
  case 0: 
    length += 3*numElecAttrVar  + 1;
    length += 3*numThermAttrVar + 1;
    break;
  case 1:
  case 2:
    length += 2*numElecAttrVar  + 2;
    length += 2*numThermAttrVar + 2;
    break;
  }

  double *optParam = static_cast<double*>(dbg_alloca(sizeof(double)*length));
  
  int i;
  for (i=0;i<length;i++) optParam[i] = 0.0;
 
  optParam[0] = (double)(optpro->optvar->numabs); 
  optParam[1] = (double)(optpro->numcrit);
  optParam[2] = (double)anagrdType;

  switch (anagrdType)
    {
    case 0:

      optParam[4] = -1.0;

      setFluidvariables(&(optParam[numBaseOptpar]),ivar);

      setFluidcriteria(icrit,&(optParam[numBaseOptpar+2*numFldAttr]));

      setElecvariables(&(optParam[numBaseOptpar+2*numFldAttr+numFldRes]),ivar,
                       anagrdType);

      setThermvariables(
			&(optParam[numBaseOptpar+2*numFldAttr+numFldRes+3*numElecAttrVar+1]),
			ivar,anagrdType);
      break;

    case 1:
      optParam[3] = (double)(ivar+1);
      optParam[4] = -1.0;

      setFluidvariables(&(optParam[numBaseOptpar]),ivar);

      setFluidcriteria(icrit,&(optParam[numBaseOptpar+2*numFldAttr]));

      setElecvariables(&(optParam[numBaseOptpar+2*numFldAttr+numFldRes]),ivar,
                       anagrdType);

      setThermvariables(
			&(optParam[numBaseOptpar+2*numFldAttr+numFldRes+2*numElecAttrVar+2]),
			ivar,anagrdType);
      break;

    case 2:
      optParam[3] = (double)(icrit+1);
      optParam[4] = (double)(activeAnalysis);
     
      setFluidvariables(&(optParam[numBaseOptpar]),ivar);

      setFluidcriteria(icrit,&(optParam[numBaseOptpar+2*numFldAttr]));

      setElecvariables(&(optParam[numBaseOptpar+2*numFldAttr+numFldRes]),ivar,
                       anagrdType);

      setThermvariables(
			&(optParam[numBaseOptpar+2*numFldAttr+numFldRes+2*numElecAttrVar+2]),
			ivar,anagrdType);
      break;
    }

  structdom->sndOptpar((double*) optParam);
  
}

//------------------------------------------------------------------------------

void Structopt_sd::setFluidvariables(double* flgList,int varFlg)
{
  double * fluidattr  = structdom->getFluidAttr();

  int i;
  for (i=0;i<numFldAttr;i++) 
    flgList[i] = fluidattr[i];
  
  if (varFlg >= 0) {
    double * gradfluidattr = structdom->getgradFluidAttr();
    for (i=0;i<numFldAttr;i++) 
      flgList[i+numFldAttr] = gradfluidattr[i];
  }
}  
 
//------------------------------------------------------------------------------

void Structopt_sd::setFluidcriteria(int icrit,double* flgList)
{
  double dummy;

  // if not adjoint sensitivity analysis activate all criteria

  if (icrit < 0) {
    int i;
    for (i=0;i<optpro->numcrit;i++) {
      int dir = optpro->opc[i]->getAeroforce(dummy);
      if(dir >= 0) flgList[dir]=1.0;
    }
    return;
  }
  
  // activate "icrit" criteria only

  if (icrit>optpro->numcrit) {
    fprintf(stderr," ... Error in  setFluidcriteria - stop\n");
    exit(-1);
  }
  
  int dir = optpro->opc[icrit]->getAeroforce(dummy);

  if(dir >= 0) flgList[dir]=1.0;
}  

//------------------------------------------------------------------------------

void Structopt_sd::setElecvariables(double *flgList,int varFlag,int anagrdType)
{

  // if anagrdType = 0, then either a function evaluation or finite differences
  // therefore send all electrostatic attribute STCVAR
  // so that all dependent electrostatic parameters can be updated

  // otherwise, send only the gradients that depend 
  // on the current optimization variable

  int i,pos;
  
  // tell the electrostatic domain how many STCVAR to expect

  flgList[0] = numElecAttrVar;

  if (!numElecAttrVar) return;

  if (anagrdType == 0) {
    for (i = 0; i < numElecAttrVar; i++) {
      pos = i*3;
      flgList[pos+1] = elecattr[pos];
      flgList[pos+2] = elecattr[pos+1];
      flgList[pos+3] = elecattr[pos+2];
    } // for i
  }
  else {
    numCurElecAttrVar = 0;
    if (varFlag >= 0) {

      int svn,ean;
      int numStcvar = optpro->optvar->asTab[varFlag]->size();

      // add stcelv to flag list
      // currently this only works for one element per variable 
      // unlike in structure
      if (optInf[varFlag]->size()) {
        //flgList[1] = double (optInf[varFlag]->getEntry(0));
        flgList[1] = -1.0;
      }
      else {
        flgList[1] = -1.0;
      }

      int numDepElcVar = 0;

      for (i = 0; i < numStcvar; i++) {
        svn = optpro->optvar->asTab[varFlag]->getEntry(i);
        if (optpro->optvar->stcvar[svn]->typ == Optvar::elecattr) {
          pos = numDepElcVar*2+2;
          ean = elecattrList[svn];
          flgList[pos]   = ean;
          flgList[pos+1] = gradelecattr[ean];
          numDepElcVar++;
        }
      } // for i
    
      // add number of dependent electrostatic STCVAR

      flgList[0]        = numDepElcVar;
      numCurElecAttrVar = numDepElcVar;

    } // if varFlag
  } 
  
}  
 
//------------------------------------------------------------------------------

int Structopt_sd::chkOptInfEL()
{

  double norm = 0;
  
  for (int i = 0; i < numGradElecAttrVar; i++) {
    norm += gradelecattr[i]*gradelecattr[i];
  }

  return (norm > 0) ? 6 : 0;
  
}

//------------------------------------------------------------------------------

void Structopt_sd::setThermvariables(double *flgList,int varFlag,int anagrdType)
{

  // if anagrdType = 0, then either a function evaluation or finite differences
  // therefore send all thermal attribute STCVAR
  // so that all dependent thermal parameters can be updated

  // otherwise, send only the gradients that depend 
  // on the current optimization variable

  int i,pos;
  
  // tell the thermal domain how many STCVAR to expect

  flgList[0] = numThermAttrVar;

  if (!numThermAttrVar) return;

  if (anagrdType == 0) {
    for (i = 0; i < numThermAttrVar; i++) {
      pos = i*3;
      flgList[pos+1] = thermattr[pos];
      flgList[pos+2] = thermattr[pos+1];
      flgList[pos+3] = thermattr[pos+2];
    } // for i
  }
  else {
    if (varFlag >= 0) {

      int svn,ean;
      int numStcvar = optpro->optvar->asTab[varFlag]->size();

      // add stcelv to flag list
      // currently this only works for one element per variable 
      // unlike in structure
      if (optInf[varFlag]->size()) {
        flgList[1] = double (optInf[varFlag]->getEntry(0));
      }
      else {
        numCurThermAttrVar = 0;
        flgList[1] = -1.0;
      }

      int numDepThermVar = 0;

      for (i = 0; i < numStcvar; i++) {
        svn = optpro->optvar->asTab[varFlag]->getEntry(i);
        if (optpro->optvar->stcvar[svn]->typ == Optvar::thermattr) {
          pos = numDepThermVar*2+2;
          ean = thermattrList[svn];
          flgList[pos]   = ean;
          flgList[pos+1] = gradthermattr[ean];
          numDepThermVar++;
        }
      } // for i
    
      // add number of dependent thermal STCVAR

      flgList[0]         = numDepThermVar;
      numCurThermAttrVar = numDepThermVar;

    } // if varFlag
  } 
  
}  
 
//------------------------------------------------------------------------------

int Structopt_sd::chkOptInfTH()
{

  double norm = 0;
  
  for (int i = 0; i < numGradThermAttrVar; i++) {
    norm += gradthermattr[i]*gradthermattr[i];
  }

  return (norm > 0) ? 7 : 0;
  
}

#endif  

//------------------------------------------------------------------------------
void Structopt_sd::addAnalysisData( anadata & adata )
{
  extern Domain* domain;

  if ( ! analysisData  && adata.typ < 3 ) {
    analysisData = new StcAnalysisData;
    analysisData->dynFlag   = 0;
    analysisData->dampFlag  = 0; 
    analysisData->transFlag = 0;
  }
    
  switch ( adata.typ )
  {  
    case 0:
       analysisData->dynFlag = 1;
       analysisData->dynTyp  = adata.ival[0];
       break;
    case 1:
       analysisData->dampFlag = 1;
       if (adata.ival[0] > 0) {
         analysisData->dampTyp    = 0;
	 analysisData->dampFrq[0] = adata.ival[0];
	 analysisData->dampFrq[1] = adata.ival[1];
       }
       else 	 
         analysisData->dampTyp    = 1;
       break;
    case 2:   	 
       analysisData->transFlag = 1;
       analysisData->transMin  = adata.ival[0];
       analysisData->transMax  = adata.ival[1];
       analysisData->initBlend = adata.ival[2];
       analysisData->designVel = adata.rval[0];
       analysisData->forceVel  = adata.rval[1];
       break; 
    default:
       fprintf(stderr,"Setting Structural Analysis Data.\n");
       exit(-1);
  }

}  

//------------------------------------------------------------------------------
void Structopt_sd::addOptInf(int ivar, int ele_1, int ele_2, int typ )
{
  extern Domain* domain;
  Domain_opt* optdom = dynamic_cast<Domain_opt*>(domain);
  assert(optdom != 0);
  
  if (! optInfset ) {
    
    int numabs = /*optdom->*/optpro->optvar->nAbs();

      if ( !numabs ) {
        fprintf(stderr,"ERROR: cannot set STCELEVAR before abstract\n");
        fprintf(stderr,"       variables are defined .... exit\n");
        exit(-1);
      }

      initOptInf(numabs);
      optInfset = 1;
   }
   
   int elex = max(ele_2,ele_1+1);
   for(int iele=ele_1; iele < elex;iele++) 
   {
     if (optInf[ivar]==0) 
       optInf[ivar] = new OptActInfo; 
  
     optInf[ivar]->add(optdom->geteleNummap(iele),typ);
   }
}


#endif
