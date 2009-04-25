#include <Structopt.d/Structopt_dec.h>
#include <Utils.d/DistHelper.h>
#include <typeinfo>

template<class Scalar> 
Structopt_dec<Scalar>::Structopt_dec(int tp, GenDecDomain_opt<Scalar>* d, Optpro *o) //: opv(0)
{
  type=tp;
  
  structdom = d;
  optpro    = o;
  
  analysisData = 0;  
  
  staticpros    = 0;
  staticsolvers = 0;
  numStaticProb = 0;

  numvar  = 0 ;
  numFunc = 0 ; 
  
  anagrdType = 0;
  optInfset  = 0;
  
  sendInitDisp  = 1;
}

//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::build(GenDecDomain_opt<Scalar>* domain, Optpro *_optpro) 
{    
  //initialize optimization procedure  
  structdom =  domain;
  optpro    = _optpro;  
  //create structural problem types
  initStructProbSol();

  //set criteria and structural variables
  optpro->optvar->setvalptr(this);
  optpro->optvar->updvar();
  //initialize arrays for analytical sensitivity analysis
  //if ( optpro->checkAnalyticalSA() ) 
    { initAnalyticalSA() ; }      
  
  //open output-files  
  geoSource->openOutputFiles();
  structdom->getDomain()->solInfo().setTimes(1.0e8, 1.0, 0.0);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::initStructProbSol()
{  
  // make subdomains and all the good stuff
  // Makes renumbering, connectivities and dofsets
  structdom->preProcess();

  numAnalysis   = structdom->getNumAnalysis();  
  numStaticProb = structdom->getNumStaticAna();

  // initialize storage arrays for problem types
  if (numStaticProb) {
    staticpros = new GenMultiDomainStatic_opt<Scalar>*[numStaticProb];
    staticsolvers = new StaticSolver_opt<Scalar, GenFetiSolver<Scalar>, GenDistrVector<Scalar>,
      GenMultiDomainPostProcessor_opt<Scalar>, GenMultiDomainStatic_opt<Scalar>,
      GenDistrVector<DComplex > >* [numStaticProb];
  }
  
  
  // allocate and initialize analysis maps  
  analysisMap = new int[numAnalysis];  
  int countStaticProb = 0;
  int countNlnstcProb = 0;
  int countEigenProb  = 0;
  int countDynamProb  = 0;
  int countNLdynProb  = 0;
  int countFVibrProb  = 0;
  
  for(int ina=0;ina<numAnalysis;ina++) 
    {      
      structdom->activateAnalysis(ina,0);  
      
      switch (structdom->probType()) 
	{
	  
	case SolverInfo::Static:
	  analysisMap[ina] = countStaticProb;
	  countStaticProb++;
	  break;
	default:
	  fprintf(stderr, "Problem Type %d is not implemented in Optimization Module\n", structdom->probType());
	  exit(-1);
	}
    }
  
  // initialize individual analysis problems  
  for(int ina=0;ina<numAnalysis;ina++) 
    {    
      structdom->activateAnalysis(ina,0);      
      int apb = analysisMap[ina];
      
      switch (structdom->probType()) 
	{      
	case SolverInfo::Static:
	  staticpros[apb] = new GenMultiDomainStatic_opt<Scalar>(structdom);
	  
	  staticsolvers[apb] = new StaticSolver_opt<Scalar, GenFetiSolver<Scalar>, GenDistrVector<Scalar>,
	    GenMultiDomainPostProcessor_opt<Scalar>, GenMultiDomainStatic_opt<Scalar>,
	    GenDistrVector<DComplex > >
	    (staticpros[apb]);
	  
	  staticsolvers[apb]->optPreSolve();
	  break;
	default:
	  fprintf(stderr, "Problem Type %d is not implemented in Optimization Module\n", structdom->probType());    	  
	  assert(0);
	}
    }
}


//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::initAnalyticalSA()
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

    for(i=0;i<numabs;i++) 
      {
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
template<class Scalar>
void Structopt_dec<Scalar>::getMidPointMass(double* totmas, double** midPoint)
{
  int numabs = optpro->optvar->nAbs(); 

  for(int i=0;i<numabs;i++) 
    {                     //Loop over all opt.var.

      setoptInf(i);          //Set Opt.Var.- Ele. Influence
      //Compute midpoint and mass of 
      //influence area   
      structdom->getMidPointMass(totmas[i], midPoint[i]);
    }
}


//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::func() 
{
  int ifail=1;  
  double funcTime = getTime();

  int ina;
  for (ina=0;ina<numAnalysis;ina++) 
    {    
      structdom->activateAnalysis(ina);      
      int apb = analysisMap[ina];
      
      switch (structdom->probType()) 
	{
	case SolverInfo::Static:
	  staticsolvers[apb]->optSolve();	//Static Analysis
	  ifail=0;
	  break;
	default:
	  assert(0);
	}    
    }  
  // Evaluation of Opt.Crit
  evaluate();                           
  // Count function calls     
  numFunc++;                             
  // Output CPU time for total function evaluation
  filePrint(stderr,"\n ... Time spent in function evaluation: %e sec\n\n",
	    (getTime()-funcTime)/1000.0);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::evaluate() 
{
  for(int i=0;i<optpro->numcrit;i++) 
    { 
      optpro->opc[i]->evaluate(this);
    }
}


//------------------------------------------------------------------------------
template<class Scalar>
double Structopt_dec<Scalar>::getstrainenergy(int* eleList, int listSize, double& time,
					      int anaId) 
{  
  double val;  

  int probType = getAnalysisData(anaId,time);
  
  switch (probType) 
    {
    case SolverInfo::Static:
      val = structdom->getStrainEnergy(*sol,eleList,listSize); 
      break;
    default:
      assert(0);
  }

  return val;
}

//------------------------------------------------------------------------------
template<class Scalar>
int Structopt_dec<Scalar>::getAnalysisData(int anaId, double& time)
{
  // initialize 
  sol = 0;
  vel = 0;
  acc = 0;

  // activate problem type and get pointers
  structdom->activateAnalysis(anaId);  

  int apb = analysisMap[anaId];  
  int lc;  
  switch (structdom->probType()) 
    {
    case SolverInfo::Static:
      lc       = (time < 0) ? 0 : static_cast<int>(time);
      sol      = staticsolvers[apb]->getpsol(lc);
      break;
    default:
      assert(0);
    }
  
  return structdom->probType();
}

//------------------------------------------------------------------------------
template<class Scalar>
double Structopt_dec<Scalar>::getmass(int* eleList, int listSize, int anaId) 
{     
  return structdom->getStructureMass(eleList,listSize);
}


//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::postProcessing(int giter) 
{
  /* FIXIT??? */
  double time=giter-1;

  for(int ina=0; ina<numAnalysis; ina++) 
    { 
      structdom->activateAnalysis(ina);      
      int apb = analysisMap[ina];

      switch (structdom->probType()) 
	{
	case SolverInfo::Static:
	  staticsolvers[apb]->optPostProcessor(time);
	  break;
	default:
	  assert(0);
	}
    }
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::gradadjoint(int icrit, double** gc) 
{
  anagrdType = 2;                       //Set type of SA
  int probType = optpro->opc[icrit]->getAnalysisType();

  switch ( probType ) 
    {
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
      assert(0);
    }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::gradadjointStress(int icrit, double** gc)
{
  // need the analysis ID in sndOptpar, to send to other domains
  int anaId = optpro->opc[icrit]->getAnalysisId();
  activeAnalysis = anaId;

  // Set design velocity = 0
  zerograd();

  // Build dq:du
  optpro->opc[icrit]->graddu(this);

  //Solve adjoint problem: K(-1) * dq:du
  int apb   = analysisMap[anaId];
  int numabs = optpro->optvar->nAbs(); 
  double time = optpro->opc[icrit]->getTime();  //Get load case of criterion
  int lc = (time < 0) ? 0 : static_cast<int>(time);
 
  switch (structdom->probType()) 
    {
    case SolverInfo::Static:  	  
      staticsolvers[apb]->optAdjSolve();
      break;
    default:
      assert(0);
      break;
  }
  // Post-multiply with pseudo-load vector
  for (int i=0;i<numabs;i++) 
    {              //Loop over all opt.var.
      optpro->optvar->updgrad(i,1.0);        //Set Variation of Abs.&Struct. Var
      setoptInf(i);                          //Set Opt.Var.- Ele. Influence
      optpro->opc[icrit]->gradpart(this);    //Determine partial(q)/partial(s)
      double grad1 = optpro->opc[icrit]->grad;
      int stcvarTyp = 0;
      double grad2 = 0.0; 
      //Determine (adjoint(T) * dK:ds * u)
      switch (structdom->probType()) 
	{   
	case SolverInfo::Static:
	  grad2 = staticsolvers[apb]->optGetAdjointPseudo(lc);
	  break;
	default:
	  assert(0);
      }
      gc[icrit][i] = grad1 + grad2;      //Store dq:ds
  }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::zerograd() 
{
  structdom->zeroGrad();
}

//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::getgraddustrainenergy(int* eleList, int listSize, double& time,
						  int anaId) 
{     
  int probType = getAnalysisAdjData(anaId, time);
  switch (structdom->probType()) 
    {
    case SolverInfo::Static:
      adj->zero();  
      structdom->getGradduStrainEnergy(*sol,*adj,eleList,listSize);
      break;
    default:
      assert(0);
    }
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
double Structopt_dec<Scalar>::getgradpartstrainenergy(int* eleList, int listSize, 
						      double& time, int anaId) 
{
  double val;
  int probType = getAnalysisData(anaId,time);
  switch (structdom->probType()) 
    {
    case SolverInfo::Static:
      val = structdom->getGradPartStrainEnergy(*sol,eleList,listSize);
      break;
    default:
      assert(0);
      break;
    }  
  return val;
} 


//------------------------------------------------------------------------------
template<class Scalar>
int Structopt_dec<Scalar>::getAnalysisAdjData(int anaId, double& time, int zeroAll)
{
  // initialize   
  sol  = 0;
  vel  = 0;
  grad = 0;
  adj  = 0;

  // zero all adjoint vectors  
  if(zeroAll) 
    {  
      for(int ina=0;ina<numAnalysis;ina++) 
	{
	  adj = 0;
	  structdom->activateAnalysis(ina);  
	  int apb = analysisMap[anaId];
	  switch (structdom->probType()) 
	    {
	    case SolverInfo::Static:
	      adj   = staticsolvers[apb]->getpadj();
	      break;
	    default:
	      assert(0);
	      break;
	    }
	  
	  if (adj) { adj->zero(); }
	}
    }

  // activate problem type and get pointers
  structdom->activateAnalysis(anaId);  
  int apb = analysisMap[anaId];  
  int lc;  
  switch (structdom->probType()) 
    {
    case SolverInfo::Static:
      lc       = (time < 0) ? 0 : static_cast<int>(time);
      sol      = staticsolvers[apb]->getpsol(lc);
      adj      = staticsolvers[apb]->getpadj();
      break;
    default:
      assert(0);
      break;
  }  
  return structdom->probType();
}


//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::gradadjointGeom(int icrit, double** gc) 
{
  zerograd();                                  //Set design velocity = 0

  int i;
  int numabs = optpro->optvar->nAbs(); 

  double time = optpro->opc[icrit]->getTime();      //Get load case of criterion
  int lc = (time < 0) ? 0 : static_cast<int>(time);
  
  for (i=0;i<numabs;i++) {                          //Loop over all opt.var.

    optpro->optvar->updgrad(i,1.0);                //Set Variation of Abs.&Struct. Var

    setoptInf(i);                                  //Set Opt.Var.- Ele. Influence
     
    optpro->opc[icrit]->gradpart(this);            //Determine partial(q)/partial(s)
       
    gc[icrit][i] = optpro->opc[icrit]->grad;       //Store dq:ds
  }
}

//------------------------------------------------------------------------------
template<class Scalar>
double Structopt_dec<Scalar>::getgradmass(int* eleList, int listSize, int anaId)
{
  return structdom->getGradStructureMass(eleList, listSize);
}


//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::cleanup() 
{
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::setoptInf(int ivar)
{
  structdom->setoptInf(optInf[ivar]);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::graddirect(int ivar, double** gc, int numAC, int* listAC) 
{
  anagrdType = 1;                       //Set type of SA
  // evaluate gradients of deterministic criteria 
  for (int ina=0;ina<numAnalysis;ina++) 
    { 
      structdom->activateAnalysis(ina); 
      int apb = analysisMap[ina]; 
      switch (structdom->probType()) 
	{
	case SolverInfo::Static:
	  staticsolvers[apb]->optSolveDuds();
	  break;
	default:
	  assert(0);
	}
    }
  gradcrit(ivar,gc,numAC,listAC);
  return;
}

//------------------------------------------------------------------------------
template<class Scalar>
double Structopt_dec<Scalar>::getgradstrainenergy(int* eleList, int listSize, double& time,
						  int anaId) 
{
  double val;
  int probType = this->getAnalysisGradData(anaId, time);
  switch (probType) 
    {      
    case SolverInfo::Static:
      val = structdom->getGradStrainEnergy(*sol,*grad,eleList,listSize);
      break;
    default:
      assert(0);
    }  
  return val;
} 

//------------------------------------------------------------------------------
template<class Scalar>
int Structopt_dec<Scalar>::getAnalysisGradData(int anaId, double& time)
{
  // initialize 
  sol  = 0;
  vel  = 0;
  grad = 0;
  adj  = 0;
  // activate problem type and get pointers
  structdom->activateAnalysis(anaId);  
  int apb = analysisMap[anaId];  
  int lc;  
  switch (structdom->probType()) 
    {
    case SolverInfo::Static:
      lc       = (time < 0) ? 0 : static_cast<int>(time);
      sol	 = staticsolvers[apb]->getpsol(lc);
      grad     = staticsolvers[apb]->getpgrad(lc);
      break;
    default:
      assert(0);
    }
  return structdom->probType();
}


//------------------------------------------------------------------------------
template<class Scalar>
double Structopt_dec<Scalar>::getdisp(int node, int typ, int dva, double& time,
				      int anaId)
{    
  double val;
  int probType = getAnalysisData(anaId,time);
  int apb = analysisMap[anaId];   
  switch (probType) 
    {
    case SolverInfo::Static:
      val = ScalarTypes::Real(structdom->getNodalDisp(*sol, node, typ, true));
      break;
    default:
      fprintf(stderr,"Error: nodal displacements|velocities|accelerations can only\n");
      fprintf(stderr,"       be evaluated for static and dynamic problems ... stop\n\n");
      exit(-1);
    }
  return val;
}				       


//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::getgraddudisp(int node, int typ, int dva, double& time,
					  int anaId) 
{
  Scalar disp;
  int probType = getAnalysisAdjData(anaId,time,1);
  switch (probType) 
    {
    case SolverInfo::Static:
      adj->zero();
      structdom->getGradduNodalDisp(*sol,*adj,node,typ, true);
      /**adj *= ScalarTypes::d_doublify(structdom->getNodalDisp(*sol, node, typ),1.0);*/
      break;
    default:
      fprintf(stderr,"Error: derivatives of nodal displacements can only be\n");
      fprintf(stderr,"       evaluated for static and dynamic problems\n");
      fprintf(stderr,"...... stop\n");
      exit(-1);
    }
}		       

//------------------------------------------------------------------------------
template<class Scalar>
double Structopt_dec<Scalar>::getgraddisp(int node, int typ, int dva, double& time,
					  int anaId)
{
  double val;
  int probType = getAnalysisGradData(anaId,time);
  switch (probType) 
    {
    case SolverInfo::Static:
      val = ScalarTypes::Real(structdom->getGradNodalDisp(*sol, *grad, node, typ, true));
      break;
    default:
      fprintf(stderr,"Error: derivatives of nodal displacements can only be\n");
      fprintf(stderr,"       evaluated for static and dynamic problems\n");
      fprintf(stderr,"...... stop\n");
      exit(-1);
    }
  return val;
}				       

//------------------------------------------------------------------------------
template<class Scalar>
double Structopt_dec<Scalar>::getdisplevel(int aFlag, int quoFlag, int size, 
					   int* nodeid, int* distyp, double* refVal, 
					   double* difVal, double& powFac, double& time, 
					   int anaId)
{
  double val;
  int probType = getAnalysisData(anaId,time);
  int apb = analysisMap[anaId];   
  switch (probType) 
    {
    case SolverInfo::Static:
      val = structdom->getDisplevel(*sol,  aFlag, quoFlag, size, nodeid, distyp, 
				    refVal, difVal, powFac);
      break;
    default:
      fprintf(stderr,"Error: nodal displacements|velocities|accelerations can only\n");
      fprintf(stderr,"       be evaluated for static and dynamic problems ... stop\n\n");
      exit(-1);
    }
  return val;
}				       

//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::getgraddudisplevel(int aFlag, int quoFlag, int size, 
					       int* nodeid, int* distyp, double* refVal, 
					       double* difVal, double& powFac, double& time,
					       int anaId)
{
  int probType = getAnalysisAdjData(anaId, time);
  int apb = analysisMap[anaId];   
  switch (probType) 
    {
    case SolverInfo::Static:
      adj->zero();  
      structdom->getGradduDisplevel(*sol, *adj, aFlag, quoFlag, size, nodeid, distyp, 
				    refVal, difVal, powFac);
      break;
    default:
      fprintf(stderr,"Error: nodal displacements|velocities|accelerations can only\n");
      fprintf(stderr,"       be evaluated for static and dynamic problems ... stop\n\n");
      exit(-1);
    }
  return;
}


//------------------------------------------------------------------------------
template<class Scalar>
double Structopt_dec<Scalar>::getgraddisplevel(int aFlag, int quoFlag, int size, 
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
      val = structdom->getGradDisplevel(*sol,  *grad, aFlag, quoFlag, size, nodeid, distyp, 
      					refVal, difVal, powFac);
      break;
    default:
      fprintf(stderr,"Error: nodal displacements|velocities|accelerations can only\n");
      fprintf(stderr,"       be evaluated for static and dynamic problems ... stop\n\n");
      exit(-1);
    }
  return val;
}				       


//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::addOptInf(int ivar, int ele_1, int ele_2, int typ )
{
  if(!optInfset ) 
    {      
      int numabs = /*optdom->*/optpro->optvar->nAbs();
      if (!numabs) 
	{
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
       {
	 optInf[ivar] = new OptActInfo;
       }     
     optInf[ivar]->add(geoSource->glToPackElem(iele),typ);
   }
}

//------------------------------------------------------------------------------
template<class Scalar>
void Structopt_dec<Scalar>::initOptInf(int numvar)
{
  optInf = new OptActInfo*[numvar];
  
  for(int ivar=0; ivar<numvar; ++ivar)
    { optInf[ivar] = new OptActInfo; }
}
