#include <Structopt.d/Structopt.h>
#include <Math.d/mathUtility.h>
#include <Element.d/Element.h>

//------------------------------------------------------------------------------

Structopt::Structopt():opc(0),opv(0) {

    structdom = 0;
    optpro    = 0;
    optsol    = 0;
    
    analysisData = 0;  
    
    staticpro    = 0;
    staticsolver = 0;
    
    eigenpro    = 0;
    eigensolver = 0;

    dynpro    = 0;
    dynsolver = 0;
    
    numcrit = 0 ;
    numvar  = 0 ;
    numfunc = 0 ; 

    anagrdType = 0;
}
 
//------------------------------------------------------------------------------

void Structopt::build( Domain *domain, Optpro *_optpro, Optsol *_optsol) {
    
    structdom =  domain;
    optpro    = _optpro;
    optsol    = _optsol;
 
    //implemented problem types: single domain static  analysis
    //                           single domain modal   analysis
    //                           single domain dynamic analysis

    //initialize optimization procedure

    //watch ! last modal analysis since this can be subproblem of
    //        static or dynamic problems 

    switch (structdom->probType()) {

      case SolverInfo::Static:
        
	staticpro = new SingleDomainStatic<double,Vector,Solver>(structdom);  
	
        staticsolver = new StaticSolver<Solver, Vector, 
 		           SingleDomainPostProcessor<double, Vector, Solver>, 
		           SingleDomainStatic<double,Vector,Solver> >
	 		   (staticpro);
      
        staticsolver->optPreSolve();      

	break;
      case SolverInfo::Dynamic:

        dynpro = new SingleDomainDynamic(structdom);

        dynsolver = new DynamicSolver <DynamMat, Vector,
                        SDDynamPostProcessor, SingleDomainDynamic,double>
                        (dynpro);
      
        dynsolver->optPreSolve(this);      
	
	break;
      case SolverInfo::Modal:
      
 	eigenpro = new SingleDomainEigen(structdom);
	
        eigensolver =  EigenSolver<DynamMat, Vector, VectorSet, 
	                           SDEigenPostProcessor,SingleDomainEigen>
			           ::buildEigenSolver(eigenpro);

        eigensolver->optPreSolve();        

        break;
      default:
      
        fprintf(stderr,"Problem Type not implemented in Optimization Module");
	exit(-1);
    }

    //structdom->openOutputFiles();                  //open output-files
    geoSource->openOutputFiles();

    //set criteria and structural variables

    numcrit = optpro->numcrit;
    numvar  = optpro->optvar->numstc;

    int i;
    for(i=0;i<numcrit;i++) { opc[i]=optpro->opc[i]; }

    for(i=0;i<numvar;i++) { opv[i]=optpro->optvar->stcvar[i]; 
                            opv[i]->setvalptr(this); }

    //initialize arrays for analytical sensitivity analysis

    int igrad = optsol->getanagrad();  
    
    if(igrad) {
    
      structdom->buildOptGrad();

      for(i=0;i<numvar;i++) { opv[i]->setgradptr(this); }
      
      //build optimization variable - elmement influence table
      
      int numabs = optpro->optvar->numabs; 
      for(i=0;i<numabs;i++) {
        optpro->optvar->updgrad(i,1.0);
	structdom->buildOptInf(i,numabs);
      }
    }

    //initialize array to save sensitivities of aeroforces
    //                     (necessary to control blending)

    aeroact  = 1;
    
    if (aeroact ) {
      aeroforce = new double[3];
      aerosens  = new double*[numvar];
      int i;
      for(i=0;i<numvar;i++) aerosens[i]=new double[3];
    }
}

//------------------------------------------------------------------------------

void Structopt::cleanup() {

 switch (structdom->probType()) {

    case SolverInfo::Static:
      break;
    case SolverInfo::Modal:
      break;
    case SolverInfo::Dynamic:      
#ifdef AEROELASTIC
      dynpro->cmdCom(-1);
#endif      
      break;
    default:
      fprintf(stderr,"Problem Type not implemented in Optimization Module");
      exit(-1);
  }
}

//------------------------------------------------------------------------------

double * Structopt::getptreleattr( int& loc1, int& loc2 ) {

      double *p;
	  
      StructProp *sprop = structdom->getStructProps();
	   	      
      switch (loc2) {
	 case 1: 
	    p=&(sprop[loc1].E); 
	    break;
	 case 3: 
	    p=&(sprop[loc1].rho); 
	    break;
         case 6: 
	    p=&(sprop[loc1].eh); 
	    break;
	 default: 
	    p=0; 
      }
      return p;
}

//------------------------------------------------------------------------------

double * Structopt::getptrnodattr( int& loc1, int& loc2 ) {

      double *p=0;
	  
      return p;
}

//------------------------------------------------------------------------------

double * Structopt::getptrnodcord( int& loc1, int& loc2 ) {

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
double * Structopt::getptrcomposite( int& loc1, int& loc2 ) 
{
  double* p=0;
  if(loc2 >= 0)
    {
      // layers
      int layer = loc2/100000;
      int typ   = loc2 - 100000*layer;
      
      LayInfo* layinfo = structdom->getLayerInfo(loc1);
      
      int numlayer = layinfo->nLayers();
      
      if (numlayer < layer) 
	{
	  fprintf(stderr," *** ERROR: Layer %d to optimize does not exist\n",
		  layer);
	  exit(-1);
	}
      
      switch (typ) 
	{      
	case 1: 
	  p=&(layinfo->data[layer-1][0]); 
	  break;
	case 2: 
	  p=&(layinfo->data[layer-1][1]); 
	  break;
	case 4: 
	  p=&(layinfo->data[layer-1][3]); 
	  break;
	case 7: 
	  p=&(layinfo->data[layer-1][6]); 
	  break;
	case 8: 
	  p=&(layinfo->data[layer-1][7]); 
	  break;
	case 9: 
	  p=&(layinfo->data[layer-1][8]); 
	  break;
	default: 
	  p=0; 
	}
    }
  else
    {
      // coefficients
      int i = (-loc2-1)/6;
      int j = (-loc2-1)%6;
      assert(loc1 >= 0 && loc1 < geoSource->numCoefData);
      CoefData* pCfDta = geoSource->coefData[loc1];
      p = pCfDta->c[std::min(i,j)+6*std::max(i,j)];
    }
  
  return p;
}

//------------------------------------------------------------------------------

double * Structopt::getptrnodalforce( int& loc1, int& loc2 ) {

      int numNBC = structdom->nNeumann();  
      BCond* nbc = structdom->getNBC();  

      double *p=0;
 
      int ibc;
      for(ibc=0;ibc<numNBC;ibc++)
        if (nbc[ibc].nnum == loc1 && nbc[ibc].dofnum)
	{
	  p = &(nbc[ibc].val);
          break;
	}
	
      if ( !p )	
      {	    	  
         fprintf(stderr,"ERROR:  Optimization variable: Force in Node %d for"
                        " DOF %d could not be set",loc1,loc2);
         exit(-1);
      }
      
      return p;	     
}

//------------------------------------------------------------------------------

double * Structopt::getptrgradeleattr( int& loc1, int& loc2 ) {

      double *p;
	  
      StructProp *gradsprop = structdom->getGradStructProps();
	   	      
      switch (loc2) {
	 case 1: 
	    p=&(gradsprop[loc1].E); 
	    break;
	 case 3: 
	    p=&(gradsprop[loc1].rho); 
	    break;
         case 6: 
	    p=&(gradsprop[loc1].eh); 
	    break;
	 default: 
	    p=0; 
      }
      return p;
}

//------------------------------------------------------------------------------

double * Structopt::getptrgradnodattr( int& loc1, int& loc2 ) {

      double *p=0;
	  
      return p;
}

//------------------------------------------------------------------------------

double * Structopt::getptrgradnodcord( int& loc1, int& loc2 ) {

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
double * Structopt::getptrgradcomposite( int& loc1, int& loc2 ) 
{
  double* p = 0;
  if(loc2 >= 0)
    {
      // layers
      int layer = loc2/100000;
      int typ   = loc2 - 100000*layer;
      
      LayInfo * layinfo = structdom->getLayerInfo(loc1);      
      switch (typ) 
	{	   	   
	case 1: 
	  p=&(layinfo->grad[layer-1][0]); 
	  break;
	case 2: 
	  p=&(layinfo->grad[layer-1][1]); 
	  break;
	case 4: 
	  p=&(layinfo->grad[layer-1][3]); 
	  break;
	case 7: 
	  p=&(layinfo->grad[layer-1][6]); 
	  break;
	case 8: 
	  p=&(layinfo->grad[layer-1][7]); 
	  break;
	case 9: 
	  p=&(layinfo->grad[layer-1][8]); 
	  break;
	default: 
	  p=0; 
	}
    }
  else
    {
      // coefficients
      int i = (-loc2-1)/6;
      int j = (-loc2-1)%6;
      assert(loc1 >= 0 && loc1 < geoSource->numCoefData);
      //std::map<int, CoefData>::iterator it = static_cast<GeoSource_opt*>(geoSource)
      //CoefData* pCfDta = geoSource->coefData[loc1];
      //p = pCfDta->c[std::min(i,j)+6*std::max(i,j)];      
    }
  return p;
}

//------------------------------------------------------------------------------

double * Structopt::getptrfluidattr( int& loc1) {

      double *p=0;

#ifdef AEROELASTIC	  
      double *fluidattr = structdom->getFluidAttr();
	   	      
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
	 default: 
            fprintf(stderr,"Fluid attribute not implemented\n");
	    exit(-1);
      }
#endif

      return p;
}

//------------------------------------------------------------------------------

double * Structopt::getptrgradfluidattr( int& loc1) {

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
	 default: 
            fprintf(stderr,"Fluid attribute not implemented\n");
	    exit(-1);
      }
#endif
      return p;
}

//------------------------------------------------------------------------------

void Structopt::func() {

#ifdef AEROELASTIC
   sndOptpar(-1,-1);  
#endif

  switch (structdom->probType()) {

    case SolverInfo::Static:
      staticsolver->optSolve();         //Structural Analysis
      break;
    case SolverInfo::Modal:
      eigensolver->optSolve();          //Eigenvalue Analysis
      break;
    case SolverInfo::Dynamic:      
      dynsolver->optSolve(this);        //Dynamic Analysis
      break;
    default:
      fprintf(stderr,"Problem Type not implemented in Optimization Module");
      exit(-1);
  }
  
  evaluate();                                    //Evaluation of Opt.Crit
}

//------------------------------------------------------------------------------

void Structopt::evaluate() {

  // get optimization results from fluid

#ifdef AEROELASTIC

  // save aerodynamic forces of fluid solution
  
  structdom->saveAeroForce();

#endif

  int i;
  for (i=0;i<numcrit;i++) { 

     opc[i]->evaluate(this);
  }
}

//------------------------------------------------------------------------------

void Structopt::postProcessing(int giter ) {

  double time=giter;
    
  switch (structdom->probType()) {

    case SolverInfo::Static:
      staticsolver->optPostProcessor(time);
      break;
    case SolverInfo::Modal:
      eigensolver->optPostProcessor(time);
      break;
    case SolverInfo::Dynamic:
      //dynsolver->optPostProcessor(time);
      break;
    default:
      fprintf(stderr,"Problem Type not implemented in Optimization Module");
      exit(-1);
  }
}

//------------------------------------------------------------------------------

double Structopt::getstrainenergy() {
  
  double * bcx;
  Vector * sol;

  SparseMatrix     * gStiff   = 0;
  FullSquareMatrix * kelArray = 0;
    
  switch (structdom->probType()) {

    case SolverInfo::Static:
      sol      = staticsolver->getpsol();
      bcx      = staticpro->getbc(); 
      kelArray = staticpro->getkelArray();     
      break;
    case SolverInfo::Dynamic:
      sol    = dynsolver->getpDis();
      gStiff = dynpro->getpK(dynsolver->getpOps());
      bcx    = dynpro->boundaryValue();      
      break;
    default:
      fprintf(stderr,"Error: strain energy can only be evaluated for\n");
      fprintf(stderr,"       static and dynamic problems ...... stop\n");
      exit(-1);
  }

  double val=structdom->getStrainEnergy(*sol, bcx, gStiff, kelArray);
     
  return val;
}

//------------------------------------------------------------------------------

double Structopt::getgradstrainenergy() {

  double * bcx;
  Vector * sol;
  Vector * grad;

  FullSquareMatrix * kelArray = 0;
     
  switch (structdom->probType()) {

    case SolverInfo::Static:
      sol      = staticsolver->getpsol();
      grad     = staticsolver->getpgrad();
      bcx      = staticpro->getbc(); 
      kelArray = staticpro->getkelArray();     
      break;
    case SolverInfo::Dynamic:
      sol      = dynsolver->getpDis();
      grad     = dynsolver->getpgDis();
      bcx      = dynpro->boundaryValue(); 
      kelArray = 0;     
      break;
    default:
      fprintf(stderr,"Error: derivatives of strain energy can only be \n");
      fprintf(stderr,"       evaluated for static and dynamic problems \n");
      fprintf(stderr,"...... stop\n");
      exit(-1);
  }

  double val=structdom->getGradStrainEnergy(*sol, *grad, bcx, kelArray);
     
  return val;
} 
//------------------------------------------------------------------------------

void Structopt::getgraddustrainenergy() {
     
  double * bcx;
  Vector * sol;
  Vector * adj;

  FullSquareMatrix * kelArray = 0;

  switch (structdom->probType()) {

    case SolverInfo::Static:
      sol      = staticsolver->getpsol();
      adj      = staticsolver->getpadj();
      bcx      = staticpro->getbc(); 
      kelArray = staticpro->getkelArray();     
      break;
    case SolverInfo::Dynamic:
      sol      = dynsolver->getpDis();
      adj      = dynsolver->getpadj();
      bcx      = dynpro->boundaryValue(); 
      kelArray = 0;     
      break;
    default:
      fprintf(stderr,"Error: derivatives of strain energy can only be \n");
      fprintf(stderr,"       evaluated for (quasi) static problems ...\n");
      fprintf(stderr,"       stop\n");
      exit(-1);
  }
    
  adj->zero();  

  structdom->getGradduStrainEnergy(*sol, *adj, bcx, kelArray);
}

//------------------------------------------------------------------------------

double Structopt::getmass() {
     
     double val = structdom->getStructureMass();
          
     return val;
}

//------------------------------------------------------------------------------

double Structopt::getgradmass() {
     
     double val = structdom->getGradStructureMass();
          
     return val;
}

//------------------------------------------------------------------------------

double Structopt::getfrequency( int ieig ) {
          
  Vector * eigval;

  switch (structdom->probType()) {

    case SolverInfo::Modal:
      eigval = eigensolver->getpeigval();
      break;
    default:
      fprintf(stderr,"Error: eigenfrequencies can only be evaluated\n");
      fprintf(stderr,"       for eigen problems ...... stop\n");
      exit(-1);
  }
     
  double val = sqrt( (*eigval)[ieig] ) / (2.0*3.14159265);     
  
  return val;
}

//------------------------------------------------------------------------------

double Structopt::getgradfrequency( int ieig ) {

  Vector * eigval;
  Vector * eigrad;

  switch (structdom->probType()) {

    case SolverInfo::Modal:
      eigval = eigensolver->getpeigval();
      eigrad = eigensolver->getpeigrad();
      break;
    default:
      fprintf(stderr,"Error: derivatives of eigenfrequencies can only be\n");
      fprintf(stderr,"       evaluated for eigen problems ...... stop\n");
      exit(-1);
  }

  double grad = (*eigrad)[ieig] / ( 4.0*3.14159265*sqrt( (*eigval)[ieig] ) );
    
  return grad;
}
    
//------------------------------------------------------------------------------

double Structopt::getcontrolcost() {
          
  double val;

  switch (structdom->probType()) {

    case SolverInfo::Modal:
      val = eigensolver->getControlCost();
      break;
    default:
      fprintf(stderr,"Error: Control cost can only be evaluated\n");
      fprintf(stderr,"       for eigen problems ...... stop\n");
      exit(-1);
  }
  
  return val;
}

//------------------------------------------------------------------------------

double Structopt::getkineticenergy() {

  Vector * vel;

  SparseMatrix * gMass = 0;

  switch (structdom->probType()) {

    case SolverInfo::Dynamic:
      vel    = dynsolver->getpVel();
      gMass  = dynpro->getpM(dynsolver->getpOps());
      break;
    default:
      fprintf(stderr,
        "Error: kinetic energy make sens for dynamic problems only. stop\n");
         exit(-1);
  }

  double val=structdom->getKineticEnergy(*vel, gMass);
     
  return val;
}

//------------------------------------------------------------------------------

double Structopt::getgradkineticenergy() {

  double val=0.0;

  fprintf(stderr,"Error: analytical derivatives of kinetic energy are not\n");
  fprintf(stderr,"	 implemented yet .... stop\n");
  exit(-1);
     
  return val;
} 

//------------------------------------------------------------------------------

double Structopt::getdampingenergy() {

  Vector * vel;

  SparseMatrix * gDamp = 0;

  switch (structdom->probType()) {

    case SolverInfo::Dynamic:
      vel    = dynsolver->getpVel();
      gDamp  = dynpro->getpC(dynsolver->getpOps());
      break;
    default:
      fprintf(stderr,
        "Error: kinetic energy make sens for dynamic problems only. stop\n");
         exit(-1);
  }

  double val=structdom->getDampingEnergy(*vel, gDamp);
     
  return val;
}

//------------------------------------------------------------------------------

double Structopt::getgraddampingenergy() {

  double val=0.0;

  fprintf(stderr,"Error: analytical derivatives of damping energy are not\n");
  fprintf(stderr,"	 implemented yet .... stop\n");
  exit(-1);
     
  return val;
} 

//------------------------------------------------------------------------------

double Structopt::getaeroforce( int idir ) {

  if ( structdom->probType() != SolverInfo::Dynamic ) {

    fprintf(stderr,
    "Error: Aeroforces only for dynamic problems in aeroelasticity. stop\n");
    exit(-1);
  }

  double val=0.0;

#ifdef AEROELASTIC
  val=structdom->getAeroForce(idir);
#else
  fprintf(stderr,
  "Error: Aeroforces only for dynamic problems in aeroelasticity. stop\n");
  exit(-1);
#endif      

  return val;
}

//------------------------------------------------------------------------------

double Structopt::getgradaeroforce( int idir ) {

  if ( structdom->probType() != SolverInfo::Dynamic ) {

    fprintf(stderr,
    "Error: Aeroforces only for dynamic problems in aeroelasticity. stop\n");
    exit(-1);
  }

  double val=0.0;

#ifdef AEROELASTIC
  if (anagrdType == 1)
   val=structdom->getgradAeroForce(idir);
#else
  fprintf(stderr,
  "Error: Derivative of Aeroforces only for dynamic problems in aeroelasticity. stop\n");
  exit(-1);
#endif      
     
  return val;
} 
//------------------------------------------------------------------------------

void Structopt::getgradduaeroforce() {

  Vector * adj;

  switch (structdom->probType()) {

    case SolverInfo::Dynamic:
      adj      = dynsolver->getpadj();
      break;
    default:
      fprintf(stderr,
      "Error: Aeroforces only for dynamic problems in aeroelasticity. stop\n");
      exit(-1);
  }
    
  adj->zero();  
} 

//------------------------------------------------------------------------------

double Structopt::getaeromoment( int node, int idir ) {

  if ( structdom->probType() != SolverInfo::Dynamic ) {
    fprintf(stderr,
    "Error: Aeroforces only for dynamic problems in aeroelasticity. stop\n");
    exit(-1);
  }

  double val=0.0;

#ifdef AEROELASTIC
  Vector * aeroForce = dynsolver->getaeroForce();
  val=structdom->getAeroMoment(*aeroForce,node,idir);
#else
  fprintf(stderr,
  "Error: Aeroforces only for dynamic problems in aeroelasticity. stop\n");
  exit(-1);
#endif      

  return val;
}

//------------------------------------------------------------------------------

double Structopt::getgradaeromoment( int inode, int idir ) {

  double val=0.0;
     
  fprintf(stderr,"Error: analytical derivatives of aero moments are not\n");
  fprintf(stderr,"	 implemented yet .... stop\n");
  exit(-1);
     
  return val;
} 

//------------------------------------------------------------------------------

double Structopt::getnodstr(int inode, int typ) {

  int surface = 1;
  if(typ == 0 || typ == 10 || typ == 20) surface = 3;
  if(typ == 1 || typ == 11 || typ == 21) surface = 2;
  if(typ == 2 || typ == 12 || typ == 22) surface = 1;
    
  int stressIndex = 6;
  if(typ < 20) stressIndex=6;
  if(typ < 10) stressIndex=6;

  double * bcx;
  Vector * sol;    

  switch (structdom->probType()) {

    case SolverInfo::Static:
      sol = staticsolver->getpsol();
      bcx = staticpro->getbc(); 
      break;
    case SolverInfo::Dynamic:
      sol = dynsolver->getpDis();
      bcx = dynpro->boundaryValue();      
      break;
    default:
      fprintf(stderr,"Error: nodal stress can only be evaluated for\n");
      fprintf(stderr,"       static and dynamic problems ..... stop\n");
      exit(-1);
  }

 
  double val=structdom->getNodalStressStrain(*sol, bcx, inode,
                                              stressIndex, surface);
  return val;
}
//------------------------------------------------------------------------------

double Structopt::getgradnodstr(int inode, int typ) {

  int surface = 1;
  if(typ == 0 || typ == 10 || typ == 20) surface = 3;
  if(typ == 1 || typ == 11 || typ == 21) surface = 2;
  if(typ == 2 || typ == 12 || typ == 22) surface = 1;
    
  int stressIndex = 6;
  if(typ < 20) stressIndex=6;
  if(typ < 10) stressIndex=6;
  
  double * bcx;
  Vector * sol;
  Vector * grad;
    
  switch (structdom->probType()) {

    case SolverInfo::Static:
      sol  = staticsolver->getpsol();
      grad = staticsolver->getpgrad();
      bcx  = staticpro->getbc(); 
      break;
   case SolverInfo::Dynamic:
      sol      = dynsolver->getpDis();
      grad     = dynsolver->getpgDis();
      bcx      = dynpro->boundaryValue(); 
      break;
    default:
      fprintf(stderr,"Error: derivatives of nodal stress can only be\n");
      fprintf(stderr,"       evaluated for static and dynamic problems\n");
      fprintf(stderr,"       .... stop\n");
      exit(-1);
  }


  double val=structdom->getGradNodalStressStrain(*sol, *grad, bcx, inode,
                                                  stressIndex, surface);
  return val;
}

//------------------------------------------------------------------------------

void Structopt::getgraddunodstr(int inode, int typ) {

  int surface = 1;
  if(typ == 0 || typ == 10 || typ == 20) surface = 3;
  if(typ == 1 || typ == 11 || typ == 21) surface = 2;
  if(typ == 2 || typ == 12 || typ == 22) surface = 1;
    
  int stressIndex = 6;
  if(typ < 20) stressIndex=6;
  if(typ < 10) stressIndex=6;
    
  double * bcx;
  Vector * sol;
  Vector * adj;  
    
  switch (structdom->probType()) {

    case SolverInfo::Static:
      sol = staticsolver->getpsol();
      adj = staticsolver->getpadj();
      bcx = staticpro->getbc(); 
      break;
    case SolverInfo::Dynamic:
      sol      = dynsolver->getpDis();
      adj      = dynsolver->getpadj();
      bcx      = dynpro->boundaryValue(); 
      break;
    default:
      fprintf(stderr,"Error: derivatives of nodal stress can only be\n");
      fprintf(stderr,"       evaluated for (quasi) static problems ...\n");
      fprintf(stderr,"       stop\n");
      exit(-1);
  }

  adj->zero();

  structdom->getGradduNodalStressStrain(*sol, *adj, bcx,
                                        inode, stressIndex, surface);
}

//------------------------------------------------------------------------------

double Structopt::getdisp( int node, int typ, int dva) {
    
  double * bcx;
  Vector * sol;  
    
  switch (structdom->probType()) {

    case SolverInfo::Static:

      if ( dva ) {
        fprintf(stderr,
	  "Error: no velocities or accerlations in statics .... stop !\n");
        exit(-1);
      }	   
      sol = staticsolver->getpsol();
      bcx = staticpro->getbc(); 

      break;
    case SolverInfo::Dynamic:

      switch (dva) {
        
	case 0:
          sol = dynsolver->getpDis();
          bcx = dynpro->boundaryValue();      
          break;
	case 1:
	  sol = dynsolver->getpVel();
	  bcx = dynpro->boundaryVeloc();
          break;
	case 2:
	  sol=dynsolver->getpAcc();
      }
      break;
    default:

      fprintf(stderr,
         "Error: nodal displacements|velocities|accerlations can only be \n");
      fprintf(stderr,
         "be evaluated for static and dynamic problems ............. stop\n");
      exit(-1);
  }

  double val=structdom->getNodalDisp(*sol, bcx, node, typ);

  return val;
}				       
				       
//------------------------------------------------------------------------------

double Structopt::getgraddisp( int node, int typ, int dva) {
    
  double * bcx;
  Vector * sol;
  Vector * grad;

  if ( dva ) {
    fprintf(stderr,
      "Error: no velocities or accerlations derivatives implemented. stop !\n");
     exit(-1);
  }	   

  switch (structdom->probType()) {

    case SolverInfo::Static:
      sol  = staticsolver->getpsol();
      grad = staticsolver->getpgrad();
      bcx  = staticpro->getbc(); 
      break;
    case SolverInfo::Dynamic:
      sol      = dynsolver->getpDis();
      grad     = dynsolver->getpgDis();
      bcx      = dynpro->boundaryValue(); 
      break;
    default:
      fprintf(stderr,"Error: derivatives of nodal displacements can only be\n");
      fprintf(stderr,"       evaluated for static and dynamic problems\n");
      fprintf(stderr,"...... stop\n");
      exit(-1);
  }

  double val=structdom->getGradNodalDisp(*sol, *grad, bcx, node, typ);

  return val;
}				       
				       
//------------------------------------------------------------------------------

void Structopt::getgraddudisp( int node, int typ, int dva) {
    
  double * bcx;
  Vector * sol;
  Vector * adj;

  if ( dva ) {
    fprintf(stderr,
      "Error: no velocities or accerlations derivatives implemented. stop !\n");
     exit(-1);
  }	   

  switch (structdom->probType()) {

    case SolverInfo::Static:
      sol = staticsolver->getpsol();
      adj = staticsolver->getpadj();
      bcx = staticpro->getbc(); 
      break;
    case SolverInfo::Dynamic:
      sol      = dynsolver->getpDis();
      adj      = dynsolver->getpadj();
      bcx      = dynpro->boundaryValue(); 
      break;
    default:
      fprintf(stderr,"Error: derivatives of nodal displacements can only be\n");
      fprintf(stderr,"       evaluated for (quasi) static problems ........\n");
      fprintf(stderr,"       stop\n");
      exit(-1);
  }

  adj->zero();
    
  structdom->getGradduNodalDisp(*sol, *adj, bcx, node, typ);
}				       
  
//------------------------------------------------------------------------------

void Structopt::getgraddumass() {

  Vector * adj;

  switch (structdom->probType()) {

    case SolverInfo::Static:
      adj = staticsolver->getpadj();
      break;
    case SolverInfo::Dynamic:
      adj = dynsolver->getpadj();
      break;
    default:
      fprintf(stderr,"Error: derivative of mass only for (quasi) static \n");
      fprintf(stderr,"       problems implemented ........\n");
      fprintf(stderr,"       stop\n");
      exit(-1);
  }

  adj->zero();
}
//------------------------------------------------------------------------------

void Structopt::graddirect(int ivar) {

  anagrdType = 1;                       //Set type of SA

#ifdef AEROELASTIC
  sndOptpar(ivar,-1);  
#endif

  switch (structdom->probType()) {

    case SolverInfo::Static:
      staticsolver->optSolveDuds();
      break;
    case SolverInfo::Modal:
      eigensolver->optSolveDwds();
      break;
    case SolverInfo::Dynamic:
       dynsolver->optDstateds(this,anagrdType);
      break;
    default:
      fprintf(stderr,
         "Error: direct analytical sensitivity analysis is implemented\n");
      fprintf(stderr,
         "       only of static, eigen and dynamic single domain problems\n");
      fprintf(stderr,
         "       ... stop\n");
      exit(-1);
  }

#ifdef AEROELASTIC
  structdom->savegradAeroForce();  
#endif
   
  gradcrit();
}
//------------------------------------------------------------------------------

void Structopt::gradcrit( ) {

 int i;
 for (i=0;i<numcrit;i++) { 

    opc[i]->gradcrit(this);
 } 
}
//------------------------------------------------------------------------------

void Structopt::gradadjoint(int icrit, double** gc) {

  anagrdType = 2;                       //Set type of SA

#ifdef AEROELASTIC
  sndOptpar(-1,icrit);  
#endif

  Vector * grad;

  zerograd();                           //Set design velocity = 0

  opc[icrit]->graddu(this);             //Build dq:ds

  switch (structdom->probType()) {      //Solve adjoint problem

    case SolverInfo::Static:            //Determine K(-1) * dq:du
      staticsolver->optAdjSolve();      
      grad = staticsolver->getpgrad();
      break;

    case SolverInfo::Dynamic:           //Determine K(-1) * dq:du
      dynsolver->optDstateds(this,anagrdType);  
      grad = dynsolver->getpgDis();
      break;

   default:
      fprintf(stderr,
         "Error: adjoint analytical sensitivity analysis is \n");
      fprintf(stderr,
         "       implemented only of static problems .. stop\n");
      exit(-1);
  }
    

  int i;
  int numabs = optsol->numvar; 
    
  for (i=0;i<numabs;i++) {              //Loop over all opt.var.

#ifndef AEROELASTIC     
     fprintf(stderr,
             " ... Analytical SA(A1): %d. optimization variable \r",i+1); 
#endif

     optpro->optvar->updgrad(i,1.0);    //Set Variation of Abs.&Struct. Var

     structdom->setoptInf(i);           //Set Opt.Var.- Ele. Influence
    
     grad->zero();                      //Initialize du:ds = 0
	
     opc[icrit]->gradcrit(this);        //Determine partial(q)/partial(s)
       
     double grad1 = opc[icrit]->grad;

     double grad2;

#ifdef AEROELASTIC
     sndOptpar(i,icrit);  
#endif
     
     switch (structdom->probType()) {   //Determine (adjoint(T) * dK:ds * u) 
       case SolverInfo::Static:
	 grad2 = staticsolver->optGetAdjointPseudo(); 
       break;
       case SolverInfo::Dynamic:
 	 grad2 = dynsolver->optGetAdjointPseudo(); 
       break;
     }  

#ifdef AEROELASTIC
    // sofar the gradients of the criteria is always saved on dir=0
    structdom->savegradAeroForce();  

    double grad3 = structdom->getgradAeroForce(0);
        
    grad2 += -grad3;
#endif

     gc[icrit][i] = grad1 + grad2;      //Store dq:ds
  }
}

//------------------------------------------------------------------------------

void Structopt::zerograd() {

    structdom->zeroGrad();
}

//------------------------------------------------------------------------------

double Structopt::initProcess(int & itrans, int & tIndex, 
                              double & stateTime, double & dt)
{
     tIndex = stateTime/dt;

     if (stateTime < dt) 
       itrans = procSetvar(-1.0,dt);
       
     return stateTime;        
}

//------------------------------------------------------------------------------

void Structopt::procEvaluate(double time, double step ) {

  FILE * optout=structdom->getOptUnit();

  double tmax = time + 0.5 * step;
  double tmin = time - 0.5 * step;
  
  int i;
  for (i=0;i<numcrit;i++) { 

     opc[i]->evaluate(this, tmax, tmin);
  }
}

//------------------------------------------------------------------------------

int Structopt::procSetvar(double time, double step ) {

  if ( ! analysisData->transFlag ) return 1;

  FILE * optout=structdom->getOptUnit();
 
  double limitdsgvel = analysisData->designVel;
  double limitaervel = analysisData->forceVel;
  
  int steplimit      = analysisData->transMax;
  int minstep        = analysisData->transMin;

  double twopi  = 6.2831853;
  double pihalf = 1.5707963;
  
  int numstep = time/step;

  int numabs = optpro->optvar->numabs;

  int iabs;
  
  if (numstep < 0) {  
    for (iabs=0;iabs<numabs;iabs++) 
      optpro->optvar->absvar[iabs]->saveval();
    maxstep=0;
    return 0;
  }

  if ( ! maxstep ) {
  
    for (iabs=0;iabs<numabs;iabs++) {  

      optpro->optvar->absvar[iabs]->setnewval();
    
      double newvar   = optpro->optvar->absvar[iabs]->newval;
      double actvar   = optpro->optvar->absvar[iabs]->val;
      double delvar   = abs(newvar - actvar);

      if ( delvar == 0.0 ) continue;

      double maxdsgvel = optpro->optvar->maxDesignVelocity(iabs);
      double maxaervel = maxSensAeroforces(iabs);     

      // design change by non-shape variation

      double one = 1.0;
      maxdsgvel = myMax(one,maxdsgvel);

      int maxdsg = ( delvar * maxdsgvel ) / ( limitdsgvel * step );
      int maxaer = ( delvar * maxaervel ) / ( limitaervel * step );
      
      fprintf(optout,"For variable : %d \n"   ,iabs+1);
      fprintf(optout,"    delvar   : %15.5e\n",delvar);
      fprintf(optout,"    maxdsgvel: %15.5e\n",maxdsgvel);
      fprintf(optout,"    maxdsg   : %d\n\n"  ,maxdsg);
      
      fprintf(optout,"For variable : %d\n"    ,iabs+1);
      fprintf(optout,"    delvar   : %15.5e\n",delvar);
      fprintf(optout,"    maxaervel: %15.5e\n",maxaervel);
      fprintf(optout,"    maxaer   : %d\n\n"  ,maxaer);
      
      maxstep=myMax(maxstep,maxdsg);
      maxstep=myMax(maxstep,maxaer);
      maxstep=myMax(maxstep,minstep);
    }
    
    if (maxstep == 0) return 1;
    
    maxstep=myMin(steplimit,maxstep);

    //    fprintf(stderr,"\n ... Transition in %d steps\n",maxstep);
    fprintf(optout,"\n ... Transition in %d steps\n",maxstep);
    
    for (iabs=0;iabs<numabs;iabs++) {  
    
      double newvar   = optpro->optvar->absvar[iabs]->newval;
      double actvar   = optpro->optvar->absvar[iabs]->val;
      double delvar   = abs(newvar - actvar);
    
      if ( delvar == 0.0 ) continue;
    
      double sum    = 0.0;
      int is;
      
      for (is=0;is<maxstep;is++) 
        sum += 1.0 + sin(twopi*(is+1.0)/maxstep - pihalf);
        
      optpro->optvar->absvar[iabs]->dynscl = (newvar - actvar) / sum;
    }
  }

  int retval=1;

  //  fprintf(stderr," ... Transition Step: %d\n",numstep);
  fprintf(optout," ... Transition Step: %d\n",numstep);

  for (iabs=0;iabs<numabs;iabs++) {
  
    double newvar = optpro->optvar->absvar[iabs]->newval;
    double actvar = optpro->optvar->absvar[iabs]->val;
    double dynscl = optpro->optvar->absvar[iabs]->dynscl;
        
    actvar += dynscl *(1.0+sin(twopi*(numstep+1.0)/maxstep - pihalf));
    
    double reach = (newvar-actvar)*dynscl;
    
    if ( reach > 0.0 ) 
      retval=0;
    else  
      actvar=newvar;
      
    optpro->optvar->absvar[iabs]->setval(actvar);

    fprintf(stderr," ... Absvar Nr. %d   Actvar: %15.5e    Newvar: %15.5e\n",
                   iabs+1,actvar,newvar);
    fprintf(optout," ... Absvar Nr. %d   Actvar: %15.5e    Newvar: %15.5e\n",
                   iabs+1,actvar,newvar);
  }     

  optpro->optvar->updvar();
  
  return retval;
}  

//------------------------------------------------------------------------------

double Structopt::maxSensAeroforces(int ivar) {

  FILE * optout=structdom->getOptUnit();

  double sensAeroforce = aerosens[ivar][0]*aerosens[ivar][0]
                       + aerosens[ivar][1]*aerosens[ivar][1]
		       + aerosens[ivar][2]*aerosens[ivar][2];
  

  if ( ! sensAeroforce ) {
    fprintf(optout,"Warning: maximum of Aeroforces is zero\n");
    return 0.0;
  }
  else {  
    double aForce = aeroforce[0]*aeroforce[0]
                  + aeroforce[1]*aeroforce[1]
		  + aeroforce[2]*aeroforce[2];
    if ( aForce ) {	
      double relForce = sqrt(sensAeroforce/aForce);	  

      fprintf(optout,"For variable: %d  aerosens[0] = %15.5e\n",
              ivar+1,aerosens[ivar][0]);
      fprintf(optout,"For variable: %d  aerosens[1] = %15.5e\n",
              ivar+1,aerosens[ivar][1]);
      fprintf(optout,"For variable: %d  aerosens[2] = %15.5e\n",
              ivar+1,aerosens[ivar][2]);
      fprintf(optout,"absolut value aeroforce = %15.5e\n\n",sqrt(aForce));
      fprintf(optout,"relative gradient = %15.5e\n\n",relForce);
      return relForce;
    }
    else {
      fprintf(optout,"Warning: Aeroforces are zero\n");
      return 0.0;      
    }
  }
}

//------------------------------------------------------------------------------
#ifdef AEROELASTIC	  

void Structopt::sndOptpar(int ivar, int icrit)
{
  double optParam[13];
  
  if(ivar < 0 & icrit < 0) anagrdType = 0;
  
  optParam[0] = (double)(optsol->numvar); 
  optParam[1] = (double)numcrit;
  optParam[2] = (double)anagrdType;

  switch (anagrdType)
  {
    case 0:
      optParam[3]  = 0.0;
      optParam[4]  = 0.0;
      optParam[5]  = 0.0;
      optParam[6]  = 0.0;
      optParam[7]  = 0.0;
      optParam[8]  = 0.0;
      optParam[9]  = 0.0;
      optParam[10] = 0.0;
      optParam[11] = 0.0;
      optParam[12] = 0.0;

      setFluidvariables(&(optParam[4]),ivar);
      break;

    case 1:
      optParam[3]  = (double)(ivar+1);
      optParam[4]  = 0.0;
      optParam[5]  = 0.0;
      optParam[6]  = 0.0;
      optParam[7]  = 0.0;
      optParam[8]  = 0.0;
      optParam[9]  = 0.0;
      optParam[10] = 0.0;
      optParam[11] = 0.0;
      optParam[12] = 0.0;

      setFluidvariables(&(optParam[4]),ivar);
      break;

    case 2:
      optParam[3]  = (double)(icrit+1);
      optParam[4]  = 0.0;
      optParam[5]  = 0.0;
      optParam[6]  = 0.0;
      optParam[7]  = 0.0;
      optParam[8]  = 0.0;
      optParam[9]  = 0.0;
      optParam[10] = 0.0; 
      optParam[11] = 0.0;
      optParam[12] = 0.0;
    
      setFluidvariables(&(optParam[4]),ivar);

      setFluidcriteria(icrit,&(optParam[10]));
      break;
  }

  structdom->sndOptpar((double*) optParam);
}

//------------------------------------------------------------------------------

void Structopt::setFluidvariables(double* flgList,int varFlg)
{
  double *     fluidattr = structdom->getFluidAttr();

  flgList[0] = fluidattr[0];
  flgList[1] = fluidattr[1];
  flgList[2] = fluidattr[2];
  

  if (varFlg >= 0) {
    double * gradfluidattr = structdom->getgradFluidAttr();
    flgList[3] = gradfluidattr[0];
    flgList[4] = gradfluidattr[1];
    flgList[5] = gradfluidattr[2];
  }
}  
 
//------------------------------------------------------------------------------

void Structopt::setFluidcriteria(int icrit,double* flgList)
{
  double dummy;

  int dir = opc[icrit]->getAeroforce(dummy);

  if(dir >= 0) flgList[dir]=1.0;
}  
#endif  
