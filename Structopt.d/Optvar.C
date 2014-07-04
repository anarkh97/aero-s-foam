#ifdef STRUCTOPT

#include <cstdlib>

#include <Structopt.d/Optvar.h>
#include <Structopt.d/Optcrit.h>
#include <Structopt.d/Optobj.h>
#include <Structopt.d/Optcon.h>
#include <Structopt.d/Optopr.h>

#include <Structopt.d/Structopt_sd.h>

#define RAND_UNIFORME static_cast<double>(random())/RAND_MAX

//------------------------------------------------------------------------------
//                            Objective 
//                                          created  9/1/98 by Kurt
//------------------------------------------------------------------------------

Stcvar * Stcvar::build ( int _num, int _typ, int _loc1, int _loc2, int _rflag) 
{
    Stcvar * stcvar = new Stcvar;

    stcvar->num   = _num;
    stcvar->typ   = _typ;
    stcvar->loc1  = _loc1;
    stcvar->loc2  = _loc2;
    stcvar->rflag = _rflag;

    stcvar->orgval = 0.0;
    
    return stcvar;
}

//------------------------------------------------------------------------------

char * Stcvar::gettypname() {

   char *typnam;

   switch (typ) {
     case 0:
       printf("Stcvar:: No type not allowed\n\n");
       exit(-1);
       break;
     case 1:
       typnam="Element Attribute";
       break;
     case 2:	 
       typnam="Nodal Coordinate";
       break;
     case 3:	 
       typnam="Element Attribute - Composite";
       break;
     case 4:	 
       typnam="Nodal Force";
       break;
     case 5:	 
       typnam="Fluid Attribute";
       break;
     case 6:	 
       typnam="Electrostatic Attribute/Variable";
       break;
     case 7:	 
       typnam="Thermal Attribute";
       break;
     default:
       printf("Stcvar:: Typ %d not implemented\n\n",typ);
       exit(-1);
   }
   
   return typnam;
}

//------------------------------------------------------------------------------

void Stcvar::setvalptr(Structopt *structopt) {
     switch (typ) {
     
        case 1:
	   pval=structopt->getptreleattr(loc1,loc2);
	   break;
        case 2:
	   pval=structopt->getptrnodcord(loc1,loc2);
	   break;
        case 3:
	   pval=structopt->getptrcomposite(loc1,loc2);
	   break;
        case 4:
	   pval=structopt->getptrnodalforce(loc1,loc2);
	   break;
        case 5:
	   pval=structopt->getptrfluidattr(loc1);
	   break;
        case 6:           
	   pval=structopt->getptrelecattr(num,loc1,loc2);           
	   break;
        case 7:
	   pval=structopt->getptrthermattr(num,loc1,loc2);
	   break;
	default:
          printf("Stcvar::setvalptr - Typ %d not implemented\n\n",typ);
	  exit(-1);
     }
}

//------------------------------------------------------------------------------

void Stcvar::setgradptr(Structopt *structopt) {
     switch (typ) {
     
        case 1:
	   pgrad=structopt->getptrgradeleattr(loc1,loc2);
	   break;
        case 2:
	   pgrad=structopt->getptrgradnodcord(loc1,loc2);
	   break;
        case 3:
	   pgrad=structopt->getptrgradcomposite(loc1,loc2);
	   break;
        case 4:
	   pgrad=structopt->getptrgradnodalforce(loc1,loc2);
	   break;
        case 5:
	   pgrad=structopt->getptrgradfluidattr(loc1);
	   break;
        case 6:
	   pgrad=structopt->getptrgradelecattr(num,loc2);
	   break;
        case 7:
	   pgrad=structopt->getptrgradthermattr(loc1,loc2);
	   break;
	default:
          printf("Stcvar::setgradptr - Typ %d not implemented\n\n",typ);
	  exit(-1);
     }
}

//------------------------------------------------------------------------------

void Stcvar::updval()
{ 
   if (!pval) {
     fprintf(stderr,"ERROR: STCVAR pointer not set\n");
     exit(-1);
   }

   if (rflag) {
     //fprintf(stderr,"updval: set orgval: %e -> %e\n",*pval,orgval);
     //fprintf(stderr,"	     new inc:	 %e\n", val);
     //fprintf(stderr,"	     new val:	 %e\n", (*pval)+val);
     orgval  =  *pval; 
     *pval  +=  val ;
   } 
   else {
     //fprintf(stderr,"updval: set val to %e\n",val);
     //fprintf(stderr,"updval: set orgval to %e\n",val);
     orgval  =  val; 
     *pval   =  val ; 
   }      
}

//------------------------------------------------------------------------------

void Stcvar::updgrad()
{ 
  if (!pgrad) {
     fprintf(stderr,"ERROR: STCVAR gradient pointer not set\n");
     exit(-1);
   }

   *pgrad = grad ;
} 

//------------------------------------------------------------------------------


void Stcvar::resetval()
{ 
   if (!pval) {
     fprintf(stderr,"ERROR: STCVAR pointer not set\n");
     exit(-1);
   }

   if (rflag) { 
     //fprintf(stderr,"reset pval to orgval: %e -> %e\n",*pval,orgval);
     *pval  =  orgval; 
   } 
   else {
     //fprintf(stderr,"reset pval: %e\n",val);
     *pval   =  val;
   }
}

//------------------------------------------------------------------------------

void Stcvar::print(FILE * optunitout) {

     char * typname=gettypname();

     fprintf(optunitout,"\n\tType of Structural Variable    : %s\n",typname);
     fprintf(optunitout,  "\tNumber of Attribute Set / Node : %d\n",loc1+1) ;
     fprintf(optunitout,  "\tType of Attribute / Coordinate : %d\n",loc2+1) ;
}     

//------------------------------------------------------------------------------

int Absvar::getDistribution() 
{
  fprintf(stderr,"WARNING: non specified virtual function ABSVAR getDistribution called\n");
  return -1;
}

double Absvar::getMeanValue()
{ 
  fprintf(stderr,"WARNING: non specified virtual function ABSVAR getMeanValue called\n");
  return -1;
} 
    
double Absvar::getStdDev()       
{ 
  fprintf(stderr,"WARNING: non specified virtual function ABSVAR getStdDev called\n");
  return -1;
}

void Absvar::setMeanValue(double v)
{
  fprintf(stderr,"WARNING: non specified virtual function ABSVAR setMean called\n");
}

void Absvar::resetMeanValue()
{
  fprintf(stderr,"WARNING: non specified virtual function ABSVAR resetMean called\n");
}

void Absvar::generateRandom()
{
  fprintf(stderr,"WARNING: non specified virtual function ABSVAR generateRandom called\n");
}

void Absvar::setFormFlag(int flag)
{
  fprintf(stderr,"WARNING: non specified virtual function ABSVAR setFormFlag called\n");
}

int Absvar::getFormFlag()
{
  fprintf(stderr,"WARNING: non specified virtual function ABSVAR getFormFlag called\n");
  return 0;
}

//------------------------------------------------------------------------------

Absvar* Absvar::create(int vartyp, absvardata & var)
{
   Absvar* absvar;

   switch (vartyp)
   {
     case 0:
       absvar = new Dsgvar(var);
       break;
     case 1:
       absvar = new Rndvar(var);
       break;
     default:
       fprintf(stderr,"ERROR: wrong type of abstract variable\n");
       exit(-1);
   }

   return absvar;
}


//------------------------------------------------------------------------------

Dsgvar::Dsgvar(absvardata & var) 
{
    typ  = 0; 

    num  = var.num;
    val  = var.val;
    scal = var.scl;
    low  = var.low;
    upp  = var.upp;

    newval = var.vin;
}


//------------------------------------------------------------------------------

void Dsgvar::print(FILE * optunitout) {

    fprintf(optunitout,"\t  %d %12.5e %12.5e %12.5e %12.5e  %12.5e\n",
	    num,val,scal,low,upp,newval);
}   

//------------------------------------------------------------------------------

Rndvar::Rndvar (absvardata & var) 
{
   typ  = 1; 

   formFlag=0;

   num  = var.num;

   dist = var.dist;

   val  = var.mean;
   scal = 1.0;
 
   mean = var.mean;
   sdev = var.sdev;
   low  = var.low;
   upp  = var.upp;

   newval = var.mean;
}

//------------------------------------------------------------------------------

void Rndvar::print(FILE * optunitout) {

    fprintf(optunitout,"\t  %d         %d  %12.5e  %12.5e %12.5e %12.5e\n",
	    num,dist,mean,sdev,low,upp);
}   

//------------------------------------------------------------------------------

double Rndvar::getScaledValue()
{
   if (!formFlag) return val*scal;

   double rc,sb,mb;

   switch (dist)
   {
      case normal:
         rc=(val-mean)/sdev; 
         break;
      case uniform:
         rc=sqrt(12.0)*(val-(upp+low)/2.0)/(upp-low);
         break;
      case lognormal:
         sb=sqrt(log(1+(sdev/mean)*(sdev/mean)));
         mb=log(mean)-0.5*sb*sb;
         rc=(log(val)-mb)/sb;
         break;
      default:
         fprintf(stderr,"RNDVAR: unspecified distribution not implemented\n");
         exit(-1);
   }

   return rc;
}

//------------------------------------------------------------------------------

double Rndvar::getScaledGrad()
{
   if (!formFlag) return grad*scal;

   double rc,sb,mb;

   switch (dist)
   {
      case normal:
         rc=grad/sdev; 
         break;
      case uniform:
         rc=sqrt(12.0)*grad/(upp-low);
         break;
      case lognormal:
         sb=sqrt(log(1+(sdev/mean)*(sdev/mean)));
         mb=log(mean)-0.5*sb*sb;
         rc=(log(grad)-mb)/sb;
      default:
         fprintf(stderr,"RNDVAR: unspecified distribution not implemented\n");
         exit(-1);
    }

   return rc;
}

//------------------------------------------------------------------------------

double Rndvar::getScaledUpperBound()
{
    if (!formFlag) {
     return upp*scal;
   }

   double rc,sb,mb;

   switch (dist)
   {
      case normal:
        rc=(upp-mean)/sdev;
        break;
      case uniform:
         rc=0.5*sqrt(12.0);
         break;
      case lognormal:
         sb=sqrt(log(1+(sdev/mean)*(sdev/mean)));
         mb=log(mean)-0.5*sb*sb;
         rc=(log(upp)-mb)/sb;
         break;
      default:
         fprintf(stderr,"RNDVAR: unspecified distribution not implemented\n");
         exit(-1);
    }

   return rc;
}

//------------------------------------------------------------------------------

double Rndvar::getScaledLowerBound()
{
     if (!formFlag) {
     return low*scal;
   }

   double rc,sb,mb;

   switch (dist)
   {
      case normal:
        rc=(low-mean)/sdev; 
        break;
      case uniform:
        rc=-0.5*sqrt(12.0);
        break;
      case lognormal:
        sb=sqrt(log(1+(sdev/mean)*(sdev/mean)));
        mb=log(mean)-0.5*sb*sb;
        rc=(log(low)-mb)/sb;
        break;
      default:
        fprintf(stderr,"RNDVAR: unspecified distribution not implemented\n");
        exit(-1);
   }

   return rc;
}

//------------------------------------------------------------------------------

void Rndvar::putScaledValue(double& u)
{
   if (!formFlag) {
     val = u/scal;
     return;
   }

   double sb,mb;

   switch (dist)
   {
      case normal:
         val = mean + sdev * u;
         break;
      case uniform:
         mb = (upp+low)/2.0;
         sb = (upp-low)/sqrt(12.0);
         val = mb + sb * u;
         break;
      case lognormal:
         sb=sqrt(log(1+(sdev/mean)*(sdev/mean)));
         mb=log(mean)-0.5*sb*sb;
         val = exp(mb + sb*u);
         break;
      default:
         fprintf(stderr,"RNDVAR: unspecified distribution not implemented\n");
         exit(-1);
    }
}

//------------------------------------------------------------------------------

void Rndvar::setgrad(double& g)
{
   if (!formFlag) {
     grad = g/scal;
     return;
   }

   double sb;

   switch (dist)
   {
      case normal:
         grad = sdev * g;
         break;
      case uniform:
         sb = (upp-low)/sqrt(12.0);
         grad = sb * g;
         break;
      case lognormal:
         sb=sqrt(log(1+(sdev/mean)*(sdev/mean)));
         grad = sb*val*g;
         break;
      default:
         fprintf(stderr,"RNDVAR: unspecified distribution not implemented\n");
         exit(-1);
    }

}

//------------------------------------------------------------------------------

void Rndvar::generateRandom()
{
    double u1,u2,z;

    switch(dist)
    {
       case Rndvar::normal:
         u1=RAND_UNIFORME;
         u2=RAND_UNIFORME;
         z=sqrt(-2.0*log(u1))*cos(2*M_PI*u2);
         val = mean+sdev*z; 
         break;  
       case Rndvar::uniform:
         val = low+RAND_UNIFORME*(upp-low);
         break;  
       case Rndvar::lognormal:
         u1=RAND_UNIFORME;
         u2=RAND_UNIFORME;
         z=sqrt(-2.0*log(u1))*cos(2*M_PI*u2);
         // Need to figure out lognormal - MPR
         val = exp(mean+sdev*z)-1; 
         //val = exp(mean+sdev*z); 
         break;  
      default:
         fprintf(stderr,"RNDVAR: unspecified distribution not implemented\n");
         exit(-1);
    }
}

//------------------------------------------------------------------------------

void Rndvar::setFormFlag(int flag)
{
  formFlag=flag;
}
   
//------------------------------------------------------------------------------

Optvar::Optvar():absvar(0),stcvar(0),opr(0) {
		 
      ivarOld=-1;  

      asTab=0;
}

//------------------------------------------------------------------------------

Optvar::~Optvar() {
		 
      if (asTab) {
        int i;
        for(i=0;i<nAbs();i++)
          if (asTab[i]) delete asTab[i];
        delete [] asTab;
      }   
}

//------------------------------------------------------------------------------

void Optvar::print(FILE * optunitout) {

     int i;

     switch (getAbsVar(0)->gettype())
     {
        case 0:
	
          fprintf(optunitout,"\n\nOptimization Variables:\n");
          fprintf(optunitout,"=======================\n\n");

          fprintf(optunitout,"Number of Optimizatin Variables: %d\n",nAbs());
          fprintf(optunitout,"Number of Structural Variables:  %d\n",nStc());
      
          fprintf(optunitout,"\nOptimization Variables:\n\n");

          fprintf(optunitout,
            "\tNumber    Value       Scaling       Lower Bd    Upper Bd.    Initial\n\n");
       
          for (i=0;i<nAbs();i++) getAbsVar(i)->print(optunitout);
          break;
     
        case 1:

          fprintf(optunitout,"\n\nRandom Variables:\n");
          fprintf(optunitout,"=======================\n\n");

          fprintf(optunitout,"Number of Random Variables:      %d\n",nAbs());
          fprintf(optunitout,"Number of Structural Variables:  %d\n",nStc());
      
          fprintf(optunitout,"\nRandom Variables:\n\n");

          fprintf(optunitout,
            "\tNumber    Type     Mean     StdDeviation    Lower Bd    Upper Bd.\n\n");
       
          for (i=0;i<nAbs();i++) getAbsVar(i)->print(optunitout);
         
	  break;
     }


     if(nStc()) {
      
       int i;
       
       fprintf(optunitout,"\nStructural Variables:\n");

       for (i=0;i<nStc();i++) {

          fprintf(optunitout,"\n\t%d. Structural Variable:\n",i+1);

          stcvar[i]->print(optunitout);

	  opr[i]->print(optunitout);
       }
     }
}

//------------------------------------------------------------------------------

void Optvar::printres(FILE * optunitout) {

     if (nAbs()) {

       fprintf(optunitout,"\n\tAbstract Variables:\n\n");
     
       int i;
       for (i=0;i<nAbs();i++)
         fprintf(optunitout,"\tVariable %d  %12.5e\n",i+1,getAbsVar(i)->getval());

     }

     if (nStc()) {
 
       fprintf(optunitout,"\n\tStructural Variables:\n\n");
     
       int i;
       for (i=0;i<nStc();i++)
         fprintf(optunitout,"\tVariable %d  %12.5e\n",i+1,stcvar[i]->getval());

     }
}

//------------------------------------------------------------------------------

void Optvar::updvar() {

      int i;
      
      for (i=0;i<nStc();i++) {
	stcvar[i]->setval(opr[i]->getval());
	stcvar[i]->updval();
      }
}

//------------------------------------------------------------------------------

void Optvar::resetvar() {

      int i;
      
      for (i=0;i<nStc();i++) {
	stcvar[i]->setval(opr[i]->getval());
	stcvar[i]->resetval();
      }
}

//------------------------------------------------------------------------------

void Optvar::initgrad() {

      if( ivarOld<0 ) {
        fprintf(stderr,"ERROR: abs-stc variable table not correctly built.\n");
        exit(-1);
      }

      int i;
      double zero=0;

      for (i=0;i<nAbs();i++) 
	{
	  getAbsVar(i)->setgrad(zero);
	}
      
      for (i=0;i<nStc();i++) 
	{
	  stcvar[i]->setgrad(opr[i]->getgrad());
	  stcvar[i]->updgrad();
	}
}

//------------------------------------------------------------------------------

void Optvar::updgrad(int ivar, double gradval) {

     double zero=0;

     getAbsVar(ivarOld)->setgrad(zero);

      getAbsVar(ivar)->setgrad(gradval);
      
      int istc,actstc,numas;

      if (!asTab[ivarOld]) {
        fprintf(stderr,
          "Abstract Optimization Variable %d not used. Stop\n",ivarOld+1);
        exit(-1);
      }

      numas = asTab[ivarOld]->size();

      for (istc=0;istc<numas;istc++) {
        actstc = asTab[ivarOld]->getEntry(istc);          
        stcvar[actstc]->setgrad(opr[actstc]->getgrad());
        stcvar[actstc]->updgrad();
      }

      if (!asTab[ivar]) {
        fprintf(stderr,
          "Abstract Optimization Variable %d not used. Stop\n",ivar+1);
        exit(-1);
      }

      numas = asTab[ivar]->size();

      for (istc=0;istc<numas;istc++) {
        actstc = asTab[ivar]->getEntry(istc);          
        stcvar[actstc]->setgrad(opr[actstc]->getgrad());
        stcvar[actstc]->updgrad();
      }

      ivarOld=ivar;
}

//------------------------------------------------------------------------------

double Optvar::maxDesignVelocity(int iabs) {

   int i;      
    
   double vx=0.0;
   double vy=0.0;
   double vz=0.0;

   double zero=0.0;
   double one=1.0;

   for (i=0;i<nAbs();i++) 
     getAbsVar(i)->setgrad(zero);

   getAbsVar(iabs)->setgrad(one);
      
   for (i=0;i<nStc();i++) { 
           
      int idir = stcvar[i]->getDesignVelocityComp();
	
      if (idir < 0) continue;
	 
      double grad = opr[i]->getgrad();
 	     grad = sqrt(grad*grad);
	 
      switch (idir) {
	 
	case 0:
	  vx=max(vx,grad);
	  break;
	case 1:
	  vy=max(vy,grad);
	  break;
	case 2:
	  vz=max(vz,grad);
      }      
   }
 
   return sqrt(vx*vx+vy*vy+vz*vz);
}

//------------------------------------------------------------------------------

void Optvar::setMeanValue(double* vals)
{
   int ivar;

   if (vals) {
     for (ivar=0;ivar<nAbs();ivar++)
       getAbsVar(ivar)->setMeanValue(vals[ivar]);
   } 
   else {
     for (ivar=0;ivar<nAbs();ivar++)
       getAbsVar(ivar)->resetMeanValue();
   } 

}

//------------------------------------------------------------------------------

void Optvar::buildAStable(int stcId)
{
   int i;

   //allocate global array

   if ( ivarOld < 0 ) {
     asTab  = new OptActInfo*[nAbs()];
     ivarOld = 0;
     for(i=0;i<nAbs();i++)
       asTab[i]=0;
   }      

   //extract abstract varibles used in func structure
   
   OptActInfo stcInfo;
   opr[stcId]->extractAbsvar(stcInfo);

   int num = stcInfo.size();

   if ( num ) {

     for(i=0;i<num;i++) {

       int iabs=stcInfo.getEntry(i);

       if (! asTab[iabs] ) asTab[iabs] = new OptActInfo;
       asTab[iabs]->add(stcId);
     }
   }
}

//------------------------------------------------------------------------------

void Optvar::generateRandom()
{
  int i;

  for (i=0;i<nAbs();i++)
    getAbsVar(i)->generateRandom();
}

//------------------------------------------------------------------------------

void Optvar::setFormFlag(int flag)
{
  int i;

  for (i=0;i<nAbs();i++)
    getAbsVar(i)->setFormFlag(flag);
}

//------------------------------------------------------------------------------

double Optvar::getFormL2norm2()
{
  double l2norm2=0.0;

  int i;

  for (i=0;i<nAbs();i++) {
    int formFlag = getAbsVar(i)->getFormFlag();
    if ( formFlag ) {
      double val=getAbsVar(i)->getScaledValue();
      l2norm2+=val*val;
    } else
      fprintf(stderr,"WARNING: forml2norm required for non-form type variable\n");
  }

  return l2norm2;

}

//------------------------------------------------------------------------------

double Optvar::getgradFormL2norm2()
{
  double gradl2norm2=0.0;

  int i;

  for (i=0;i<nAbs();i++) {
    int formFlag = getAbsVar(i)->getFormFlag();
    if ( formFlag ) {
      double val=getAbsVar(i)->getScaledValue();
      double grad=getAbsVar(i)->getScaledGrad();
      gradl2norm2+=2.0*val*grad;
    } else
      fprintf(stderr,"WARNING: forml2norm required for non-form type variable\n");
  }

  return gradl2norm2;
}

//------------------------------------------------------------------------------

OptActInfo::OptActInfo():elements(0,5) 
{ 
  lsize = 0; 
  typ=Optvar::notyp; 
}

//------------------------------------------------------------------------------

void OptActInfo::add( int ele, int _typ ) 
{ 

   if ( typ != Optvar::notyp && typ != _typ && _typ != Optvar::mpattr) 
     typ = Optvar::mixedvar;
   else
     typ = _typ;

   if (ele < 0 )  return;
   
   elements[lsize]=ele; 
   lsize++; 
}

//------------------------------------------------------------------------------

int OptActInfo::getEntry( int id ) 
{ 
   if (id > lsize-1) {
     fprintf(stderr,"ERROR in OptActInfo::getEntry\n");
     exit(-1);
   }
   
   return elements[id];
}

//------------------------------------------------------------------------------

int OptActInfo::check( int num ) 
{ 
   int i;
   for (i=0;i<lsize;i++)
      if (num == elements[i])
        return 1;
   
   return 0;
}                             

//------------------------------------------------------------------------------

char* OptActInfo::typname() 
{  
   char *typnam;

   switch (typ) {
   
     case Optvar::notyp:   
       typnam="Notyp";
       break;
     case Optvar::attribute:   
       typnam="ATT";
       break;
     case Optvar::coordinate:   
       typnam="COO";
       break;
     case Optvar::composite:   
       typnam="COM";
       break;
     case Optvar::nodalforce:   
       typnam="NODF";
       break;
     case Optvar::fluidparam:   
       typnam="FLUID";
       break;
     case Optvar::elecattr:   
       typnam="ELCATT";
       break;
     case Optvar::thermattr:   
       typnam="THERMATT";
       break;
     case Optvar::mixedvar:   
       typnam="MIX";
       break;
     default:      
       typnam="ERROR";
  }

  return typnam;
}                             

//------------------------------------------------------------------------------
void Optvar::addAbsVar(size_t varno, Absvar* varptr)
{
  if(varno < absvar.size() && absvar[varno] != 0)
    {
      fprintf(stderr,
	      "Error generating abs.var.: %d. variable already defined\n",
	      varno+1);
      exit(-1);
    }
  if(absvar.size() <= varno)
    { absvar.resize(varno+1, 0); }
  absvar[varno] = varptr;
  return;
}

//------------------------------------------------------------------------------
void Optvar::addStcVar(size_t varno, Stcvar* varptr)
{
  if(varno < stcvar.size() && stcvar[varno] != 0)
    {
      fprintf(stderr,
	      "Error generating stc.var.: %d. variable already defined\n",
	      varno+1);
      exit(-1);
    }
  if(stcvar.size() <= varno)
    { stcvar.resize(varno+1, 0); }
  stcvar[varno] = varptr;
  return;
}

//------------------------------------------------------------------------------
void Optvar::addOpr(size_t varno, Optopr* varptr)
{
  if(varno < opr.size() && opr[varno] != 0)
    {
      fprintf(stderr,
	      "Error generating opt.opr.: %d. already defined\n",
	      varno+1);
      exit(-1);
    }
  if(opr.size() <= varno)
    { opr.resize(varno+1, 0); }
  opr[varno] = varptr;
  return;
}

//------------------------------------------------------------------------------
void Optvar::setvalptr(Structopt* sopt)
{
  for(std::vector<Stcvar*>::iterator i=stcvar.begin();
      i != stcvar.end(); ++i)
    {
      (*i)->setvalptr(sopt);
    }
  return;
}


//------------------------------------------------------------------------------
void Optvar::setgradptr(Structopt* sopt)
{
  for(std::vector<Stcvar*>::iterator i=stcvar.begin();
      i != stcvar.end(); ++i)
    {
      (*i)->setgradptr(sopt);
    }
  return;
}

#endif
