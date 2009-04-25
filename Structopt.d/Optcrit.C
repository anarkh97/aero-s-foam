#ifdef STRUCTOPT

#include <stdio.h>
#include <Structopt.d/Optcrit.h>
#include <Structopt.d/Structopt_sd.h>
#include <Math.d/mathUtility.h>

class Domain;

//------------------------------------------------------------------------------
//                          Optimization Criteria 
//                      Base Class and Derived Classes
//                                                    created  9/1/98 by Kurt
//------------------------------------------------------------------------------
char * Criteriatyp::gettypname() {

   char *typnam;

   switch (typ) {
     case 0:
       typnam="StrainEnergy";
       break;
     case 1:
       typnam="Mass";
       break;
     case 2:
       typnam="Frequency";
       break;
     case 3:	 
       typnam="Nodal Stress";
       break;
     case 4:
       typnam="Displacement";
       break;
     case 5:
       typnam="Velocity";
       break;
     case 6:
       typnam="Acceleration";
       break;
     case 7:
       typnam="KineticEnergy";
       break;
     case 8:
       typnam="DissipationEnergy";
       break;
     case 9:
       typnam="AeroForce/Moments/SonicBoom";
       break;
     case 10:
       typnam="Internal Force";
       break;
     case 11:
      typnam="ControlCost";
      break;
     case 12:
      typnam="StressIntegral";
      break;
     case 13:
      typnam="Displacement Leveling";
      break;
     case 14:
      typnam="FORM L2-norm";
      break;
     case 15:
      typnam="Failure Probability";
      break;
     case 16:
      typnam="Electrostatic Criteria";
      break;
     case 17:
      typnam="Reliability Index";
      break;
     case 18:
      typnam="Performance Measure";
      break;
     case 19:
      typnam="Moment of Inertia";
      break;
     case 20:
      typnam="Stress Leveling";
      break;
     case 21:
      typnam="Nodal Temperature";
      break;
     case 22:
      typnam="Electro-mechanical pull-in instability";
      break;
     case 23:
       typnam="Change of Nodal Position";
       break;
     case 24:
       typnam="Nodal Voltage";
       break;
     case 25:
       typnam="Electrical RHS";
       break;
     case 26:
       typnam="Electrical Power";   
       break;   
     case 27:
       typnam="Thermal Power Loss from Conduction";   
       break;    
     case 28:
       typnam="Single Domain Temperature";
       break;
     case 29:
       typnam="Nonlinear Load Control Parameter";
       break;      
     case 30:
       typnam="Mass in deformed configuration";
       break;      
    default:
       fprintf(stderr,"Criteria Type not defined\n");
       exit(-1);
   }
   
   return typnam;
}

//------------------------------------------------------------------------------
char * Stresstyp::gettypname() {

   char *typnam;

   switch (typ) {
     case 0:
      typnam="von Mises Stress on Lower Surface";
       break;
     case 1:
      typnam="von Mises Stress on Middle Surface";
       break;
     case 2:
      typnam="von Mises Stress on Upper Surface";
       break;
     case 10:
       typnam="1. Principal Stress on Lower Surface";
       break;
     case 11:
       typnam="1. Principal Stress on Middle Surface";
       break;
     case 12:	 
       typnam="1. Principal Stress on Upper Surface";
       break;
     case 20:
       typnam="2. Principal Stress on Lower Surface";
       break;
     case 21:
       typnam="2. Principal Stress on Middle Surface";
       break;
     case 22:
       typnam="2. Principal Stress on Upper Surface";
       break;
     case 30:
      typnam="Normal Stress in X-Direction on Lower Surface";
       break;
     case 31:
      typnam="Normal Stress in X-Direction on Middle Surface";
       break;
     case 32:
      typnam="Normal Stress in X-Direction on Upper Surface";
       break;
     case 40:
      typnam="Normal Stress in Y-Direction on Lower Surface";
       break;
     case 41:
      typnam="Normal Stress in Y-Direction on Middle Surface";
       break;
     case 42:
      typnam="Normal Stress in Y-Direction on Upper Surface";
       break;
     case 50:
      typnam="Normal Stress in Z-Direction on Lower Surface";
       break;
     case 51:
      typnam="Normal Stress in Z-Direction on Middle Surface";
       break;
     case 52:
      typnam="Normal Stress in Z-Direction on Upper Surface";
       break;
     default:
       fprintf(stderr,"Stress Type not specified\n");
       exit(-1);
   }
   
   return typnam;
}

//------------------------------------------------------------------------------
char * Disptyp::gettypname() {

   char *typnam;

   switch (typ) {
     case 0:
       typnam="x-Direction";
       break;
     case 1:
       typnam="y-Direction";
       break;
     case 2:
       typnam="z-Direction";
       break;
     case 3:	 
       typnam="x-Rotation";
       break;
     case 4:
       typnam="y-Rotation";
       break;
     case 5:
       typnam="z-Rotation";
       break;
     case 6:
       typnam="Length of Displacement/Force Vector";
       break;
     case 7:
       typnam="Temperature";
       break;
     case 8:
       typnam="x-Direction real part";
       break;
     case 9:
       typnam="y-Direction real part";
       break;
     case 10:
       typnam="z-Direction real part";
       break;
     case 11:
       typnam="Length of real part vector";
       break;
     case 12:
       typnam="x-Direction imaginary part";
       break;
     case 13:
       typnam="y-Direction imaginary part";
       break;
     case 14:
       typnam="z-Direction imaginary part";
       break;
     case 15:
       typnam="Length of imaginary part vector";
       break;
     default:
       fprintf(stderr,"Displacement/Force Type not specified\n");
       exit(-1);
   }
   
   return typnam;
}

//------------------------------------------------------------------------------
Optcrit * Optcrit::build(int _num, critdata& crit, double& _time, int _anaId )
{
   Optcrit *optcr=0;
   
   switch ( crit.typ ) {
     case 0:
       optcr = new OptcritGlobal(0,crit.ip,crit.ipsize);
       optcr->crittyp.typ=0;
       break;
     case 1:
       optcr = new OptcritGlobal(0,crit.ip,crit.ipsize);
       optcr->crittyp.typ=1;
       break;
     case 2:
       optcr = new OptcritGlobal(crit.i[0]);
       optcr->crittyp.typ=2;
       break;
     case 3:	 
       optcr = new OptcritNodstr(crit.i[0],crit.i[1]);    
       optcr->crittyp.typ=3;
       break;
     case 4:
       optcr = new OptcritDisp(crit.i[0],crit.i[1]);
       optcr->crittyp.typ=4;
       break;
     case 5:
       optcr = new OptcritDisp(crit.i[0],crit.i[1]);
       optcr->crittyp.typ=5;
       break;
     case 6:
       optcr = new OptcritDisp(crit.i[0],crit.i[1]);
       optcr->crittyp.typ=6;
       break;
     case 7:
       optcr = new OptcritGlobal;
       optcr->crittyp.typ=7;
       break;
     case 8:
       optcr = new OptcritGlobal;
       optcr->crittyp.typ=8;
       break;
     case 9:
       optcr = new OptcritGlobal(crit.i[0]);
       optcr->crittyp.typ=9;
       break;
     case 10:
       optcr = new OptcritDisp(crit.i[0],crit.i[1]);
       optcr->crittyp.typ=10;
       break;
     case 11:
       optcr = new OptcritGlobal;
       optcr->crittyp.typ=11;
       break;
     case 12:
       optcr = new OptcritGlobal(crit.i[0],crit.d[0],crit.d[1],crit.i[1],
                                 crit.i[2],crit.ip,crit.ipsize);
       optcr->crittyp.typ=12;
       break;
     case 13:
       optcr = new OptcritDispLevel(crit.i,crit.d,crit.dl,crit.dlsize);
       optcr->crittyp.typ=13;
       break;
     case 14:
       optcr = new OptcritGlobal;
       optcr->crittyp.typ=14;
       break;
     case 15:
       optcr = new OptcritGlobal(crit.i[0]);
       optcr->crittyp.typ=15;
       break;
     case 16:
       optcr = new OptcritGlobal(crit.i[0]);
       optcr->crittyp.typ=16;
     case 17:
       optcr = new OptcritGlobal(crit.i[0]);
       optcr->crittyp.typ=17;
       break;
     case 18:
       optcr = new OptcritGlobal(crit.i[0]);
       optcr->crittyp.typ=18;
       break;
     case 19:
       optcr = new OptcritGlobal(crit.i[0]);
       optcr->crittyp.typ=19;
       break;
     case 20:
       optcr = new OptcritStressLevel(crit.i,crit.d,crit.dl,crit.dlsize);
       optcr->crittyp.typ=20;
       break;
     case 21:
       optcr = new OptcritExt(crit.i[0]);
       optcr->crittyp.typ=21;
       break;
     case 22:
       optcr = new OptcritElecMech();
       optcr->crittyp.typ=22;
       break;
     case 23:
       optcr = new OptcritDisp(crit.i[0],crit.i[1]);
       optcr->crittyp.typ=23;
       break;
     case 24:
       optcr = new OptcritExt(crit.i[0]);
       optcr->crittyp.typ=24;
       break;
     case 25:
       optcr = new OptcritExt(crit.i[0]);
       optcr->crittyp.typ=25;
       break;    
     case 26:
       optcr = new OptcritExt();
       optcr->crittyp.typ=26;
       break;  
     case 27:
       optcr = new OptcritExt();
       optcr->crittyp.typ=27;
       break;    
     case 28:
       optcr = new OptcritDisp(crit.i[0],0);
       optcr->crittyp.typ=28;
       break;
     case 29:
       optcr = new OptcritDisp();
       optcr->crittyp.typ=29;
       break;
     case 30:
       optcr = new OptcritGlobal(crit.d[0],crit.ip,crit.ipsize);
       optcr->crittyp.typ=30;
       break;
     default:
       fprintf(stderr,"Building Opt.Criteria: Type not specified\n");
       exit(-1);
   }

   optcr->num   = _num;
   optcr->time  = _time;
   optcr->anaId = _anaId;

   // Initialize value and gradients
   
   optcr->val  = 0.0;
   optcr->grad = 0.0;

   return optcr;
}

//------------------------------------------------------------------------------

int Optcrit::checkTime( double tmax, double tmin)
{
   // return 1 if tmin = -1 (time is load case or of no importance)

   if (tmin == -1 ) return 1;

   // return 0 if time > tmax or time < tmin

   if ( time > tmax || time < tmin ) return 0;

   // otherwise return 1

   return 1;
}

//------------------------------------------------------------------------------

int Optcrit::getAnalysisType()
{
  switch (crittyp.typ) 
  {
    case 1: 
      return GeometricAnalysis;
      break;
    case 2: 
      return ModalAnalysis;
      break;
    case 14: 
      return GeometricAnalysis;
      break;
    case 15: 
      return GeometricAnalysis;
      break;
    case 17: 
      return GeometricAnalysis;
      break;
    case 18: 
      return GeometricAnalysis;
      break;
    case 19: 
      return GeometricAnalysis;
      break;
    case 21:
      return ExtCritAnalysis;
      break;  
    case 23: 
      return GeometricAnalysis;
      break;
    case 24:
      return ExtCritAnalysis;
      break; 
    case 25:
      return ExtCritAnalysis;
      break;     
    case 26:
      return ExtCritAnalysis;
      break;  
    case 27:
      return ExtCritAnalysis;
      break;   
    default:
      return StressAnalysis;
  }
}

//------------------------------------------------------------------------------

int Optcrit::getEigvalId()
{
  fprintf(stderr," ERROR: Freqency ID required from a non-modal criteria \n");
  exit(-1);
  return -1;
}

//------------------------------------------------------------------------------

// strain energy, mass, frequency, aeroforces, control costs

OptcritGlobal::OptcritGlobal ( int _iloc1, int* _list, int _listsize)
{
   iloc1    = _iloc1; 
   listsize = _listsize;
   list     =  0;

   if(listsize) {
  
     list = new int[listsize];
 
     int i;
     for (i=0;i<listsize;i++)
        list[i] = _list[i];
   }
}

//------------------------------------------------------------------------------

// mass in deformed configuration

OptcritGlobal::OptcritGlobal ( double& _d0, int* _list, int _listsize)
{
   powFac   = _d0; 
   listsize = _listsize;
   list     =  0;

   if(listsize) {
  
     list = new int[listsize];
 
     int i;
     for (i=0;i<listsize;i++)
        list[i] = _list[i];
   }
}

//------------------------------------------------------------------------------

// stress integral

OptcritGlobal::OptcritGlobal( int _strtyp, double _d1, double _d2, int _i1, 
                              int _i2, int* _list, int _listsize) 
{ 
   strtyp.typ =_strtyp; 

   qFlag = _i1;

   nodeleTyp = _i2;
  
   refVal=_d1; 
   powFac=_d2; 

   iloc1= 0; 
   iloc2= _strtyp; 

   listsize = _listsize;
   list     =  0;

   if(listsize) {
  
     list = new int[listsize];
 
     int i;
     for (i=0;i<listsize;i++)
        list[i] = _list[i];
   }
} 

//------------------------------------------------------------------------------

void OptcritGlobal::evaluate(Structopt *structopt, double tmax, double tmin) 
{
    // Check Time
    if ( ! checkTime(tmax,tmin) ) return; 

    int typ=crittyp.typ;
     
    switch (typ) {
     
    case 0:
       val = structopt->getstrainenergy(list,listsize,time,anaId);
       break;
    case 1:
       val = structopt->getmass(list,listsize,anaId);
       break;
    case 2:
       val = structopt->getfrequency(iloc1,anaId);
       break;
    case 7:
       val = structopt->getkineticenergy(time,anaId);
       break;
    case 8:
       val = structopt->getdampingenergy(time,anaId);
       break;
    case 9:
       val = structopt->getaeroforce(iloc1,time,anaId);
       break;
    case 11:
       val = structopt->getcontrolcost(anaId);
       break;
    case 12:
       val = structopt->getstressint(strtyp.typ,refVal,powFac,qFlag,nodeleTyp,
                                     list,listsize,time,anaId);
       break;
    case 14:
       val = structopt->getforml2norm2();
       break;
    case 15:
       val = structopt->getFailprob(0,iloc1);
       break;
    case 16:
       val = structopt->getestat(iloc1,list,listsize,time,anaId);
       break;
    case 17:
       val = structopt->getFailprob(1,iloc1);
       break;
    case 18:
       val = structopt->getFailprob(2,iloc1);
       break;
    case 19:
       val = structopt->getMomOfInertia(list,listsize,iloc1,anaId);
       break;
    case 30:
       val = structopt->getDCmass(powFac,list,listsize,time,anaId);
       break;
   default:
       fprintf(stderr,"Evaluate Global Opt.Criteria: Type not specified !\n");
       exit(-1);
   }           
}

//------------------------------------------------------------------------------

void OptcritGlobal::gradcrit(Structopt *structopt, int actvar) {

    int typ=crittyp.typ;
     
    switch (typ) {
     
    case 0:
       grad = structopt->getgradstrainenergy(list,listsize,time,anaId);
       break;
    case 1:
       grad = structopt->getgradmass(list,listsize,anaId);
       break;
    case 2:
       grad = structopt->getgradfrequency(iloc1,anaId);
       break;
    case 7:
       grad = structopt->getgradkineticenergy(time,anaId);
       break;
    case 8:
       grad = structopt->getgraddampingenergy(time,anaId);
       break;
    case 9:
       grad = structopt->getgradaeroforce(iloc1,time,anaId);
       break;
    case 12:
       grad = structopt->getgradstressint(strtyp.typ,refVal,powFac,qFlag,nodeleTyp,
                                          list,listsize,time,anaId);
       break;
    case 14:
       grad = structopt->getgradforml2norm2();
       break;
    case 15:
       grad = structopt->getgradFailprob(0,iloc1,actvar);
       break;
    case 16:
       grad = structopt->getgradestat(iloc1,list,listsize,time,anaId);
       break;
    case 17:
       grad = structopt->getgradFailprob(1,iloc1,actvar);
       break;
    case 18:
       grad = structopt->getgradFailprob(2,iloc1,actvar);
       break;
    case 19:
       grad = structopt->getgradMomOfInertia(list,listsize,iloc1,anaId);
       break;
    case 30:
       grad = structopt->getgradDCmass(powFac,list,listsize,time,anaId);
       break;
    default:
       fprintf(stderr,"Evaluate Derivative of Global Opt.Criteria: Type not specified !\n");
       exit(-1);
   }           
}

//------------------------------------------------------------------------------

void OptcritGlobal::gradpart(Structopt *structopt) {

    int typ=crittyp.typ;
     
    switch (typ) {
     
    case 0:
       grad = structopt->getgradpartstrainenergy(list,listsize,time,anaId);
       break;
    case 1:
       grad = structopt->getgradmass(list,listsize,anaId);
       break;
    case 2:
       fprintf(stderr,"Frequency: no gradpart ! \n");
       exit(-1);
    case 7:
       fprintf(stderr,"Kinetic Energy: no gradpart ! \n");
       exit(-1);
      break;
    case 8:
       fprintf(stderr,"Dissipation Energy: no gradpart ! \n");
       exit(-1);
       break;
    case 9:
       grad = structopt->getgradaeroforce(iloc1,time,anaId);
       break;
    case 12:
       grad = structopt->getgradpartstressint(strtyp.typ,refVal,powFac,qFlag,nodeleTyp,
                                              list,listsize,time,anaId);
       break;
    case 14:
       grad = structopt->getgradforml2norm2();
       break;
    case 15:
       grad = structopt->getgradFailprob(0,iloc1);
       break;
    case 16:
       grad = structopt->getgradpartestat(iloc1,list,listsize,time,anaId);
       break;
    case 17:
       grad = structopt->getgradFailprob(1,iloc1);
       break;
    case 18:
       grad = structopt->getgradFailprob(2,iloc1);
       break;
    case 19:
       grad = structopt->getgradMomOfInertia(list,listsize,iloc1,anaId);
       break;
    case 30:
       grad = structopt->getgradpartDCmass(powFac,list,listsize,time,anaId);
       break;
    default:
       fprintf(stderr,"Evaluate Derivative of Global Opt.Criteria: Type not specified !\n");
       exit(-1);
   }           
}

//------------------------------------------------------------------------------

void OptcritGlobal::graddu(Structopt *structopt) {
 
    int typ=crittyp.typ;
     
    switch (typ) {
     
    case 0:
       structopt->getgraddustrainenergy(list,listsize,time,anaId);
       break;
    case 1:
       fprintf(stderr,"Structural Mass: no graddu ! \n");
       exit(-1);
       break;
    case 2:
       fprintf(stderr,"Frequency: no graddu ! \n");
       exit(-1);
       break;
    case 7:
       fprintf(stderr,"Kinetic Energy: no graddu ! \n");
       exit(-1);
      break;
    case 8:
       fprintf(stderr,"Dissipation Energy: no graddu ! \n");
       exit(-1);
       break;
    case 9:
       structopt->getgradduaeroforce(iloc1,time,anaId);
       break;
    case 11:
       fprintf(stderr,"ControlCost: no graddu ! \n");
       exit(-1);
       break;
    case 12:
       structopt->getgraddustressint(strtyp.typ,refVal,powFac,qFlag,nodeleTyp,
                                     list,listsize,time,anaId);
       break;
    case 14:
       fprintf(stderr,"FORM L2-norm: no graddu ! \n");
       exit(-1);
       break;
    case 15:
       fprintf(stderr,"Failure Probability: no graddu ! \n");
       exit(-1);
       break;
    case 16:
       fprintf(stderr,"Electrostatic Criteria: no graddu ! \n");
       exit(-1);
    case 17:
       fprintf(stderr,"Reliability Index: no graddu ! \n");
       exit(-1);
       break;
    case 18:
       fprintf(stderr,"Performance Measure: no graddu ! \n");
       exit(-1);
       break;
    case 19:
       fprintf(stderr,"Moment of Inertia: no graddu ! \n");
       exit(-1);
       break;
    case 30:
       structopt->getgradduDCmass(powFac,list,listsize,time,anaId);
       break;
    default:
       fprintf(stderr,"Evaluate Global Opt.Criteria: Type not specified !\n");
       exit(-1);
   }           
}

//------------------------------------------------------------------------------

int OptcritGlobal::getEigvalId (){

  if ( crittyp.typ == 2 ) {
    return iloc1;
  }  
  else {
    fprintf(stderr," ERROR: Freqency ID required from a non-modal criteria \n");
    exit(-1);
  }
  
  return -1;
}

//------------------------------------------------------------------------------

int OptcritGlobal::getAeroforce ( double &aforce )  {

  int irc;

  switch(crittyp.typ) {
     case 9:
        aforce = val;
        irc    = iloc1;
        break;
     default:
        aforce = 0.0;
        irc    = -1;
  }

  return irc;
}

//------------------------------------------------------------------------------

void OptcritDisp::evaluate(Structopt *structopt, double tmax, double tmin) 
{
    // Check Time
    if ( ! checkTime(tmax,tmin) ) return; 

    enum {DISP,VEL,ACC};

    int typ=crittyp.typ;
     
    switch (typ) {

    case 4:
       val=structopt->getdisp(node,disptyp.typ,DISP,time,anaId);
       break;
    case 5:
       val=structopt->getdisp(node,disptyp.typ,VEL,time,anaId);
       break;
    case 6:
       val=structopt->getdisp(node,disptyp.typ,ACC,time,anaId);
       break;
    case 10:
       val=structopt->getinternalforce(node,disptyp.typ,time,anaId);
       break;
    case 23:
       val=structopt->getnodalpos(node,disptyp.typ,time,anaId);
       break;
    case 28:
       val=structopt->getdisp(node,7,0,time,anaId);
       break;
    case 29:
       val=structopt->getLambda(time,anaId);  
       break;    
    default:
       fprintf(stderr,"Error Evaluating Displacement|Velocity|Accelertion! \n");
       exit(-1);
   }           
}

//------------------------------------------------------------------------------

void OptcritDisp::gradcrit(Structopt *structopt, int actvar) {

    enum {DISP,VEL,ACC};

    int typ=crittyp.typ;
     
    switch (typ) {

    case 4:
       grad=structopt->getgraddisp(node,disptyp.typ,DISP,time,anaId);
       break;
    case 5:
       grad=structopt->getgraddisp(node,disptyp.typ,VEL,time,anaId);
       break;
    case 6:
       grad=structopt->getgraddisp(node,disptyp.typ,ACC,time,anaId);
       break;
    case 10:
       grad=structopt->getgradinternalforce(node,disptyp.typ,time,anaId);
       break;
    case 23:
       grad=structopt->getgradnodalpos(node,disptyp.typ,time,anaId);
       break;
    case 28:
       grad=structopt->getgraddisp(node,7,0,time,anaId);
       break;
    case 29:
       grad=structopt->getgradLambda(time,anaId);  
       break;      
    default:
       fprintf(stderr,"Error Deriving Displacement|Velocity|Accelertion|FORCE! \n");
       exit(-1);
   }           
}

//------------------------------------------------------------------------------

void OptcritDisp::gradpart(Structopt *structopt) {

    enum {DISP,VEL,ACC};

    int typ=crittyp.typ;
     
    switch (typ) {

    case 4:
    case 5:
    case 6:
    case 28:
    case 29:
       grad=0.0;
       break;
    case 10:
       grad=structopt->getgradpartinternalforce(node,disptyp.typ,time,anaId);
       break;
    case 23:
       grad=structopt->getgradnodalpos(node,disptyp.typ,time,anaId);
       break;
    default:
       fprintf(stderr,"Error Deriving Displacement|Velocity|Accelertion|FORCE! \n");
       exit(-1);
   }           
}

//------------------------------------------------------------------------------

void OptcritDisp::graddu(Structopt *structopt) {

    enum {DISP,VEL,ACC};

    int typ=crittyp.typ;

    switch (typ) {

    case 4:
       structopt->getgraddudisp(node,disptyp.typ,DISP,time,anaId);
       break;
    case 5:
       fprintf(stderr,"Velocities: no graddu ! \n");
       exit(-1);
       break;
    case 6:
       fprintf(stderr,"Accelerations: no graddu ! \n");
       exit(-1);
       break;
    case 10:
       structopt->getgradduinternalforce(node,disptyp.typ,time,anaId);
       break;
    case 23:
       fprintf(stderr,"Change of Nodal Position: no graddu ! \n");
       exit(-1);
       break;
    case 28:
       structopt->getgraddudisp(node,7,0,time,anaId);
       break;
    case 29:
       structopt->getgraddulambda(time,anaId);
       break;      
    default:
       fprintf(stderr,"Error Evaluating Displacement|Velocity|Accelertion|FORCE! \n");
       exit(-1);
   }           
}

//------------------------------------------------------------------------------

void OptcritNodstr::evaluate(Structopt *structopt, double tmax, double tmin)
{
    // Check Time
    if ( ! checkTime(tmax,tmin) ) return; 

    val=structopt->getnodstr(node,strtyp.typ,time,anaId);
}

//------------------------------------------------------------------------------

void OptcritNodstr::gradcrit(Structopt *structopt, int actvar) {

    grad=structopt->getgradnodstr(node,strtyp.typ,time,anaId);
}

//------------------------------------------------------------------------------

void OptcritNodstr::gradpart(Structopt *structopt) {

    grad=structopt->getgradpartnodstr(node,strtyp.typ,time,anaId);
}

//------------------------------------------------------------------------------

void OptcritNodstr::graddu(Structopt *structopt) {

    structopt->getgraddunodstr(node,strtyp.typ,time,anaId);
}

//------------------------------------------------------------------------------

void OptcritGlobal::print(FILE * optunitout) {
     
     int typ=crittyp.typ;
     
     char *critname=crittyp.gettypname();

     fprintf(optunitout,"\t%d. Design Criteria - Type: Global | %s \n",
             num,critname);

     fprintf(optunitout,"\tAnalysisId: %d   time: %e\n",anaId,time);
     
     if ( typ == 2 ) fprintf(optunitout,"\tEigenvalue Nr.: %d \n",iloc1+1);

     if ( typ == 9 ) {
       switch (iloc1) {
          case 0: 
          case 1: 
          case 2:
            fprintf(optunitout,"\tForce - Direction : %d \n",iloc1+1);
            break;
          case 3: 
          case 4: 
          case 5:
            fprintf(optunitout,"\tMoment - Direction : %d \n",iloc1-2);
            break;
          case 6:
            fprintf(optunitout,"\tSonic Boom :\n");
            break;
       }
     }
     
     if ( typ == 12 ) {    
       char *stressname = strtyp.gettypname();
       fprintf(optunitout,"\n");
       fprintf(optunitout,"\tStress           : %s \n",stressname);
       fprintf(optunitout,"\tReference Stress : %e \n",refVal);
       fprintf(optunitout,"\tPower Factor     : %e \n",powFac);
       fprintf(optunitout,"\tDividing by Area : %d \n",qFlag);
       fprintf(optunitout,"\tSelected elements: %d \n",listsize);
       fprintf(optunitout,"\n");
       if (nodeleTyp)
       fprintf(optunitout,"\t%d Nodal stress values used.\n",listsize);
       else
       fprintf(optunitout,"\t%d Elemental stress values used.\n",listsize);

     }
 
     if ( typ == 14 ) fprintf(optunitout,"\tCriteria Nr.: %d \n",iloc1+1);

     if ( typ == 19 ) {
       fprintf(optunitout,"\n");
       fprintf(optunitout,"\tabout axis: %d \n",iloc1+1);
     }

     if ( typ == 30 ) {
       fprintf(optunitout,"\tPower Factor     : %e \n",powFac);
     }
     
     fprintf(optunitout,"\n");
}

//------------------------------------------------------------------------------

void OptcritNodstr::print(FILE * optunitout) {

     char *critname   = crittyp.gettypname();
     char *stressname = strtyp.gettypname();

     fprintf(optunitout,"\t%d. Design Criteria  - Type: Local | %s \n",
             num,critname);
     fprintf(optunitout,"\tAnalysisId: %d   time: %e\n",anaId,time);
     fprintf(optunitout,"\tStress: %s at Node: %i\n\n",stressname,node+1);
}

//------------------------------------------------------------------------------

void OptcritDisp::print(FILE * optunitout) {

     char *critname = crittyp.gettypname();
     char *dispname = disptyp.gettypname();

     fprintf(optunitout,"\t%d. Design Criteria  - Type: Local | %s \n ",
             num,critname);
     fprintf(optunitout,"\tAnalysisId: %d   time: %e\n",anaId,time);
     fprintf(optunitout,"\tType: %s at Node: %i\n\n",
               dispname,node+1);
}


//------------------------------------------------------------------------------

// displacement leveling

OptcritDispLevel::OptcritDispLevel ( int* id, double * xd, dldata * dl, int dlsize )
{ 
   size   = dlsize;

   aFlag  = id[0];         //0: no averaging, 1: averaging, 2: reference node displ.
   qFlag  = id[1];

   powFac = xd[0];

   nodeid = new int[dlsize];
   distyp = new int[dlsize];

   refVal = new double[dlsize];
   difVal = new double[dlsize];

   int i;
   for (i=0;i<dlsize;i++) {
     nodeid[i] = dl[i].node;
     distyp[i] = dl[i].typ;
     refVal[i] = dl[i].refVal;
     difVal[i] = dl[i].difVal; 
   }
} 

//------------------------------------------------------------------------------

void OptcritDispLevel::evaluate(Structopt *structopt, double tmax, double tmin) 
{
    // Check Time
    if ( ! checkTime(tmax,tmin) ) return; 

    val = structopt->getdisplevel(aFlag,qFlag,size,nodeid,distyp,
                                  refVal,difVal,powFac,time,anaId);
}

//------------------------------------------------------------------------------

void OptcritDispLevel::gradcrit(Structopt *structopt, int actvar) 
{
    grad = structopt->getgraddisplevel(aFlag,qFlag,size,nodeid,distyp,
                                       refVal,difVal,powFac,time,anaId);
}


//------------------------------------------------------------------------------

void OptcritDispLevel::gradpart(Structopt *structopt) 
{
    grad = structopt->getgradpartdisplevel(aFlag,qFlag,size,nodeid,distyp,
                                           refVal,difVal,powFac,time,anaId);
}

//------------------------------------------------------------------------------

void OptcritDispLevel::graddu(Structopt *structopt) 
{ 
    structopt->getgraddudisplevel(aFlag,qFlag,size,nodeid,distyp,
                                  refVal,difVal,powFac,time,anaId);
}

//------------------------------------------------------------------------------

void OptcritDispLevel::print(FILE * optunitout) 
{
     fprintf(optunitout,
             "\t%d. Design Criteria  - Type: Displacement Leveling\n",num);

     fprintf(optunitout,"\tAnalysisId: %d   time: %e\n",anaId,time);

     if (aFlag < 2)
       fprintf(optunitout,"\tAverageing      : %d \n",aFlag);
     else
       fprintf(optunitout,"\tNodal difference\n");
          
     fprintf(optunitout,"\tDividing by N   : %d \n",qFlag);
     fprintf(optunitout,"\tExponent        : %e \n",powFac);
     fprintf(optunitout,"\tNumber of nodes : %d \n",size);
     fprintf(optunitout,"\n");

     Disptyp disptyp;
     char *dispname;

     int i;
     for (i=0;i<size;i++) {

       disptyp.typ = distyp[i];
       dispname    = disptyp.gettypname();

       fprintf(optunitout,"\t\tNode                   : %d \n",nodeid[i]+1);      
       fprintf(optunitout,"\t\tDisplacement           : %s \n",dispname);

       fprintf(optunitout,"\t\tReference (Numerator)  : %e \n",refVal[i]);
       
       if (aFlag < 2)
       fprintf(optunitout,"\t\tReference (Denumerator): %e \n",difVal[i]);
       else 
       fprintf(optunitout,"\t\tReference Node         : %d \n",int(difVal[i]));
    }

    fprintf(optunitout,"\n");
}
//------------------------------------------------------------------------------

// stress leveling

OptcritStressLevel::OptcritStressLevel ( int* id, double * xd, dldata * dl, int dlsize )
{ 
   size   = dlsize;

   aFlag  = id[0];
   qFlag  = id[1];

   eFlag  = id[2];

   powFac = xd[0];

   nodeid = new int[dlsize];
   strtyp = new int[dlsize];

   refVal = new double[dlsize];
   difVal = new double[dlsize];

   int i;
   for (i=0;i<dlsize;i++) {
     nodeid[i] = dl[i].node;
     strtyp[i] = dl[i].typ;
     refVal[i] = dl[i].refVal;
     difVal[i] = dl[i].difVal; 
   }
} 

//------------------------------------------------------------------------------

void OptcritStressLevel::evaluate(Structopt *structopt, double tmax, double tmin) 
{
    // Check Time
    if ( ! checkTime(tmax,tmin) ) return; 

    val = structopt->getstresslevel(aFlag,qFlag,eFlag,size,nodeid,strtyp,
                                  refVal,difVal,powFac,time,anaId);
}

//------------------------------------------------------------------------------

void OptcritStressLevel::gradcrit(Structopt *structopt, int actvar) 
{
    grad = structopt->getgradstresslevel(aFlag,qFlag,eFlag,size,nodeid,strtyp,
                                       refVal,difVal,powFac,time,anaId);
}


//------------------------------------------------------------------------------

void OptcritStressLevel::gradpart(Structopt *structopt) 
{
    grad = structopt->getgradpartstresslevel(aFlag,qFlag,eFlag,size,nodeid,strtyp,
                                           refVal,difVal,powFac,time,anaId);
}

//------------------------------------------------------------------------------

void OptcritStressLevel::graddu(Structopt *structopt) 
{ 
    structopt->getgraddustresslevel(aFlag,qFlag,eFlag,size,nodeid,strtyp,
                                  refVal,difVal,powFac,time,anaId);
}

//------------------------------------------------------------------------------

void OptcritStressLevel::print(FILE * optunitout) 
{
    fprintf(optunitout,
   	    "\t%d. Design Criteria  - Type: Stress Leveling\n",num);

    fprintf(optunitout,"\tAnalysisId: %d   time: %e\n",anaId,time);

    fprintf(optunitout,"\tAverageing	     : %d \n",aFlag);
    fprintf(optunitout,"\tDividing by N      : %d \n",qFlag);
    fprintf(optunitout,"\tExponent	     : %e \n",powFac);
    if (!eFlag) {
      fprintf(optunitout,"\tNumber of elements : %d \n",size);
      fprintf(optunitout,"\tElemental stress values are used.\n");
    }
    else {
      fprintf(optunitout,"\tNumber of nodes    : %d \n",size);
      fprintf(optunitout,"\tNodal stress values are used.\n");
    }

    fprintf(optunitout,"\n");

    Stresstyp stresstyp;
    char *stressname;

    int i;
    for (i=0;i<size;i++) {

      stresstyp.typ = strtyp[i];
      stressname    = stresstyp.gettypname();

      if (!eFlag) {
   	fprintf(optunitout,"\t\tElement 	       : %d \n",nodeid[i]+1);	   
   	fprintf(optunitout,"\t\tStress type	       : %s \n",stressname);
   	fprintf(optunitout,"\t\tReference (Numerator)  : %e \n",refVal[i]);
   	fprintf(optunitout,"\t\tReference (Denumerator): %e \n",difVal[i]);
      }
      else {
   	fprintf(optunitout,"\t\tNode		       : %d \n",nodeid[i]+1);	   
   	fprintf(optunitout,"\t\tStress type	       : %s \n",stressname);
   	fprintf(optunitout,"\t\tReference (Numerator)  : %e \n",refVal[i]);
   	fprintf(optunitout,"\t\tReference (Denumerator): %e \n",difVal[i]);
      }
    }

   fprintf(optunitout,"\n");
}

//------------------------------------------------------------------------------

void OptcritExt::evaluate(Structopt *structopt,double tmax,double tmin) 
{
      // Check Time
    if ( ! checkTime(tmax,tmin) ) return; 
    
  val = structopt->getExtVal(num-1,1,anaId,time);  

}

//------------------------------------------------------------------------------

void OptcritExt::gradcrit(Structopt *structopt,int actvar)
{

   grad = structopt->getExtVal(num-1,2,anaId,time);  

} 

//------------------------------------------------------------------------------

void OptcritExt::gradpart(Structopt *structopt)
{
   grad = 0.0;

}

//------------------------------------------------------------------------------

void OptcritExt::graddu(Structopt *structopt)
{
  structopt->getgradduExtVal(num-1,anaId,time); 
 
}

//------------------------------------------------------------------------------

void OptcritExt::print(FILE * optunitout)
{

     char *critname = crittyp.gettypname();

     fprintf(optunitout,"\t%d. Design Criteria  - Type: Local | %s \n ",
             num,critname);

     fprintf(optunitout,"\tAnalysisId: %d   time: %e\n",anaId,time);

}

//------------------------------------------------------------------------------

void OptcritElecMech::evaluate(Structopt *structopt,double tmax,double tmin) 
{

  // Check Time
  if ( ! checkTime(tmax,tmin) ) return; 
  
  int typ = crittyp.typ;

  switch (typ) {
    case 22:
      val = structopt->getPullIn(anaId);
      break;
    default:
      fprintf(stderr,"\n Criteria type %d is not implemented for electro-mechanical criteria. Exiting!!",typ);
      exit(-1);
  }

}

//------------------------------------------------------------------------------

void OptcritElecMech::gradcrit(Structopt *structopt,int actvar)
{

  int typ=crittyp.typ;

  switch (typ) {
    case 22:
      grad = structopt->getgradPullIn(actvar,anaId);
      break;
  }

} 

//------------------------------------------------------------------------------

void OptcritElecMech::gradpart(Structopt *structopt)
{

  int typ=crittyp.typ;

  switch (typ) {
    case 22:
      grad = structopt->getgradpartPullIn(anaId);
      break;
  }

}

//------------------------------------------------------------------------------

void OptcritElecMech::graddu(Structopt *structopt)
{

  grad = 0.0;

}

//------------------------------------------------------------------------------

void OptcritElecMech::print(FILE * optunitout)
{

  char *critname = crittyp.gettypname();

  fprintf(optunitout,"\t%d. Design Criteria  - Type: Local | %s \n ",
          num,critname);

  fprintf(optunitout,"\tAnalysisId: %d   time: %e\n",anaId,time);

}

#endif
