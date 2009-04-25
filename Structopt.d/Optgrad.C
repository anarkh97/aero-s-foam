#ifdef STRUCTOPT

#include <stdio.h>
#include <stdlib.h>

#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optinp.h>
#include <Structopt.d/Optsol.h>
#include <Structopt.d/Optcrit.h>
#include <Structopt.d/Optobj.h>
#include <Structopt.d/Optcon.h>
#include <Structopt.d/Optfilter.h>
#include <Structopt.d/Structopt_sd.h>

#include <Timers.d/GetTime.h>

//------------------------------------------------------------------------------

void Optgrad::buildgrad ( graddata & gd) 
{
    typ       = gd.typ;

    switch (gd.mth) {
    
      case 0:
        anatyp = 0;
	break;
      case 1:
        anatyp = 1;
	break;
      case 2: 
        numtyp = 0;
	break;
      case 3:
        numtyp = 1;
	break;
      default:
        fprintf(stderr,"Error: Gradient Method not correctly specified.\n");          
        exit(-1);
    }	

    epstyp      = gd.epstyp;
    epsval      = gd.epsval;

    icount      = 0;
    actCriteria = 0;
    filterOp    = 0;

    // set up gradient filtering

    if (gd.filter) {
      filterOp   = new Optfilter( gd.filterTyp, gd.filterScale, gd.radius,
                                  gd.maxCount, gd.minExp, gd.maxExp,
				  gd.numGroups, gd.numFilGrps,  gd.filGrpsList);

      numFilCrit = gd.numFilCrit;
      filterCrit = new int[numFilCrit];
    
      int i;
      for(i=0;i<numFilCrit;i++) {
        filterCrit[i] = gd.filcritList[i];
      }
    }
}
     
//------------------------------------------------------------------------------

void Optgrad::setactpro ( Optpro *_optpro, Optsol *_optsol, 
                          Structopt *_structopt) {

    optpro    = _optpro;
    optsol    = _optsol;
    structopt = _structopt;

    if (!actCriteria) actCriteria = new int[optpro->numcrit];

    // build objective/constraint <-> criteria table

    optpro->buildOCC();

    // initialize filter operator

    if (filterOp) {
      filterOp->initialize(structopt,optsol->getVariables(),optsol->numvar);

      //build list for criteria to be filtered

      int* tmp = new int[optpro->numcrit];
      int i;

      for (i=0;i<optpro->numcrit;i++) tmp[i]=0;

      for(i=0;i<numFilCrit;i++) 
        tmp[abs(filterCrit[i])-1] = filterCrit[i]/abs(filterCrit[i]);

      delete [] filterCrit;
      filterCrit = tmp;
    }
}			  

//------------------------------------------------------------------------------

void Optgrad::grad(int* active) {

   double gradTime = getTime();

   switch (typ) {
   
     case 0:
       gradanalytic(active);
       break;
     case 1:
       gradnumeric();
       break;
     default:
       fprintf(stderr,"Type of Sensitivity Analysis not specified");
       exit(-1);
   }

   filePrint(stderr,"\n ... Time spent in gradient evaluation: %e sec\n\n",
	     (getTime()-gradTime)/1000.0);
}  

//------------------------------------------------------------------------------

void Optgrad::gradanalytic(int* active) {

   icount++;

   //Get number of abstract variables and criteria

   int numcrit = optpro->numcrit;
   int numvar  = optsol->numvar;

   //Allocate 2-D gc array for gradients of criteria
   //pgc contains the same information but in a 1-D array

   int i,j;       
     
   double * pgc = new double [numcrit*numvar];
   double ** gc = new double*[numcrit]; // Points to array of size numcrit x numvar 

   for (i=0;i<numcrit;i++) gc[i]=pgc+i*numvar;  

   //Get active criteria to be differentiated 

   optpro->getActiveCriteria(active,actCriteria);  

   //Direct Method

   if (anatyp == 0 ) {

     structopt->zerograd();                //Initialize Arrays for analyt. SA 

     int i;
     for (i=0;i<numvar;i++) {
     
        filePrint(stderr,
		  " ... Analytical SA(D): %d. optimization variable \n",i+1); 

        optpro->optvar->updgrad(i,1.0);    //Set Variation of Abs.&Struct. Var
	 
	structopt->setoptInf(i);           //Set Opt.Var.- Ele. Influence
	
	structopt->graddirect(i,gc);       //Direct Method
     }  
   }
   else {
           
   //Adjoint Method
 
     for (i=0;i<numcrit;i++) {

        if (!actCriteria[i]) continue;

#ifndef AEROELASTIC     
        filePrint(stderr,
		  " ... Analytical SA (A): %d. optimization criteria \r",i+1); 
#endif

        structopt->gradadjoint(i,gc); 
     }	       
   }		                 

   // filter gradients of criteria

   filePrint(stderr,"\n");

   if (filterOp) {

     for (i=0;i<numcrit;i++) { 

       if (actCriteria[i] && filterCrit[i] != 0 ) {
         filePrint(stderr,
		   "\n ... Filtering %d. optimization criteria \n",i+1); 
         filterOp->project(gc[i],filterCrit[i]);
       }
     }
     filterOp->stepCounter();
   }

   // computing derivatives of objectives & constraints taking 
   // into account variable scaling

   for (i=0;i<numvar;i++) {

      for (j=0;j<numcrit;j++)              //Store Deriv. of Criteria
        if (actCriteria[j]) 
          optpro->opc[j]->grad = gc[j][i];

      optpro->optvar->updgrad(i,1.0);      //Set Variation of Abs.&Struct. Var
                                           //for considering explicit functions
	
      optpro->optobj->grad();		   //Evaluation of Deriv. of Objective
  
      optpro->optcon->grad();		   //Evaluation of Deriv. of Constraints
   
      optsol->scalgrad(i);		   //Scaling Deriv. of Object. & Const.
   }  

   //delete local arrays

   delete [] pgc;
   delete [] gc;

#ifndef AEROELASTIC     
   filePrint(stderr,"\n"); 
#endif
}

//------------------------------------------------------------------------------

void Optgrad::gradanalytic(int* actcrit,double** gradmat,int numact,int loctyp) {

   int i,j;       
     
   //Get number of abstract variables and criteria

   int numcrit = optpro->numcrit;
   int numvar  = optsol->numvar;

   //Determine type of approach (direct/adjoint)

   if (loctyp < 0) loctyp=anatyp;

   //Direct Method

   if (loctyp == 0 ) {

     structopt->zerograd();                //Initialize Arrays for analyt. SA 

     for (i=0;i<numvar;i++) {
     
#ifndef AEROELASTIC     
        filePrint(stderr,
		  " ... Analytical SA(D) - multiple criterion: %d. optimization variable \r",i+1); 
#endif

        optpro->optvar->updgrad(i,1.0);     //Set Variation of Abs.&Struct. Var
	
	structopt->setoptInf(i);            //Set Opt.Var.- Ele. Influence
	
	structopt->graddirect(i,gradmat,numact,actcrit);        //Direct Method
     }  
   }
   else {

   //Adjoint Method

     double * pgc = new double [numcrit*numvar];
     double ** gc = new double*[numcrit]; 

     for (i=0;i<numcrit;i++) gc[i]=pgc+i*numvar;  
           
     for (i=0;i<numact;i++) {
 
#ifndef AEROELASTIC     
       filePrint(stderr,
		 " ... Analytical SA (A) - multiple criterion: %d. optimization criteria \r",actcrit[i]+1); 
#endif

        structopt->gradadjoint(actcrit[i],gc); 
     }	
       
     for (i=0;i<numact;i++) 
       for (j=0;j<numvar;j++) 
         gradmat[i][j]=gc[actcrit[i]][j];

     delete [] pgc;
     delete [] gc;
   }		                 

   // filter gradients of criteria

   if (filterOp) 
     filePrint(stderr,"\n ... Filtering for multiple criterion SA ignored \n"); 

}

//------------------------------------------------------------------------------

void Optgrad::gradnumeric() {

  icount++;
  
  int i,j;
  double oldvar;      
  double incvar;
  double oldgrd;

  const int numFluidRes = 7;

  int numvar  = optsol->numvar;   
  int numcon  = optsol->numcon;

  // save current values of optimization problem

  double oldobj;
  double *oldcon = new double[numcon];

  optpro->saveProblem(oldobj,oldcon,numcon);

  for (j=0;j<numcon;j++) { oldcon[j] = optsol->con[j] ; }

#ifdef AEROELASTIC

  int numcrit = optpro->numcrit;

  int    aeroact  = structopt->aeroact; 
  double oldaero [numFluidRes]; 

  if (aeroact) {     
    double aforce;
    for (j=0;j<numFluidRes;j++) oldaero[j]=0.0;
    int icr;

    for (icr=0;icr<numcrit;icr++) {
      int idir=optpro->opc[icr]->getAeroforce(aforce);
      if ( idir >= 0 ) oldaero[idir] = aforce;
    }
  }

#endif  

  for (i=0;i<numvar;i++) {
  
    filePrint(stderr," ... Numerical SA: %d. optimization variable \n",i+1); 

    oldvar  = optsol->var[i];
    
    if (epstyp) {
      incvar  = epsval;
    }
    else {
      incvar  = oldvar * epsval;
    }
          
    optsol->var[i] = oldvar + incvar ; 
    
    optsol->func(0);

    optsol->gradobj[i]= ( optsol->obj - oldobj ) / incvar;

    for (j=0;j<numcon;j++) { 
       optsol->gradcon[i][j]=( optsol->con[j] - oldcon[j] ) / incvar;
    }
    
    if ( numtyp ) {
    
      optsol->var[i] = oldvar - incvar ; 

      optsol->func(0);

      oldgrd = optsol->gradobj[i];
      
      optsol->gradobj[i] = 0.5 * ( oldgrd + ( oldobj - optsol->obj ) / incvar );

      for (j=0;j<numcon;j++) { 

        oldgrd = optsol->gradcon[i][j];
       
        optsol->gradcon[i][j]= 0.5 * ( oldgrd  
	                               + ( oldcon[j] - optsol->con[j] ) / incvar );
      }
    }
    optsol->var[i] = oldvar;
  }
  
  // reset objective and constraints

  optpro->restoreProblem(oldobj,oldcon,numcon);

  delete [] oldcon;
    
  filePrint(stderr,"\n"); 
}
  
//------------------------------------------------------------------------------

void Optgrad::print(FILE * optunitout) {

   char * typname;
   
   switch (typ) {
   
     case 0:
       if (anatyp) { 
         typname="Analytical | Adjoint Method"; }
       else {
         typname="Analytical | Direct Method";  }
       break;
     case 1:
       if (!numtyp)
         typname="Numerical | Forward Difference";
       else
         typname="Numerical | Central Difference";
       break;
     default:
       typname="not specified";
   }    

   filePrint(optunitout,"\n\tGradient Method: %s\n",typname);
   
   if (typ) {
     if (epstyp) {   
       filePrint(optunitout,"\tAbsolute Distrubation: %12.5f\n\n",epsval);   
     }
     else {
       filePrint(optunitout,"\tRelative Distrubation: %12.5f\n\n",epsval);   
     }
   }

   if (filterOp) filterOp->print(filterCrit,optpro->numcrit,optunitout);

}

#endif
