#ifdef STRUCTOPT

#include <cstdio>
#include <dlfcn.h>

#include <Structopt.d/Structopt_sd.h>

#include <Structopt.d/Optopr.h>
#include <Structopt.d/Optpro.h>
#include <Structopt.d/Optvar.h>
#include <Structopt.d/Optinp.h>

//------------------------------------------------------------------------------
//                   Optimization Problem 
//                                          created  9/1/98 by Kurt
//------------------------------------------------------------------------------

Optopr * Optopr::create(int _typ) {

   Optopr *opr;
   
   switch (_typ) {
     case 0:
       opr = new OptoprSum;
       opr->typ=0;
       break;
     case 1:
       opr = new OptoprKS;
       opr->typ=1;
       break;
     case 2:
       opr = new OptoprMul;
       opr->typ=2;
       break;
     case 3:
       opr = new OptoprTrig;
       opr->typ=3;
       break;
     case 4:
       opr = new OptoprTrig;
       opr->typ=4;
       break;
     case 5:
       opr = new OptoprExp;
       opr->typ=5;
       break;       
     case 6:
       opr = new OptoprUdef;
       opr->typ=6;
       break;       
     case 7:
       // log
       opr = new OptoprTrig;
       opr->typ=7;
       break;       
    default:
       fprintf(stderr,"Operator Type not implemented: %d\n", _typ );
       exit(-1);
   }
   
   return opr;
}

//------------------------------------------------------------------------------

char * Optopr::gettypname() {

   char *typnam;

   switch (typ) {
     case 0:
       typnam="Sum over i of ( a[i] * x[i]^p[i] + b[i] )";
       break;
     case 1:
       typnam="KS-Function: Log ( Sum over i of e^( a[i] * x[i]^p[i] + b[i] )";
       break;
     case 2:
       typnam="Multiplication over i of ( a[i] * x[i]^p[i] + b[i] )";
       break;
     case 3:
       typnam="Sinus of ( a * x^p + b )";
       break;
     case 4:
       typnam="Cosinus of ( a * x^p + b )";
       break;
     case 5:
       typnam="Exponential of ( a * x^p + b )";
       break;     
     case 6:
       typnam="User-defined Function";
       break;     
     default:
       typnam="not specified";
   }
   
   return typnam;
}

//------------------------------------------------------------------------------

OptoprSum::OptoprSum():opc(0,5),opr(0,5),opv(0,5),
                       indopc(0,5),indopr(0,5),indopv(0,5),
                       a(0,5),b(0,5),p(0,5) {

  numsum=0;

  numopc=0;
  numopr=0;
  numopv=0;
  
}
    
//------------------------------------------------------------------------------

void OptoprSum::build( Optcrit *_opc, double _a, double _p, double _b) {
  
  opc[numopc]=_opc;

  indopc[numopc]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopc++;
  numsum++;

}

//------------------------------------------------------------------------------

void OptoprSum::build( Optopr *_opr, double _a, double _p, double _b) {
  
  opr[numopr]=_opr;

  indopr[numopr]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopr++;
  numsum++;

}

//------------------------------------------------------------------------------

void OptoprSum::build( Absvar *_opv, double _a, double _p, double _b) {
  
  opv[numopv]=_opv;

  indopv[numopv]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopv++;
  numsum++;

}

//------------------------------------------------------------------------------

double OptoprSum::getval() {
  
  double sum=0.0;
  double dummy;
  
  int i,j;

  if (numopr) {
    
    for (i=0;i<numopr;i++) {
       
       dummy = opr[i]->getval();
       j     = indopr[i];
       sum  += a[j]*pow(dummy,p[j]) + b[j];
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
 
       dummy = opc[i]->getval();
       j     = indopc[i];
       sum  += a[j]*pow(dummy,p[j]) + b[j];
    }
  }

  if (numopv) { 
   
    for (i=0;i<numopv;i++) {

       dummy = opv[i]->getval();
       j     = indopv[i];
       sum  += a[j]*pow(dummy,p[j]) + b[j];
    }
  }
  
  return sum;
}

//------------------------------------------------------------------------------

double OptoprSum::getgrad() {
  
  double sum=0.0;
  double dummy,gummy;
  
  int i,j;

  if (numopr) {
    
    for (i=0;i<numopr;i++) {
       
       dummy = opr[i]->getval();
       gummy = opr[i]->getgrad();
       j     = indopr[i];
       sum  += a[j]*p[j]*gummy*pow(dummy,p[j]-1);
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
 
       dummy = opc[i]->getval();
       gummy = opc[i]->getgrad();
       j     = indopc[i];
       sum  += a[j]*p[j]*gummy*pow(dummy,p[j]-1);
    }
  }

  if (numopv) { 
   
    for (i=0;i<numopv;i++) {

       dummy = opv[i]->getval();
       gummy = opv[i]->getgrad();
       j     = indopv[i];
       sum  += a[j]*p[j]*gummy*pow(dummy,p[j]-1);
    }
  }
  
  return sum;
}

//------------------------------------------------------------------------------

void OptoprSum::print(FILE * optunitout) {

  char *typname=gettypname();      
  fprintf(optunitout,"\n\tOperator Type: %s\n",typname);  

  if (numopr) {

    int i;
    int j;
    
    for (i=0;i<numopr;i++) {
   
       j = indopr[i];
       fprintf(optunitout,"\n\t%d. Summand:                 ",j+1);
       fprintf(optunitout,
              "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);

       opr[i]->print(optunitout);
    }
  } 
  
  if (numopc) { 

    int i;
    int j;
    int num;

    for (i=0;i<numopc;i++) {
    
      num=opc[i]->getnumber();
  
      j = indopc[i];
      fprintf(optunitout,"\n\t%d. Criteria Nr %d:          ",j+1,num);
      fprintf(optunitout,
             "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);
    }
  }

  if (numopv) { 
  
    int i;
    int j;
    int num;

    for (i=0;i<numopv;i++) {

      num=opv[i]->getnumber();
  
      j = indopv[i];
      fprintf(optunitout,"\n\t%d. Abstract Variable Nr. %d:",j+1,num);
      fprintf(optunitout,
             "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);
    }
  }
}

//------------------------------------------------------------------------------

void OptoprSum::printres(FILE * optunitout) {

  if (numopr) {

    int i;
    int j;

    for (i=0;i<numopr;i++) {
       j = indopr[i];
       fprintf(optunitout,"\n\t%d. Summand:\n",j+1);
       opr[i]->printres(optunitout);
    }
  } 
  
  if (numopc) { 

    int i;
    int num;
 
    for (i=0;i<numopc;i++) {
      num=opc[i]->getnumber();
      fprintf(optunitout,"\tCriterion Nr.%d: %20.10e\n",num,opc[i]->getval());
    }
  }

  if (numopv) { 
  
    int i;
    int num;

    for (i=0;i<numopv;i++) {
      num=opv[i]->getnumber();
      fprintf(optunitout,"\tAbstract Variable Nr.%d: %12.5e\n",num,opv[i]->getval());
    }
  }
}

//------------------------------------------------------------------------------

funcall & OptoprSum::extractFuncall() {

  funcall * fca = new funcall;

  funcdata * fda = buildFuncdata(); 

  fca->typ   = typ;
  fca->fdata = fda;

  fda->numopr = numopr + numopc + numopv;
  fda->numgen = 0;

  int iopr = 0;
  int i,j;

  if (numopr) {

    for (i=0;i<numopr;i++) {
 
       j = indopr[i];
  
       fda->oprtyp[iopr] = 1;

       fda->oprnum[iopr] = i;
       fda->subfunc[i]   = opr[i]->extractFuncall();

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
    
       j = indopc[i];
  
       fda->oprtyp[iopr] = 0;

       fda->oprnum[iopr] = opc[i]->getnumber();
       fda->oprnum[iopr]--;

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  }

  if (numopv) { 
  
    for (i=0;i<numopv;i++) {

       j = indopv[i];
 
       fda->oprtyp[iopr] = 2;
       fda->oprnum[iopr] = opv[i]->getnumber();
       fda->oprnum[iopr]--;

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  }

  return *fca;
}

//------------------------------------------------------------------------------

void OptoprSum::extractAbsvar(OptActInfo& varInfo) 
{
  int i,j;

  if (numopr)
    for (i=0;i<numopr;i++)
       opr[i]->extractAbsvar(varInfo);
 

  if (numopv) 
    for (i=0;i<numopv;i++) {
       j = opv[i]->getnumber();
       varInfo.add(j-1);   
    }
}

//------------------------------------------------------------------------------

void OptoprSum::extractCriteria(int* critList) 
{
  int i,j;

  if (numopr)
    for (i=0;i<numopr;i++)
       opr[i]->extractCriteria(critList);
 

  if (numopc) 
    for (i=0;i<numopc;i++) {
       j = opc[i]->getnumber();
       critList[j-1]=1;   
    }
}

//------------------------------------------------------------------------------

OptoprKS::OptoprKS():opc(0),opr(0),opv(0),
                     indopc(0),indopr(0),indopv(0),
                     a(0),b(0),p(0) {

  numsum=0;

  numopc=0;
  numopr=0;
  numopv=0;
  
}
    
//------------------------------------------------------------------------------

void OptoprKS::build( Optcrit *_opc, double _a, double _p, double _b) {
  
  opc[numopc]=_opc;

  indopc[numopc]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopc++;
  numsum++;

}

//------------------------------------------------------------------------------

void OptoprKS::build( Optopr *_opr, double _a, double _p, double _b) {
  
  opr[numopr]=_opr;

  indopr[numopr]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopr++;
  numsum++;

}

//------------------------------------------------------------------------------

void OptoprKS::build( Absvar *_opv, double _a, double _p, double _b) {
  
  opv[numopv]=_opv;

  indopv[numopv]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopv++;
  numsum++;

}

//------------------------------------------------------------------------------

double OptoprKS::getval() {
  
  double sum=0.0;
  double dummy;
  
  int i,j;

  if (numopr) {
    
    for (i=0;i<numopr;i++) {
       
       dummy = opr[i]->getval();
       j     = indopr[i];
       sum  += exp( a[j]*pow(dummy,p[j]) + b[j] );
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
 
       dummy = opc[i]->getval();
       j     = indopc[i];
       sum  += exp( a[j]*pow(dummy,p[j]) + b[j] );
    }
  }

  if (numopv) { 
   
    for (i=0;i<numopv;i++) {

       dummy = opv[i]->getval();
       j     = indopv[i];
       sum  += exp( a[j]*pow(dummy,p[j]) + b[j] );
    }
  }
  
  return log(sum);
}

//------------------------------------------------------------------------------

double OptoprKS::getgrad() {
  
  double sum  = 0.0;
  double gsum = 0.0;
  double dummy,gummy,hummy;
  
  int i,j;

  if (numopr) {
    
    for (i=0;i<numopr;i++) {
       
       dummy = opr[i]->getval();
       gummy = opr[i]->getgrad();
       j     = indopr[i];
       hummy = exp( a[j]*pow(dummy,p[j]) + b[j] );
       gsum  += hummy * a[j]*p[j]*gummy*pow(dummy,p[j]-1);
       sum   += hummy;
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
 
       dummy = opc[i]->getval();
       gummy = opc[i]->getgrad();
       j     = indopc[i];
       hummy = exp( a[j]*pow(dummy,p[j]) + b[j] );
       gsum  += hummy * a[j]*p[j]*gummy*pow(dummy,p[j]-1);
       sum   += hummy;
    }
  }

  if (numopv) { 
   
    for (i=0;i<numopv;i++) {

       dummy = opv[i]->getval();
       gummy = opv[i]->getgrad();
       j     = indopv[i];
       hummy = exp( a[j]*pow(dummy,p[j]) + b[j] );
       gsum  += hummy * a[j]*p[j]*gummy*pow(dummy,p[j]-1);
       sum   += hummy;
    }
  }
  
  return (gsum/sum);
}

//------------------------------------------------------------------------------

void OptoprKS::print(FILE * optunitout) {

  char *typname=gettypname();      
  fprintf(optunitout,"\n\tOperator Type: %s\n",typname);  

  if (numopr) {

    int i;
    int j;

    for (i=0;i<numopr;i++) {
   
       j = indopr[i];
       fprintf(optunitout,"\n\t%d. KS-term:                 ",j+1);
       fprintf(optunitout,
              "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);

       opr[i]->print(optunitout);
    }
  } 
  
  if (numopc) { 

    int i;
    int j;
    int num;

    for (i=0;i<numopc;i++) {
    
      num=opc[i]->getnumber();
  
      j = indopc[i];
      fprintf(optunitout,"\n\t%d. Criteria Nr %d:          ",j+1,num);
      fprintf(optunitout,
             "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);
    }
  }

  if (numopv) { 
  
    int i;
    int j;
    int num;

    for (i=0;i<numopv;i++) {

      num=opv[i]->getnumber();
  
      j = indopv[i];
      fprintf(optunitout,"\n\t%d. Abstract Variable Nr. %d:",j+1,num);
      fprintf(optunitout,
             "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);
    }
  }
}

//------------------------------------------------------------------------------

void OptoprKS::printres(FILE * optunitout) {

  if (numopr) {

    int i;
    int j;

    for (i=0;i<numopr;i++) {
       j = indopr[i];
       fprintf(optunitout,"\n\t%d. KS-term:\n",j+1);
       opr[i]->printres(optunitout);
    }
  } 
  
  if (numopc) { 

    int i;
    int num;
 
    for (i=0;i<numopc;i++) {
      num=opc[i]->getnumber();
      fprintf(optunitout,"\tCriterion Nr.%d: %20.10e\n",num,opc[i]->getval());
    }
  }

  if (numopv) { 
  
    int i;
    int num;

    for (i=0;i<numopv;i++) {
      num=opv[i]->getnumber();
      fprintf(optunitout,"\tAbstract Variable Nr.%d: %12.5e\n",num,opv[i]->getval());
    }
  }
}

//------------------------------------------------------------------------------

funcall & OptoprKS::extractFuncall() {

  funcall * fca = new funcall;

  funcdata * fda = buildFuncdata(); 

  fca->typ   = typ;
  fca->fdata = fda;

  fda->numopr = numopr + numopc + numopv;
  fda->numgen = 0;

  int iopr = 0;
  int i,j;

  if (numopr) {

    for (i=0;i<numopr;i++) {
 
       j = indopr[i];
  
       fda->oprtyp[iopr] = 1;

       fda->oprnum[iopr] = i;
       fda->subfunc[i]   = opr[i]->extractFuncall();

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
    
       j = indopc[i];
  
       fda->oprtyp[iopr] = 0;

       fda->oprnum[iopr] = opc[i]->getnumber();
       fda->oprnum[iopr]--;

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  }

  if (numopv) { 
  
    for (i=0;i<numopv;i++) {

       j = indopv[i];
 
       fda->oprtyp[iopr] = 2;
       fda->oprnum[iopr] = opv[i]->getnumber();
       fda->oprnum[iopr]--;

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  }

  return *fca;
}

//------------------------------------------------------------------------------

void OptoprKS::extractAbsvar(OptActInfo& varInfo) 
{
  int i,j;

  if (numopr)
    for (i=0;i<numopr;i++)
       opr[i]->extractAbsvar(varInfo);
 

  if (numopv)
    for (i=0;i<numopv;i++) {
       j = opv[i]->getnumber();
       varInfo.add(j-1);   
    }
}

//------------------------------------------------------------------------------

void OptoprKS::extractCriteria(int* critList) 
{
  int i,j;

  if (numopr)
    for (i=0;i<numopr;i++)
       opr[i]->extractCriteria(critList);
 

  if (numopc) 
    for (i=0;i<numopc;i++) {
       j = opc[i]->getnumber();
       critList[j-1]=1;   
    }
}
//------------------------------------------------------------------------------

OptoprMul::OptoprMul():opc(0),opr(0),opv(0),
                       indopc(0),indopr(0),indopv(0),
                       a(0),b(0),p(0) {

  numsum=0;

  numopc=0;
  numopr=0;
  numopv=0;
  
}
    
//------------------------------------------------------------------------------

void OptoprMul::build( Optcrit *_opc, double _a, double _p, double _b) {
  
  opc[numopc]=_opc;

  indopc[numopc]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopc++;
  numsum++;
  
  if (numopr || numopv) {
    fprintf(stderr,
            "Error: variables in function <mul> have to be of same type\n"); 
    exit(-1);
  }  

}

//------------------------------------------------------------------------------

void OptoprMul::build( Optopr *_opr, double _a, double _p, double _b) {
  
  opr[numopr]=_opr;

  indopr[numopr]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopr++;
  numsum++;

  if (numopc || numopv) {
    fprintf(stderr,
            "Error: variables in function <mul> have to be of same type\n"); 
    exit(-1);
  }  

}

//------------------------------------------------------------------------------

void OptoprMul::build( Absvar *_opv, double _a, double _p, double _b) {
  
  opv[numopv]=_opv;

  indopv[numopv]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopv++;
  numsum++;

  if (numopc || numopr) {
    fprintf(stderr,
            "Error: variables in function <mul> have to be of same type\n"); 
    exit(-1);
  }  
}

//------------------------------------------------------------------------------

double OptoprMul::getval() {
  
  double sum=1.0;
  double dummy;
  
  int i,j;

  if (numopr) {
    
    for (i=0;i<numopr;i++) {
       
       dummy = opr[i]->getval();
       j     = indopr[i];
       sum  *= a[j]*pow(dummy,p[j]) + b[j];
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
 
       dummy = opc[i]->getval();
       j     = indopc[i];
       sum  *= a[j]*pow(dummy,p[j]) + b[j];
    }
  }

  if (numopv) { 
   
    for (i=0;i<numopv;i++) {

       dummy = opv[i]->getval();
       j     = indopv[i];
       sum  *= a[j]*pow(dummy,p[j]) + b[j];
    }
  }
  
  return sum;
}

//------------------------------------------------------------------------------

double OptoprMul::getgrad() {
  
  double sum=0.0;
  double dummy,gummy,summy;
  
  int i,j,k;

  if (numopr) {
    
    for (i=0;i<numopr;i++) {
    
      summy = 1.0;

      for (k=0;k<numopr;k++) {
       
        if ( k==i ) {    
       
          dummy  = opr[k]->getval();
          gummy  = opr[k]->getgrad();
          j      = indopr[k];
          summy *= a[j]*p[j]*gummy*pow(dummy,p[j]-1);
	}
	else {
          dummy  = opr[k]->getval();
          j      = indopr[k];
          summy *= a[j]*pow(dummy,p[j]) + b[j];
	}
      } 
      sum  += summy;
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
 
      summy = 1.0;

      for (k=0;k<numopc;k++) {
       
        if ( k==i ) {    

          dummy  = opc[k]->getval();
          gummy  = opc[k]->getgrad();
          j      = indopc[k];
          summy *= a[j]*p[j]*gummy*pow(dummy,p[j]-1);
	}
	else {
 
          dummy  = opc[k]->getval();
          j      = indopc[k];
          summy *= a[j]*pow(dummy,p[j]) + b[j];
	}
      } 
      sum  += summy;
    }
  } 

  if (numopv) { 
   
    for (i=0;i<numopv;i++) {

      summy = 1.0;

      for (k=0;k<numopv;k++) {
       
        if ( k==i ) {    

          dummy  = opv[k]->getval();
          gummy  = opv[k]->getgrad();
          j      = indopv[k];
          summy *= a[j]*p[j]*gummy*pow(dummy,p[j]-1);
 	}
	else {
 
         dummy  = opv[k]->getval();
         j      = indopv[k];
         summy *= a[j]*pow(dummy,p[j]) + b[j];
	}
      } 
      sum  += summy;
    }
  } 
  
  return sum;
}

//------------------------------------------------------------------------------

void OptoprMul::print(FILE * optunitout) {

  char *typname=gettypname();	   
  fprintf(optunitout,"\n\tOperator Type: %s\n",typname);  
    
  if (numopr) {

    int i;
    int j;

    for (i=0;i<numopr;i++) {
   
       j = indopr[i];
       fprintf(optunitout,"\n\t%d. Term:                    ",j+1);
       fprintf(optunitout,
              "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);

       opr[i]->print(optunitout);
    }
  } 
  
  if (numopc) { 

    int i;
    int j;
    int num;

    for (i=0;i<numopc;i++) {
    
      num=opc[i]->getnumber();
  
      j = indopc[i];
      fprintf(optunitout,"\n\t%d. Criteria Nr %d:          ",j+1,num);
      fprintf(optunitout,
             "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);
    }
  }

  if (numopv) { 
  
    int i;
    int j;
    int num;

    for (i=0;i<numopv;i++) {

      num=opv[i]->getnumber();
  
      j = indopv[i];
      fprintf(optunitout,"\n\t%d. Abstract Variable Nr. %d:",j+1,num);
      fprintf(optunitout,
             "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);
    }
  }
}


//------------------------------------------------------------------------------

void OptoprMul::printres(FILE * optunitout) {

  if (numopr) {

    int i;
    int j;

    for (i=0;i<numopr;i++) {
       j = indopr[i];
       fprintf(optunitout,"\n\t%d. term:\n",j+1);
       opr[i]->printres(optunitout);
    }
  } 
  
  if (numopc) { 

    int i;
    int num;
 
    for (i=0;i<numopc;i++) {
      num=opc[i]->getnumber();
      fprintf(optunitout,"\tCriterion Nr.%d: %20.10e\n",num,opc[i]->getval());
    }
  }

  if (numopv) { 
  
    int i;
    int num;

    for (i=0;i<numopv;i++) {
      num=opv[i]->getnumber();
      fprintf(optunitout,"\tAbstract Variable Nr.%d: %12.5e\n",num,opv[i]->getval());
    }
  }
}

//------------------------------------------------------------------------------

funcall & OptoprMul::extractFuncall() {

  funcall * fca = new funcall;

  funcdata * fda = buildFuncdata(); 

  fca->typ   = typ;
  fca->fdata = fda;

  fda->numopr = numopr + numopc + numopv;
  fda->numgen = 0;

  int iopr = 0;
  int i,j;

  if (numopr) {

    for (i=0;i<numopr;i++) {
 
       j = indopr[i];
  
       fda->oprtyp[iopr] = 1;

       fda->oprnum[iopr] = i;
       fda->subfunc[i]   = opr[i]->extractFuncall();

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
    
       j = indopc[i];
  
       fda->oprtyp[iopr] = 0;

       fda->oprnum[iopr] = opc[i]->getnumber();
       fda->oprnum[iopr]--;

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  }

  if (numopv) { 
  
    for (i=0;i<numopv;i++) {

       j = indopv[i];
 
       fda->oprtyp[iopr] = 2;
       fda->oprnum[iopr] = opv[i]->getnumber();
       fda->oprnum[iopr]--;

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  }

  return *fca;
}

//------------------------------------------------------------------------------

void OptoprMul::extractAbsvar(OptActInfo& varInfo) 
{
  int i,j;

  if (numopr)
    for (i=0;i<numopr;i++)
       opr[i]->extractAbsvar(varInfo);
 

  if (numopv)
    for (i=0;i<numopv;i++) {
       j = opv[i]->getnumber();
       varInfo.add(j-1);   
    }
}

//------------------------------------------------------------------------------

void OptoprMul::extractCriteria(int* critList) 
{
  int i,j;

  if (numopr)
    for (i=0;i<numopr;i++)
       opr[i]->extractCriteria(critList);
 

  if (numopc) 
    for (i=0;i<numopc;i++) {
       j = opc[i]->getnumber();
       critList[j-1]=1;   
    }
}
//------------------------------------------------------------------------------

OptoprTrig::OptoprTrig() {

   opc = 0;
   opr = 0;
   opv = 0;  
}
    
//------------------------------------------------------------------------------

void OptoprTrig::build( Optcrit *_opc, double _a, double _p, double _b) {
  
  opc = _opc;

  a=_a;
  b=_b;
  p=_p;
}

//------------------------------------------------------------------------------

void OptoprTrig::build( Optopr *_opr, double _a, double _p, double _b) {
  
  opr=_opr;

  a=_a;
  b=_b;
  p=_p;
}

//------------------------------------------------------------------------------

void OptoprTrig::build( Absvar *_opv, double _a, double _p, double _b) {
  
  opv=_opv;

  a=_a;
  b=_b;
  p=_p;
}

//------------------------------------------------------------------------------

double OptoprTrig::getval() 
{  
  double sum,dummy;
  if (opr) dummy = opr->getval();  
  if (opc) dummy = opc->getval();
  if (opv) dummy = opv->getval();
  switch (typ) 
    {  
    case 3:
      sum = sin(a*pow(dummy,p) + b); 
      break;    
    case 4:
      sum = cos(a*pow(dummy,p) + b);
      break;
    case 7: 
      dummy = a*pow(dummy,p) + b;
      assert(dummy>0);
      sum = log(dummy);
      break;
    default:
      fprintf(stderr,"OptoprTrig: Operator Type not implemented");
      exit(-1);
    }
  return sum;
}

//------------------------------------------------------------------------------

double OptoprTrig::getgrad() {
  
  double sum=0.0;
  double dummy,gummy,summy;
  
  if (opr) {
    dummy = opr->getval();
    gummy = opr->getgrad();
  }
  
  if (opc) {
    dummy = opc->getval();
    gummy = opc->getgrad();
  }

  if (opv) {
    dummy = opv->getval();
    gummy = opv->getgrad();
  }

  summy = a*p*pow(dummy,p-1)*gummy;

  switch (typ) 
    {  
    case 3:
      sum = summy*cos(a*pow(dummy,p)+ b); 
      break;    
    case 4:
      sum = -summy*sin(a*pow(dummy,p)+ b);
      break;
    case 7:
      dummy = a*pow(dummy,p) + b;
      assert(dummy>0);
      sum = summy/dummy;      
      break;
    default:
      fprintf(stderr,"OptoprTrig: Operator Type not implemented");
      exit(-1);
  }
    
  return sum;
}

//------------------------------------------------------------------------------

void OptoprTrig::print(FILE * optunitout) {

  char *typname=gettypname();      
  fprintf(optunitout,"\n\tOperator Type: %s\n",typname);  
    
  if (opr) {
   
    fprintf(optunitout,"\n\tTerm:                        ");
    fprintf(optunitout,
            "a = %8.3f  p= %8.3f  b= %8.3f\n",a,p,b);
     opr->print(optunitout);
  } 
  
  if (opc) {

    int num=opc->getnumber();
    
    fprintf(optunitout,"\n\tCriteria Nr %d:              ",num);
    fprintf(optunitout,
           "a = %8.3f  p= %8.3f  b= %8.3f\n",a,p,b);
  }

  if (opv) { 
  
    int num=opv->getnumber();
  
    fprintf(optunitout,"\n\tAbstract Variable Nr. %d:    ",num);
    fprintf(optunitout,
           "a = %8.3f  p= %8.3f  b= %8.3f\n",a,p,b);
  }
}


//------------------------------------------------------------------------------

void OptoprTrig::printres(FILE * optunitout) {

  if (opr) {
    fprintf(optunitout,"\n\tterm:\n");
    opr->printres(optunitout);
  } 
  
  if (opc) { 
    int num=opc->getnumber();
    fprintf(optunitout,"\tCriterion Nr.%d: %20.10e\n",num,opc->getval());
  }

  if (opv) {
    int num=opv->getnumber();
    fprintf(optunitout,"\tAbstract Variable Nr.%d: %12.5e\n",num,opv->getval());
  }
}

//------------------------------------------------------------------------------

funcall & OptoprTrig::extractFuncall() {

  funcall * fca = new funcall;

  fprintf(stderr,"ERROR: OptoprTrig::extractFuncall not implemented yet\n");
  exit(-1);

  return *fca;
}

//------------------------------------------------------------------------------

void OptoprTrig::extractAbsvar(OptActInfo& varInfo) 
{
  int j;

  if (opr)
    opr->extractAbsvar(varInfo);
 
  if (opv) {
    j = opv->getnumber();
    varInfo.add(j-1);   
  }
}

//------------------------------------------------------------------------------

void OptoprTrig::extractCriteria(int* critList) 
{
  if (opr)
    opr->extractCriteria(critList);
 
  int j;
  if (opc) { 
    j = opc->getnumber();
    critList[j-1]=1;   
  }
}

//-------------------------------------------------------------------------------

OptoprExp::OptoprExp():opc(0),opr(0),opv(0),
                       indopc(0),indopr(0),indopv(0),
                       a(0),b(0),p(0) {

  numsum=0;

  numopc=0;
  numopr=0;
  numopv=0;
  
}

//-------------------------------------------------------------------------------


void OptoprExp::build( Optcrit *_opc, double _a, double _p, double _b) {
  
  opc[numopc]=_opc;

  indopc[numopc]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopc++;
  numsum++;
  
  if (numopr || numopv) {
    fprintf(stderr,
            "Error: variables in function <exp> have to be of same type\n"); 
    exit(-1);
  }  
}

//------------------------------------------------------------------------------

void OptoprExp::build( Optopr *_opr, double _a, double _p, double _b) {
  
  opr[numopr]=_opr;

  indopr[numopr]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopr++;
  numsum++;

  if (numopc || numopv) {
    fprintf(stderr,
            "Error: variables in function <exp> have to be of same type\n"); 
    exit(-1);
  }  
}

//------------------------------------------------------------------------------

void OptoprExp::build( Absvar *_opv, double _a, double _p, double _b) {
  
  opv[numopv]=_opv;

  indopv[numopv]=numsum;

  a[numsum]=_a;
  b[numsum]=_b;
  p[numsum]=_p;

  numopv++;
  numsum++;

  if (numopc || numopr) {
    fprintf(stderr,
            "Error: variables in function <exp> have to be of same type\n"); 
    exit(-1);
  }  
}

//------------------------------------------------------------------------------

double OptoprExp::getval() {

  double sum=1.0;
  double dummy;
  
  int i,j;

  if (numopr) {
    
    for (i=0;i<numopr;i++) {
       
       dummy = opr[i]->getval();
       j     = indopr[i];
       sum  *= exp( a[j]*pow(dummy,p[j]) + b[j] );
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
 
       dummy = opc[i]->getval();
       j     = indopc[i];
       sum  *= exp( a[j]*pow(dummy,p[j]) + b[j] );
    }
  }

  if (numopv) { 
   
    for (i=0;i<numopv;i++) {

       dummy = opv[i]->getval();
       j     = indopv[i];
       sum  *= exp( a[j]*pow(dummy,p[j]) + b[j] );
    }
  }
  
  return sum;  

}

//------------------------------------------------------------------------------

double OptoprExp::getgrad() {

  double sum=0.0;
  double dummy,gummy,summy;
  
  int i,j,k;

  if (numopr) {
    
    for (i=0;i<numopr;i++) {
    
      summy = 1.0;

      for (k=0;k<numopr;k++) {
       
        if ( k==i ) {    
       
          dummy  = opr[k]->getval();
          gummy  = opr[k]->getgrad();
          j      = indopr[k];
          summy *= (a[j]*p[j]*gummy*pow(dummy,p[j]-1))*exp( a[j]*pow(dummy,p[j]) + b[j] );
	}
	else {
          dummy  = opr[k]->getval();
          j      = indopr[k];
          summy *= exp( a[j]*pow(dummy,p[j]) + b[j] );
	}
      } 
      sum  += summy;
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
 
      summy = 1.0;

      for (k=0;k<numopc;k++) {
       
        if ( k==i ) {    

          dummy  = opc[k]->getval();
          gummy  = opc[k]->getgrad();
          j      = indopc[k];
          summy *= (a[j]*p[j]*gummy*pow(dummy,p[j]-1))*exp( a[j]*pow(dummy,p[j]) + b[j] );
	}
	else {
 
          dummy  = opc[k]->getval();
          j      = indopc[k];
          summy *= exp( a[j]*pow(dummy,p[j]) + b[j] );
	}
      } 
      sum  += summy;
    }
  } 

  if (numopv) { 
   
    for (i=0;i<numopv;i++) {

      summy = 1.0;

      for (k=0;k<numopv;k++) {
       
        if ( k==i ) {    

          dummy  = opv[k]->getval();
          gummy  = opv[k]->getgrad();
          j      = indopv[k];
          summy *= (a[j]*p[j]*gummy*pow(dummy,p[j]-1))*exp( a[j]*pow(dummy,p[j]) + b[j] );
 	}
	else {
 
         dummy  = opv[k]->getval();
         j      = indopv[k];
         summy *= exp( a[j]*pow(dummy,p[j]) + b[j] );
	}
      } 
      sum  += summy;
    }
  } 
  
  return sum;
 
}

//------------------------------------------------------------------------------

void OptoprExp::print(FILE * optunitout) {

  char *typname=gettypname();	   
  fprintf(optunitout,"\n\tOperator Type: %s\n",typname);  
    
  if (numopr) {

    int i;
    int j;

    for (i=0;i<numopr;i++) {
   
       j = indopr[i];
       fprintf(optunitout,"\n\t%d. Term:                    ",j+1);
       fprintf(optunitout,
              "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);

       opr[i]->print(optunitout);
    }
  } 
  
  if (numopc) { 

    int i;
    int j;
    int num;

    for (i=0;i<numopc;i++) {
    
      num=opc[i]->getnumber();
  
      j = indopc[i];
      fprintf(optunitout,"\n\t%d. Criteria Nr %d:          ",j+1,num);
      fprintf(optunitout,
             "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);
    }
  }

  if (numopv) { 
  
    int i;
    int j;
    int num;

    for (i=0;i<numopv;i++) {

      num=opv[i]->getnumber();
  
      j = indopv[i];
      fprintf(optunitout,"\n\t%d. Abstract Variable Nr. %d:",j+1,num);
      fprintf(optunitout,
             "a = %8.3f  p= %8.3f  b= %8.3f\n",a[j],p[j],b[j]);
    }
  }
}

//---------------------------------------------------------------------------

void OptoprExp::printres(FILE * optunitout) {

  if (numopr) {

    int i;
    int j;

    for (i=0;i<numopr;i++) {
       j = indopr[i];
       fprintf(optunitout,"\n\t%d. term:\n",j+1);
       opr[i]->printres(optunitout);
    }
  } 
  
  if (numopc) { 

    int i;
    int num;
 
    for (i=0;i<numopc;i++) {
      num=opc[i]->getnumber();
      fprintf(optunitout,"\tCriterion Nr.%d: %20.10e\n",num,opc[i]->getval());
    }
  }

  if (numopv) { 
  
    int i;
    int num;

    for (i=0;i<numopv;i++) {
      num=opv[i]->getnumber();
      fprintf(optunitout,"\tAbstract Variable Nr.%d: %12.5e\n",num,opv[i]->getval());
    }
  }
}

//------------------------------------------------------------------------------

funcall & OptoprExp::extractFuncall() {

  funcall * fca = new funcall;

  funcdata * fda = buildFuncdata(); 

  fca->typ   = typ;
  fca->fdata = fda;

  fda->numopr = numopr + numopc + numopv;
  fda->numgen = 0;

  int iopr = 0;
  int i,j;

  if (numopr) {

    for (i=0;i<numopr;i++) {
 
       j = indopr[i];
  
       fda->oprtyp[iopr] = 1;

       fda->oprnum[iopr] = i;
       fda->subfunc[i]   = opr[i]->extractFuncall();

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
    
       j = indopc[i];
  
       fda->oprtyp[iopr] = 0;

       fda->oprnum[iopr] = opc[i]->getnumber();
       fda->oprnum[iopr]--;

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  }

  if (numopv) { 
  
    for (i=0;i<numopv;i++) {

       j = indopv[i];
 
       fda->oprtyp[iopr] = 2;
       fda->oprnum[iopr] = opv[i]->getnumber();
       fda->oprnum[iopr]--;

       fda->a[iopr]      = a[j];
       fda->p[iopr]      = p[j];
       fda->b[iopr]      = b[j];
       iopr++;
    }
  }

  return *fca;
}

//-----------------------------------------------------------------------------

void OptoprExp::extractAbsvar(OptActInfo& varInfo) 
{
  int i,j;

  if (numopr)
    for (i=0;i<numopr;i++)
       opr[i]->extractAbsvar(varInfo);
 

  if (numopv)
    for (i=0;i<numopv;i++) {
       j = opv[i]->getnumber();
       varInfo.add(j-1);   
    }
}

//------------------------------------------------------------------------------

void OptoprExp::extractCriteria(int* critList) 
{
  int i,j;

  if (numopr)
    for (i=0;i<numopr;i++)
       opr[i]->extractCriteria(critList);
 

  if (numopc) 
    for (i=0;i<numopc;i++) {
       j = opc[i]->getnumber();
       critList[j-1]=1;   
    }
}

//-------------------------------------------------------------------------------

OptoprUdef::OptoprUdef():opc(0),opr(0),opv(0),
                       indopc(0),indopr(0),indopv(0),
                       indpar(0) 
{
  // initialize counters

  idUdef=-1;

  numsum=0;

  numopc=0;
  numopr=0;
  numopv=0;

  // open shared object
  
  char* udefname = "./userfunctions.so";

  handle = dlopen(udefname, RTLD_NOW);

  char *msg = dlerror();

  if (msg) {
    fprintf(stderr,"*** Error: dynamic loading of \'%s\': %s\n", udefname, msg);
    exit(-1);
  }

}

//-------------------------------------------------------------------------------


void OptoprUdef::build( Optcrit *_opc, double _a, double _p, double _b) {
  
  opc[numopc]=_opc;

  indopc[numopc]=numsum;

  indpar[numsum]=static_cast<int> (_a);

  if (!numsum) {
    idUdef=static_cast<int> (_b);
    getUdef(); 
  }
  else {
    if (idUdef != static_cast<int> (_b)) {
      fprintf(stderr,
            "Error: all operands need to belong to same user function\n"); 
      exit(-1);
    }
  }  
    
  numopc++;
  numsum++;
  
}

//------------------------------------------------------------------------------

void OptoprUdef::build( Optopr *_opr, double _a, double _p, double _b) {
  
  opr[numopr]=_opr;

  indopr[numopr]=numsum;

  indpar[numsum]=static_cast<int> (_a);

  if (!numsum) {
    idUdef=static_cast<int> (_b);
    getUdef(); 
  }
  else {
    if (idUdef != static_cast<int> (_b)) {
      fprintf(stderr,
            "Error: all operands need to belong to same user function\n"); 
      exit(-1);
    }
  }  

  numopr++;
  numsum++;

}

//------------------------------------------------------------------------------

void OptoprUdef::build( Absvar *_opv, double _a, double _p, double _b) {
  
  opv[numopv]=_opv;

  indopv[numopv]=numsum;

  indpar[numsum]=static_cast<int> (_a);

  if (!numsum) {
    idUdef=static_cast<int> (_b);
    getUdef(); 
  }
  else {
    if (idUdef != static_cast<int> (_b)) {
      fprintf(stderr,
            "Error: all operands need to belong to same user function\n"); 
      exit(-1);
    }
  }  

  numopv++;
  numsum++;

}

//------------------------------------------------------------------------------

double OptoprUdef::getval() {


  //build parameter list
  
  double* vals = static_cast<double*>(dbg_alloca(sizeof(double)*numsum));
  double* grds = static_cast<double*>(dbg_alloca(sizeof(double)*numsum));
  
  int i,j;

  if (numopr) {

    for (i=0;i<numopr;i++) {
       j       = indpar[indopr[i]]-1;
       vals[j] = opr[i]->getval();
       grds[j] = 0;
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
 
       j       = indpar[indopc[i]]-1;
       vals[j] = opc[i]->getval();
       grds[j] = 0;
    }
  }

  if (numopv) { 
   
    for (i=0;i<numopv;i++) {

       j       = indpar[indopv[i]]-1;
       vals[j] = opv[i]->getval();
       grds[j] = 0;
    }
  }

  // call user-defined function

  int needGrad = 0;
  
  double sum = fct(vals,grds,numsum,needGrad);
  
  return sum;  

}

//------------------------------------------------------------------------------

double OptoprUdef::getgrad() {

  //build parameter list
  
  double* vals = static_cast<double*>(dbg_alloca(sizeof(double)*numsum));
  double* grds = static_cast<double*>(dbg_alloca(sizeof(double)*numsum));
  
  int i,j;

  if (numopr) {

    for (i=0;i<numopr;i++) {
       j       = indpar[indopr[i]]-1;
       vals[j] = opr[i]->getval();
       grds[j] = opr[i]->getgrad();
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
 
       j       = indpar[indopc[i]]-1;
       vals[j] = opc[i]->getval();
       grds[j] = opc[i]->getgrad();
    }
  }

  if (numopv) { 
   
    for (i=0;i<numopv;i++) {

       j       = indpar[indopv[i]]-1;
       vals[j] = opv[i]->getval();
       grds[j] = opv[i]->getgrad();
    }
  }

  // call user-defined function

  int needGrad = 1;

  double sum = fct(vals,grds,numsum,needGrad);
  
  return sum;  

}

//------------------------------------------------------------------------------

void OptoprUdef::print(FILE * optunitout) {

  char *typname=gettypname();	   
  fprintf(optunitout,"\n\tOperator Type: %s\n",typname);  
    
  if (numopr) {

    int i;
    int j;

    for (i=0;i<numopr;i++) {
   
       j = indopr[i];
       fprintf(optunitout,"\n\t%d. Term:                    ",j+1);
       fprintf(optunitout,
              "User-function Nr. %d  param-id = %d\n",idUdef,indpar[j]);

       opr[i]->print(optunitout);
    }
  } 
  
  if (numopc) { 

    int i;
    int j;
    int num;

    for (i=0;i<numopc;i++) {
    
      num=opc[i]->getnumber();
  
      j = indopc[i];
      fprintf(optunitout,"\n\t%d. Criteria Nr %d:          ",j+1,num);
      fprintf(optunitout,
              "User-function Nr. %d  param-id = %d\n",idUdef,indpar[j]);
    }
  }

  if (numopv) { 
  
    int i;
    int j;
    int num;

    for (i=0;i<numopv;i++) {

      num=opv[i]->getnumber();
  
      j = indopv[i];
      fprintf(optunitout,"\n\t%d. Abstract Variable Nr. %d:",j+1,num);
      fprintf(optunitout,
              "User-function Nr. %d  param-id = %d\n",idUdef,indpar[j]);
    }
  }
}

//---------------------------------------------------------------------------

void OptoprUdef::printres(FILE * optunitout) {

  if (numopr) {

    int i;
    int j;

    for (i=0;i<numopr;i++) {
       j = indopr[i];
       fprintf(optunitout,"\n\t%d. term:\n",j+1);
       opr[i]->printres(optunitout);
    }
  } 
  
  if (numopc) { 

    int i;
    int num;
 
    for (i=0;i<numopc;i++) {
      num=opc[i]->getnumber();
      fprintf(optunitout,"\tCriterion Nr.%d: %20.10e\n",num,opc[i]->getval());
    }
  }

  if (numopv) { 
  
    int i;
    int num;

    for (i=0;i<numopv;i++) {
      num=opv[i]->getnumber();
      fprintf(optunitout,"\tAbstract Variable Nr.%d: %12.5e\n",num,opv[i]->getval());
    }
  }
}

//------------------------------------------------------------------------------

funcall & OptoprUdef::extractFuncall() {

  funcall * fca = new funcall;

  funcdata * fda = buildFuncdata(); 

  fca->typ   = typ;
  fca->fdata = fda;

  fda->numopr = numopr + numopc + numopv;
  fda->numgen = 0;

  int iopr = 0;
  int i,j;

  if (numopr) {

    for (i=0;i<numopr;i++) {
 
       j = indopr[i];
  
       fda->oprtyp[iopr] = 1;

       fda->oprnum[iopr] = i;
       fda->subfunc[i]   = opr[i]->extractFuncall();

       fda->a[iopr]      = static_cast<double> (indpar[j]);
       fda->p[iopr]      = 0.0;
       fda->b[iopr]      = static_cast<double> (idUdef);
       iopr++;
    }
  } 
  
  if (numopc) { 
 
    for (i=0;i<numopc;i++) {
    
       j = indopc[i];
  
       fda->oprtyp[iopr] = 0;

       fda->oprnum[iopr] = opc[i]->getnumber();
       fda->oprnum[iopr]--;

       fda->a[iopr]      = static_cast<double> (indpar[j]);
       fda->p[iopr]      = 0.0;
       fda->b[iopr]      = static_cast<double> (idUdef);
       iopr++;
    }
  }

  if (numopv) { 
  
    for (i=0;i<numopv;i++) {

       j = indopv[i];
 
       fda->oprtyp[iopr] = 2;
       fda->oprnum[iopr] = opv[i]->getnumber();
       fda->oprnum[iopr]--;

       fda->a[iopr]      = static_cast<double> (indpar[j]);
       fda->p[iopr]      = 0.0;
       fda->b[iopr]      = static_cast<double> (idUdef);
       iopr++;
    }
  }

  return *fca;
}

//-----------------------------------------------------------------------------

void OptoprUdef::extractAbsvar(OptActInfo& varInfo) 
{
  int i,j;

  if (numopr)
    for (i=0;i<numopr;i++)
       opr[i]->extractAbsvar(varInfo);
 

  if (numopv)
    for (i=0;i<numopv;i++) {
       j = opv[i]->getnumber();
       varInfo.add(j-1);   
    }
}

//------------------------------------------------------------------------------

void OptoprUdef::extractCriteria(int* critList) 
{
  int i,j;

  if (numopr)
    for (i=0;i<numopr;i++)
       opr[i]->extractCriteria(critList);
 

  if (numopc) 
    for (i=0;i<numopc;i++) {
       j = opc[i]->getnumber();
       critList[j-1]=1;   
    }
}
//------------------------------------------------------------------------------

void OptoprUdef::getUdef() 
{
  char* udefname = "./userfunctions.so";

  char * routine  = static_cast<char*>(malloc(sizeof(char*)*(10)));

  sprintf(routine,"userfunc_%1d",idUdef); 

  //fct = reinterpret_cast<double (*)(double*,double*,int,int)>(dlsym(handle,routine));
  fct = (double (*)(double*,double*,int,int)) (dlsym(handle,routine));
  if (!fct) {
    fprintf(stderr,"*** Error: could not find \'%s\' in \'%s\'\n", routine, udefname);
    exit(-1);
  }
}

#endif
