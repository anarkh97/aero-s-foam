#ifndef _MMAGC_H_
#define _MMAGC_H_

#ifdef STRUCTOPT

#include <stdio.h>

#include <Structopt.d/Optsol.h>
//class OptalgNlpmma;

class MMAgc {

   OptalgNlpmma* optalg;

   FILE* ofile;

   double saFac;
   double sbFac;
   double scFac;

   double dstep;

   double acc;

   int    numvar;
   int    numcon;

   int    maxiter;
   int    maxinnerit;

   double epsimin;
   
   double  raa0;
   double  raa0eps;
   double  f0app;
   double  a0;
   double  r0;
   double  f0val;
   double  f0valnew;

   double  z;
   double  zold; 
   double  zmma;
   double  dz;
   double  delz;

   double  zet;
   double  zetold; 
   double  zetmma;
   double  dzet;

   double*  xval;
   double*  xold1;
   double*  xold2;
   double*  xmax;
   double*  xmin;
   
   double*  fval;
   double*  fvalnew;
   double*  df0dx;
   double*  dfdat;
   double** dfdx;
   
   double* raa;
   
   double*  a;
   double*  b;
   double*  c;
   double*  d;
   double*  r;
   
   double*  low;
   double*  upp;

   double*  alfa;
   double*  beta;

   double*  p0;
   double*  q0;
   
   double** P;
   double*  Pd;

   double** Q;
   double*  Qd;

   double** GG;
   double*  GGd;

   double** AA;
   double*  AAd;

   double*  ymma;
   double*  lam;
   double*  lamold;
   double*  dlam;
   double*  lamma;
   double*  xsi;
   double*  xsiold;
   double*  dxsi;
   double*  xsimma;
   double*  eta;
   double*  etaold;
   double*  deta;
   double*  etamma;
   double*  mu;
   double*  muold;
   double*  dmu;
   double*  mumma;
   double*  s;
   double*  sold;
   double*  ds;
   double*  smma;
   double*  fapp;

   double*  x;
   double*  xold;
   double*  dx;
   double*  diagx;	
   double*  diagxinv;
   double*  xmma;
   double*  delx;
   double*  y;
   double*  yold;
   double*  dy;
   double*  dely;
   double*  gvec;
   double*  diagy;
   double*  diagyinv;
   double*  dellam;
   double*  dellamyi;
   double*  diaglam;
   double*  diaglamyi;
   double*  diaglamyiinv;
   double*  bb;
   
   int* active;

  public:
  
    MMAgc(OptalgNlpmma*, 
          double*, double*, double*, int, int,
	  int, int, 
	  double&, double&, double&, double&, double&,
	  FILE*);

    void initialize(double*,double*,double*);
    void asymp(int);
    void gcmmasub(int);
    void raaupdate();
    void subsolve();
    void linsolve(double**,double*,double*,int);
    void cleanup();
    
    int solve();
    int kktcheck();
    
    void func(int,double*,double&,double*);
    void grad(int*,double*,double*,double**);
    
    void printsol(int);

    inline double max(double a, double b) { return a > b ? a : b; }
    inline double min(double a, double b) { return a < b ? a : b; }
    inline double abs(double a) { return a < 0 ? -a : a; }

};

#endif

#endif
