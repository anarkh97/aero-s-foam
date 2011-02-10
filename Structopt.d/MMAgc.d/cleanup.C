#ifdef STRUCTOPT

#include "mma.h"

#include <cstdio>

void MMAgc::cleanup()
{
 delete [] xval;
 delete [] xmax;
 delete [] xmin;
 delete [] xold1;
 delete [] xold2;
 delete [] low;
 delete [] upp;

 delete [] fval;
 delete [] fvalnew;
 delete [] df0dx;
 delete [] dfdat;
 delete [] dfdx;

 delete [] a;
 delete [] b;
 delete [] c;
 delete [] d;
 delete [] r;

 delete [] raa;

 delete [] alfa;
 delete [] beta;

 delete [] p0;
 delete [] q0;

 delete [] Pd;
 delete [] P;

 delete [] Qd;
 delete [] Q;
 
 delete [] GGd;
 delete [] GG;
 
 delete [] AAd;
 delete [] AA;
 delete [] bb;
 
 delete [] x;
 delete [] xold;
 delete [] dx;
 delete [] delx;
 delete [] xmma;
 delete [] diagx;
 delete [] diagxinv;
 delete [] xsi;
 delete [] xsiold;
 delete [] dxsi;
 delete [] xsimma; 
 delete [] eta;
 delete [] etaold;
 delete [] deta;
 delete [] etamma;

 delete [] y;
 delete [] yold;
 delete [] dy;
 delete [] dely;
 delete [] diagy;
 delete [] ymma;
 delete [] diagyinv;
 delete [] lam; 
 delete [] lamold; 
 delete [] dlam; 
 delete [] dellam;
 delete [] dellamyi;
 delete [] lamma;
 delete [] diaglam;
 delete [] diaglamyi;
 delete [] diaglamyiinv;
 delete [] mu;
 delete [] muold;
 delete [] dmu;
 delete [] mumma;
 delete [] s; 
 delete [] sold; 
 delete [] ds; 
 delete [] smma; 
 delete [] gvec;
 
 delete [] fapp; 

 delete [] active;
}

#endif
