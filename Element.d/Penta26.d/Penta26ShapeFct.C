// ---------------------------------------------------------------------
// HB - 05-22-05
// ---------------------------------------------------------------------
// Shape fcts & derivatives for 26 nodes wedge element 
// Serendipity finite element basis
// Iso-parametric formulation
// ---------------------------------------------------------------------
// EXPERIMENTAL ...
// ---------------------------------------------------------------------
// Std C/C++ lib
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <Utils.d/linkfc.h>

#define USE_NAVY_PENTA26  //HB: force using the shape function routines provided by the Navy ...

// HB (05-21-05): for given reference shape fct derivatives dShape & nodes coordinates (X.Y,Z), 
//                compute the jacobian & shape fcts derivative w.r.t to the global coordinate system.
//                assuming an iso-parametric formulation.
double
computePenta26DShapeFct(double dShape[26][3], double X[26], double Y[26], double Z[26], double (*DShape)[3]=0)
{
  double xd1,xd2,xd3;
  double yd1,yd2,yd3;
  double zd1,zd2,zd3;

  xd1 = 0.0 ; xd2 = 0.0 ; xd3 = 0.0 ;
  yd1 = 0.0 ; yd2 = 0.0 ; yd3 = 0.0 ;
  zd1 = 0.0 ; zd2 = 0.0 ; zd3 = 0.0 ;

  for(int i=0; i<26; i+=2){
    xd1 += dShape[i][0]*X[i] + dShape[i+1][0]*X[i+1];
    xd2 += dShape[i][1]*X[i] + dShape[i+1][1]*X[i+1];
    xd3 += dShape[i][2]*X[i] + dShape[i+1][2]*X[i+1];
    
    yd1 += dShape[i][0]*Y[i] + dShape[i+1][0]*Y[i+1];
    yd2 += dShape[i][1]*Y[i] + dShape[i+1][1]*Y[i+1];
    yd3 += dShape[i][2]*Y[i] + dShape[i+1][2]*Y[i+1];
    
    zd1 += dShape[i][0]*Z[i] + dShape[i+1][0]*Z[i+1];
    zd2 += dShape[i][1]*Z[i] + dShape[i+1][1]*Z[i+1];
    zd3 += dShape[i][2]*Z[i] + dShape[i+1][2]*Z[i+1];
  }
  //fprintf(stderr," %e  %e  %e \n %e  %e  %e \n %e  %e  %e \n",xd1,xd2,xd3,yd1,yd2,yd3,zd1,zd2,zd3);
                                                                                                                       
  double a11 = yd2*zd3 - yd3*zd2 ; double a12 = yd3*zd1 - yd1*zd3 ; double a13 = yd1*zd2 - yd2*zd1 ;
  double a21 = xd3*zd2 - xd2*zd3 ; double a22 = xd1*zd3 - zd1*xd3 ; double a23 = xd2*zd1 - xd1*zd2 ;
  double a31 = xd2*yd3 - xd3*yd2 ; double a32 = yd1*xd3 - xd1*yd3 ; double a33 = xd1*yd2 - yd1*xd2 ;
 
  /* -------> DETERMINANT OF THE JACOBIAN MATRIX <--- */

  double J = xd1*a11 + yd1*a21 + zd1*a31 ;
  if ( J == 0.0 ) {
    fprintf(stderr," *** WARNING: NULL JACOBIAN IN computePenta26DShapeFct.C ROUTINE.\n");
  }
  if ( J < 0.0 ) {
    fprintf(stderr," *** WARNING: NEGATIVE JACOBIAN IN computePenta26DShapeFct.C ROUTINE.\n");
    //J *= -1.; // BAD !!! BUT WORKS !!!
  }
  double cdet = 1.0 / J ;

  /* -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REAL" ELEMENT <--- */
  if(DShape){
    for(int i=0; i<26; i+=2){
      DShape[i  ][0] = cdet * ( a11*dShape[i  ][0] + a12*dShape[i  ][1] + a13*dShape[i  ][2] );
      DShape[i  ][1] = cdet * ( a21*dShape[i  ][0] + a22*dShape[i  ][1] + a23*dShape[i  ][2] );
      DShape[i  ][2] = cdet * ( a31*dShape[i  ][0] + a32*dShape[i  ][1] + a33*dShape[i  ][2] );
                                                                                                            
      DShape[i+1][0] = cdet * ( a11*dShape[i+1][0] + a12*dShape[i+1][1] + a13*dShape[i+1][2] );
      DShape[i+1][1] = cdet * ( a21*dShape[i+1][0] + a22*dShape[i+1][1] + a23*dShape[i+1][2] );
      DShape[i+1][2] = cdet * ( a31*dShape[i+1][0] + a32*dShape[i+1][1] + a33*dShape[i+1][2] );
    } 
  }

  return(J);
}

#ifdef USE_NAVY_PENTA26
static int GMap_cubicLagrangeRegion_wedge(const double r,const double s,
					  const double zeta, double N[26])
{
  N[ 0] = (3.0*r-1.0)*(3.0*r-2.0)*r*(1.0-zeta)/4.0-9.0/8.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)*r+9.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0)*r;
  N[ 1] = (3.0*s-1.0)*(3.0*s-2.0)*s*(1.0-zeta)/4.0-9.0/8.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)*s+9.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0)*s;
  N[ 2] = (2.0-3.0*r-3.0*s)*(1.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0-zeta)/4.0-9.0/8.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)*(1.0-r-s)+9.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0)*(1.0-r-s);
  N[ 3] = (3.0*r-1.0)*(3.0*r-2.0)*r*(1.0+zeta)/4.0-9.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)*r+9.0/8.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0)*r;
  N[ 4] = (3.0*s-1.0)*(3.0*s-2.0)*s*(1.0+zeta)/4.0-9.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)*s+9.0/8.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0)*s;
  N[ 5] = (2.0-3.0*r-3.0*s)*(1.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0+zeta)/4.0-9.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)*(1.0-r-s)+9.0/8.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0)*(1.0-r-s);
  N[ 6] = 9.0/4.0*r*s*(3.0*r-1.0)*(1.0-zeta);
  N[ 7] = 9.0/4.0*r*s*(3.0*s-1.0)*(1.0-zeta);
  N[ 8] = 9.0/4.0*(1.0-r-s)*s*(3.0*s-1.0)*(1.0-zeta);
  N[ 9] = 9.0/4.0*(1.0-r-s)*s*(2.0-3.0*r-3.0*s)*(1.0-zeta);
  N[10] = 9.0/4.0*r*(1.0-r-s)*(2.0-3.0*r-3.0*s)*(1.0-zeta);
  N[11] = 9.0/4.0*r*(1.0-r-s)*(3.0*r-1.0)*(1.0-zeta);
  N[12] = 9.0/4.0*r*s*(3.0*r-1.0)*(1.0+zeta);
  N[13] = 9.0/4.0*r*s*(3.0*s-1.0)*(1.0+zeta);
  N[14] = 9.0/4.0*(1.0-r-s)*s*(3.0*s-1.0)*(1.0+zeta);
  N[15] = 9.0/4.0*(1.0-r-s)*s*(2.0-3.0*r-3.0*s)*(1.0+zeta);
  N[16] = 9.0/4.0*r*(1.0-r-s)*(2.0-3.0*r-3.0*s)*(1.0+zeta);
  N[17] = 9.0/4.0*r*(1.0-r-s)*(3.0*r-1.0)*(1.0+zeta);
  N[18] = 27.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)*r;
  N[19] = -27.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0)*r;
  N[20] = 27.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)*s;
  N[21] = -27.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0)*s;
  N[22] = 27.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)*(1.0-r-s);
  N[23] = -27.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0)*(1.0-r-s);
  N[24] = 27.0/2.0*r*s*(1.0-r-s)*(1.0-zeta);
  N[25] = 27.0/2.0*r*s*(1.0-r-s)*(1.0+zeta);
  return(26);
}

static int GMap_cubicLagrangeDrvRegion_wedge(const double r,const double s,
					     const double zeta,
					     double dN[26][3]){
  dN[ 0][0] = 3.0/4.0*(3.0*r-2.0)*r*(1.0-zeta)+3.0/4.0*(3.0*r-1.0)*r*(1.0-zeta)+(3.0*r-1.0)*(3.0*r-2.0)*(1.0-zeta)/4.0-9.0/8.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)+9.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[ 0][1] = 0.0;
  dN[ 0][2] = -(3.0*r-1.0)*(3.0*r-2.0)*r/4.0-9.0/8.0*(zeta-1.0/3.0)*(zeta-1.0)*r-9.0/16.0*(1.0+zeta)*(zeta-1.0)*r-9.0/8.0*(1.0+zeta)*(zeta-1.0/3.0)*r+9.0/16.0*(zeta+1.0/3.0)*(zeta-1.0)*r+9.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*r;
  dN[ 1][0] = 0.0;
  dN[ 1][1] = 3.0/4.0*(3.0*s-2.0)*s*(1.0-zeta)+3.0/4.0*(3.0*s-1.0)*s*(1.0-zeta)+(3.0*s-1.0)*(3.0*s-2.0)*(1.0-zeta)/4.0-9.0/8.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)+9.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[ 1][2] = -(3.0*s-1.0)*(3.0*s-2.0)*s/4.0-9.0/8.0*(zeta-1.0/3.0)*(zeta-1.0)*s-9.0/16.0*(1.0+zeta)*(zeta-1.0)*s-9.0/8.0*(1.0+zeta)*(zeta-1.0/3.0)*s+9.0/16.0*(zeta+1.0/3.0)*(zeta-1.0)*s+9.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*s;
  dN[ 2][0] = -3.0/4.0*(1.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0-zeta)-3.0/4.0*(2.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0-zeta)-(2.0-3.0*r-3.0*s)*(1.0-3.0*r-3.0*s)*(1.0-zeta)/4.0+9.0/8.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)-9.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[ 2][1] = -3.0/4.0*(1.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0-zeta)-3.0/4.0*(2.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0-zeta)-(2.0-3.0*r-3.0*s)*(1.0-3.0*r-3.0*s)*(1.0-zeta)/4.0+9.0/8.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)-9.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[ 2][2] = -(2.0-3.0*r-3.0*s)*(1.0-3.0*r-3.0*s)*(1.0-r-s)/4.0-9.0/8.0*(zeta-1.0/3.0)*(zeta-1.0)*(1.0-r-s)-9.0/16.0*(1.0+zeta)*(zeta-1.0)*(1.0-r-s)-9.0/8.0*(1.0+zeta)*(zeta-1.0/3.0)*(1.0-r-s)+9.0/16.0*(zeta+1.0/3.0)*(zeta-1.0)*(1.0-r-s)+9.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(1.0-r-s);
  dN[ 3][0] = 3.0/4.0*(3.0*r-2.0)*r*(1.0+zeta)+3.0/4.0*(3.0*r-1.0)*r*(1.0+zeta)+(3.0*r-1.0)*(3.0*r-2.0)*(1.0+zeta)/4.0-9.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)+9.0/8.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[ 3][1] = 0.0;
  dN[ 3][2] = (3.0*r-1.0)*(3.0*r-2.0)*r/4.0-9.0/16.0*(zeta-1.0/3.0)*(zeta-1.0)*r+9.0/16.0*(1.0+zeta)*(zeta-1.0)*r-9.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*r+9.0/8.0*(zeta+1.0/3.0)*(zeta-1.0)*r+9.0/8.0*(1.0+zeta)*(zeta+1.0/3.0)*r;
  dN[ 4][0] = 0.0;
  dN[ 4][1] = 3.0/4.0*(3.0*s-2.0)*s*(1.0+zeta)+3.0/4.0*(3.0*s-1.0)*s*(1.0+zeta)+(3.0*s-1.0)*(3.0*s-2.0)*(1.0+zeta)/4.0-9.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)+9.0/8.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[ 4][2] = (3.0*s-1.0)*(3.0*s-2.0)*s/4.0-9.0/16.0*(zeta-1.0/3.0)*(zeta-1.0)*s+9.0/16.0*(1.0+zeta)*(zeta-1.0)*s-9.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*s+9.0/8.0*(zeta+1.0/3.0)*(zeta-1.0)*s+9.0/8.0*(1.0+zeta)*(zeta+1.0/3.0)*s;
  dN[ 5][0] = -3.0/4.0*(1.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0+zeta)-3.0/4.0*(2.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0+zeta)-(2.0-3.0*r-3.0*s)*(1.0-3.0*r-3.0*s)*(1.0+zeta)/4.0+9.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)-9.0/8.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[ 5][1] = -3.0/4.0*(1.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0+zeta)-3.0/4.0*(2.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0+zeta)-(2.0-3.0*r-3.0*s)*(1.0-3.0*r-3.0*s)*(1.0+zeta)/4.0+9.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0)-9.0/8.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[ 5][2] = (2.0-3.0*r-3.0*s)*(1.0-3.0*r-3.0*s)*(1.0-r-s)/4.0-9.0/16.0*(zeta-1.0/3.0)*(zeta-1.0)*(1.0-r-s)+9.0/16.0*(1.0+zeta)*(zeta-1.0)*(1.0-r-s)-9.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(1.0-r-s)+9.0/8.0*(zeta+1.0/3.0)*(zeta-1.0)*(1.0-r-s)+9.0/8.0*(1.0+zeta)*(zeta+1.0/3.0)*(1.0-r-s);
  dN[ 6][0] = 9.0/4.0*s*(3.0*r-1.0)*(1.0-zeta)+27.0/4.0*r*s*(1.0-zeta);
  dN[ 6][1] = 9.0/4.0*(3.0*r-1.0)*r*(1.0-zeta);
  dN[ 6][2] = -9.0/4.0*r*s*(3.0*r-1.0);
  dN[ 7][0] = 9.0/4.0*(3.0*s-1.0)*s*(1.0-zeta);
  dN[ 7][1] = 9.0/4.0*r*(3.0*s-1.0)*(1.0-zeta)+27.0/4.0*r*s*(1.0-zeta);
  dN[ 7][2] = -9.0/4.0*r*s*(3.0*s-1.0);
  dN[ 8][0] = -9.0/4.0*(3.0*s-1.0)*s*(1.0-zeta);
  dN[ 8][1] = -9.0/4.0*(3.0*s-1.0)*s*(1.0-zeta)+9.0/4.0*(1.0-r-s)*(3.0*s-1.0)*(1.0-zeta)+27.0/4.0*(1.0-r-s)*s*(1.0-zeta);
  dN[ 8][2] = -9.0/4.0*(1.0-r-s)*s*(3.0*s-1.0);
  dN[ 9][0] = -9.0/4.0*s*(2.0-3.0*r-3.0*s)*(1.0-zeta)-27.0/4.0*(1.0-r-s)*s*(1.0-zeta);
  dN[ 9][1] = -9.0/4.0*s*(2.0-3.0*r-3.0*s)*(1.0-zeta)+9.0/4.0*(2.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0-zeta)-27.0/4.0*(1.0-r-s)*s*(1.0-zeta);
  dN[ 9][2] = -9.0/4.0*(1.0-r-s)*s*(2.0-3.0*r-3.0*s);
  dN[10][0] = 9.0/4.0*(2.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0-zeta)-9.0/4.0*r*(2.0-3.0*r-3.0*s)*(1.0-zeta)-27.0/4.0*r*(1.0-r-s)*(1.0-zeta);
  dN[10][1] = -9.0/4.0*r*(2.0-3.0*r-3.0*s)*(1.0-zeta)-27.0/4.0*r*(1.0-r-s)*(1.0-zeta);
  dN[10][2] = -9.0/4.0*r*(1.0-r-s)*(2.0-3.0*r-3.0*s);
  dN[11][0] = 9.0/4.0*(1.0-r-s)*(3.0*r-1.0)*(1.0-zeta)-9.0/4.0*(3.0*r-1.0)*r*(1.0-zeta)+27.0/4.0*r*(1.0-r-s)*(1.0-zeta);
  dN[11][1] = -9.0/4.0*(3.0*r-1.0)*r*(1.0-zeta);
  dN[11][2] = -9.0/4.0*r*(1.0-r-s)*(3.0*r-1.0);
  dN[12][0] = 9.0/4.0*s*(3.0*r-1.0)*(1.0+zeta)+27.0/4.0*r*s*(1.0+zeta);
  dN[12][1] = 9.0/4.0*(3.0*r-1.0)*r*(1.0+zeta);
  dN[12][2] = 9.0/4.0*r*s*(3.0*r-1.0);
  dN[13][0] = 9.0/4.0*(3.0*s-1.0)*s*(1.0+zeta);
  dN[13][1] = 9.0/4.0*r*(3.0*s-1.0)*(1.0+zeta)+27.0/4.0*r*s*(1.0+zeta);
  dN[13][2] = 9.0/4.0*r*s*(3.0*s-1.0);
  dN[14][0] = -9.0/4.0*(3.0*s-1.0)*s*(1.0+zeta);
  dN[14][1] = -9.0/4.0*(3.0*s-1.0)*s*(1.0+zeta)+9.0/4.0*(1.0-r-s)*(3.0*s-1.0)*(1.0+zeta)+27.0/4.0*(1.0-r-s)*s*(1.0+zeta);
  dN[14][2] = 9.0/4.0*(1.0-r-s)*s*(3.0*s-1.0);
  dN[15][0] = -9.0/4.0*s*(2.0-3.0*r-3.0*s)*(1.0+zeta)-27.0/4.0*(1.0-r-s)*s*(1.0+zeta);
  dN[15][1] = -9.0/4.0*s*(2.0-3.0*r-3.0*s)*(1.0+zeta)+9.0/4.0*(2.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0+zeta)-27.0/4.0*(1.0-r-s)*s*(1.0+zeta);
  dN[15][2] = 9.0/4.0*(1.0-r-s)*s*(2.0-3.0*r-3.0*s);
  dN[16][0] = 9.0/4.0*(2.0-3.0*r-3.0*s)*(1.0-r-s)*(1.0+zeta)-9.0/4.0*r*(2.0-3.0*r-3.0*s)*(1.0+zeta)-27.0/4.0*r*(1.0-r-s)*(1.0+zeta);
  dN[16][1] = -9.0/4.0*r*(2.0-3.0*r-3.0*s)*(1.0+zeta)-27.0/4.0*r*(1.0-r-s)*(1.0+zeta);
  dN[16][2] = 9.0/4.0*r*(1.0-r-s)*(2.0-3.0*r-3.0*s);
  dN[17][0] = 9.0/4.0*(1.0-r-s)*(3.0*r-1.0)*(1.0+zeta)-9.0/4.0*(3.0*r-1.0)*r*(1.0+zeta)+27.0/4.0*r*(1.0-r-s)*(1.0+zeta);
  dN[17][1] = -9.0/4.0*(3.0*r-1.0)*r*(1.0+zeta);
  dN[17][2] = 9.0/4.0*r*(1.0-r-s)*(3.0*r-1.0);
  dN[18][0] = 27.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0);
  dN[18][1] = 0.0;
  dN[18][2] = 27.0/16.0*(zeta-1.0/3.0)*(zeta-1.0)*r+27.0/16.0*(1.0+zeta)*(zeta-1.0)*r+27.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*r;
  dN[19][0] = -27.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[19][1] = 0.0;
  dN[19][2] = -27.0/16.0*(zeta+1.0/3.0)*(zeta-1.0)*r-27.0/16.0*(1.0+zeta)*(zeta-1.0)*r-27.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*r;
  dN[20][0] = 0.0;
  dN[20][1] = 27.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0);
  dN[20][2] = 27.0/16.0*(zeta-1.0/3.0)*(zeta-1.0)*s+27.0/16.0*(1.0+zeta)*(zeta-1.0)*s+27.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*s;
  dN[21][0] = 0.0;
  dN[21][1] = -27.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[21][2] = -27.0/16.0*(zeta+1.0/3.0)*(zeta-1.0)*s-27.0/16.0*(1.0+zeta)*(zeta-1.0)*s-27.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*s;
  dN[22][0] = -27.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0);
  dN[22][1] = -27.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(zeta-1.0);
  dN[22][2] = 27.0/16.0*(zeta-1.0/3.0)*(zeta-1.0)*(1.0-r-s)+27.0/16.0*(1.0+zeta)*(zeta-1.0)*(1.0-r-s)+27.0/16.0*(1.0+zeta)*(zeta-1.0/3.0)*(1.0-r-s);
  dN[23][0] = 27.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[23][1] = 27.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(zeta-1.0);
  dN[23][2] = -27.0/16.0*(zeta+1.0/3.0)*(zeta-1.0)*(1.0-r-s)-27.0/16.0*(1.0+zeta)*(zeta-1.0)*(1.0-r-s)-27.0/16.0*(1.0+zeta)*(zeta+1.0/3.0)*(1.0-r-s);
  dN[24][0] = 27.0/2.0*(1.0-r-s)*s*(1.0-zeta)-27.0/2.0*r*s*(1.0-zeta);
  dN[24][1] = 27.0/2.0*r*(1.0-r-s)*(1.0-zeta)-27.0/2.0*r*s*(1.0-zeta);
  dN[24][2] = -27.0/2.0*r*s*(1.0-r-s);
  dN[25][0] = 27.0/2.0*(1.0-r-s)*s*(1.0+zeta)-27.0/2.0*r*s*(1.0+zeta);
  dN[25][1] = 27.0/2.0*r*(1.0-r-s)*(1.0+zeta)-27.0/2.0*r*s*(1.0+zeta);
  dN[25][2] = 27.0/2.0*r*s*(1.0-r-s);
  return(26);
}

#else // "a priori" same shape fcts as above but directly in correct order ... BUGGY ... NEED SOME TESTING ...

// z in [-1,1]
void myPenta26ShapeFct(double Shape[26], double m[3])
{
  fprintf(stderr," *** DO NOT USE ROUTINE myPenta26ShapeFct: BUGGY. ABORT.\n"); exit(-1);
  double r = m[0];
  double s = m[1];
  double z = m[2];

  double t = 1.-r-s;
  
  double zzm = 1.-z*z;
  double zm  = 1.-z;
  double zp  = 1.+z;
  double z3m = 1.-3.*z;
  double z3p = 1.+3.*z;
  
  double c0  = 2./9.;
  double c1  = 9./16.;
  double c2  = 4.*c1; // = 9./4.
  double c3  = 27./2.;
  
  // Lower & upper wegde's corners  
  Shape[ 0] = c1*zm*t*( 4.*(t*(t-1.)+c0) - zzm );
  Shape[ 1] = c1*zm*r*( 4.*(r*(r-1.)+c0) - zzm );
  Shape[ 2] = c1*zm*s*( 4.*(s*(s-1.)+c0) - zzm );
  Shape[ 3] = c1*zp*t*( 4.*(t*(t-1.)+c0) - zzm );
  Shape[ 4] = c1*zp*r*( 4.*(r*(r-1.)+c0) - zzm );
  Shape[ 5] = c1*zp*s*( 4.*(s*(s-1.)+c0) - zzm );

  // Lower & upper Wegde's triangular edges'nodes
  Shape[ 6] = c2*r*t*(3.*t-1.)*zm;
  Shape[ 7] = c2*r*t*(3.*r-1.)*zm;
  Shape[ 8] = c2*r*s*(3.*r-1.)*zm;
  Shape[ 9] = c2*r*s*(3.*s-1.)*zm;
  Shape[10] = c2*t*s*(3.*s-1.)*zm;
  Shape[11] = c2*t*s*(3.*t-1.)*zm;

  Shape[12] = c2*r*t*(3.*t-1.)*zp;
  Shape[13] = c2*r*t*(3.*r-1.)*zp;
  Shape[14] = c2*r*s*(3.*r-1.)*zp;
  Shape[15] = c2*r*s*(3.*s-1.)*zp;
  Shape[16] = c2*t*s*(3.*s-1.)*zp;
  Shape[17] = c2*t*s*(3.*t-1.)*zp;

  // Wegde vertical edges'nodes
  Shape[18] = c1*t*zzm*z3m;
  Shape[19] = c1*t*zzm*z3p;
  Shape[20] = c1*r*zzm*z3m;
  Shape[21] = c1*r*zzm*z3p;
  Shape[22] = c1*s*zzm*z3m;
  Shape[23] = c1*s*zzm*z3p;

  // Lower & upper wegde triangular bubble'nodes
  Shape[24] = c3*r*s*t*zm;
  Shape[25] = c3*r*s*t*zp;
}

// z in [-1,1]
void myPenta26dShapeFct(double dShape[26][3], double m[3])
{
  fprintf(stderr," *** DO NOT USE ROUTINE myPenta26dShapeFct: BUGGY. ABORT.\n"); exit(-1);
  double r = m[0];
  double s = m[1];
  double z = m[2];

  double t = 1.-r-s;
  
  double zzm = 1.-z*z;
  double zm  = 1.-z;
  double zp  = 1.+z;
  double z3m = 1.-3.*z;
  double z3p = 1.+3.*z;
  
  double c0  = 2./9.;
  double c1  = 9./16.;
  double c2  = 4.*c1; // = 9./4.
  double c3  = 27./2.;
 
  // Lower & upper wegde's corners   
  dShape[ 0][0] = -c1*zm* ( 4.*(t*(t-1.)+c0) + 4.*t*(2.*t-1.) - zzm  );
  dShape[ 0][1] = -c1*zm* ( 4.*(t*(t-1.)+c0) + 4.*t*(2.*t-1.) - zzm  );
  dShape[ 0][2] = -c1*t * ( 4.*(t*(t-1.)+c0) - zzm - 2.*zm*z ) ;

  dShape[ 1][0] =  c1*zm* ( 4.*(r*(r-1.)+c0) + 4.*r*(2.*r-1) - zzm   );
  dShape[ 1][1] =  0.0;
  dShape[ 1][2] = -c1*r * ( 4.*(r*(r-1.)+c0) - zzm - 2.*zm*z ) ;

  dShape[ 2][0] =  0.0;
  dShape[ 2][1] =  c1*zm* ( 4.*(s*(s-1.)+c0) + 4.*s*(2.*s-1.) - zzm );
  dShape[ 2][2] = -c1*s * ( 4.*(s*(s-1.)+c0) - zzm - 2.*zm*z );

  dShape[ 3][0] = -c1*zp * ( 4.*(t*(t-1.)+c0) + 4.*t*(2.*t - 1.) - zzm );
  dShape[ 3][1] = -c1*zp * ( 4.*(t*(t-1.)+c0) + 4.*t*(2.*t - 1.) - zzm );
  dShape[ 3][2] =  c1*t  * ( 4.*(t*(t-1.)+c0) - zzm  + 2.*zp*z); 
  
  dShape[ 4][0] =  c1*zp* ( 4.*(r*(r-1.)+c0) + 4.*r*(2.*r-1) - zzm   ); 
  dShape[ 4][1] =  0.0; 
  dShape[ 4][2] =  c1*r * ( 4.*(r*(r-1.)+c0) - zzm + 2.*zp*z );
  
  dShape[ 5][0] =  0.0;
  dShape[ 5][1] =  c1*zp* ( 4.*(s*(s-1.)+c0) + 4.*s*( 2.*s -1.) - zzm );
  dShape[ 5][2] =  c1*s * ( 4.*(s*(s-1.)+c0) - zzm + 2.*zp*z );

  // Lower & upper Wegde's triangular edges'nodes  
  dShape[ 6][0] =  c2*zm*   ( t*(3.*t-1.) - r*(6.*t-1.) );
  dShape[ 6][1] = -c2*zm*r *(6.*t-1.);
  dShape[ 6][2] = -c2*r *t *(3.*t-1.);
  
  dShape[ 7][0] =  c2*zm*   (-r*(3.*r-1.) + t*(6.*r-1.) ); 
  dShape[ 7][1] = -c2*zm*r *(3.*r-1.);
  dShape[ 7][2] = -c2*r *t *(3.*r-1.);
  
  dShape[ 8][0] =  c2*zm*s*(6.*r-1.);
  dShape[ 8][1] =  c2*zm*r*(3.*r-1.);
  dShape[ 8][2] = -c2*r *s*(3.*r-1.);
   
  dShape[ 9][0] =  c2*zm*s*(3.*s-1.);
  dShape[ 9][1] =  c2*zm*r*(6.*s-1.);
  dShape[ 9][2] = -c2*r *s*(3.*s-1.);
  
  dShape[10][0] = -c2*zm*s*(3.*s-1.);
  dShape[10][1] =  c2*zm*t*(6.*s-1.);
  dShape[10][2] = -c2*t *s*(3.*s-1.);
  
  dShape[11][0] = -c2*zm*s*(6.*t-1.);
  dShape[11][1] =  c2*zm  *( t*(3.*t-1.) + s*(6.*t-1.));
  dShape[11][2] = -c2*t*s*(3.*t-1.);
  
  dShape[12][0] =  c2*zp*   ( t*(3.*t-1.) - r*(6.*t-1.) );;
  dShape[12][1] = -c2*zp*r *(6.*t-1.);
  dShape[12][2] =  c2*r *t *(3.*t-1.);
  
  dShape[13][0] =  c2*zp*   (-r*(3.*r-1.) + t*(6.*r-1.) );
  dShape[13][1] = -c2*zp*r *(3.*r-1.);
  dShape[13][2] =  c2*r *t *(3.*r-1.); 
  
  dShape[14][0] =  c2*zp*s*(6.*r-1.);
  dShape[14][1] =  c2*zp*r*(3.*r-1.);
  dShape[14][2] =  c2*r *s*(3.*r-1.);
  
  dShape[15][0] =  c2*zp*s*(3.*s-1.);
  dShape[15][1] =  c2*zp*r*(6.*s-1.);
  dShape[15][2] =  c2*r *s*(3.*s-1.);
  
  dShape[16][0] = -c2*zp*s*(3.*s-1.);
  dShape[16][1] =  c2*zp*  ( -s*(3.*s-1.) + t*(6.*s-1.));
  dShape[16][2] =  c2*t*s*(3.*s-1.);
  
  dShape[17][0] = -c2*zp*s*(6.*t-1.);
  dShape[17][1] =  c2*zp*  ( t*(3.*t-1.) - s*(6.*t-1.));
  dShape[17][2] =  c2*t *s*(3.*t-1.);

  // Wegde vertical edges'nodes    
  dShape[18][0] = -c1*zzm*z3m;
  dShape[18][1] = -c1*zzm*z3m;
  dShape[18][2] =  c1*t*(-2.*z*z3m - 3.*zzm);
  
  dShape[19][0] = -c1*zzm*z3p;
  dShape[19][1] = -c1*zzm*z3p;
  dShape[19][2] =  c1*t*(-2.*z*z3p + 3.*zzm);
  
  dShape[20][0] =  c1*zzm*z3m;
  dShape[20][1] =  0.0;
  dShape[20][2] =  c1*r*(-2.*z*z3m - 3.*zzm);
  
  dShape[21][0] =  c1*zzm*z3p;
  dShape[21][1] =  0.0;
  dShape[21][2] =  c1*r*(-2.*z*z3p + 3.*zzm);
  
  dShape[22][0] =  0.0;
  dShape[22][1] =  c1*zzm*z3m;
  dShape[22][2] =  c1*s*(-2.*z*z3m - 3.*zzm);
  
  dShape[23][0] =  0.0;
  dShape[23][1] =  c1*zzm*z3p;
  dShape[23][2] =  c1*s*(-2.*z*z3p + 3.*zzm);

  // Lower & upper wegde triangular bubble'nodes      
  dShape[24][0] =  c3*s*zm*(t-r);
  dShape[24][1] =  c3*r*zm*(t-s);
  dShape[24][2] = -c3*r*s*t;
  
  dShape[25][0] =  c3*s*zp*(t-r);
  dShape[25][1] =  c3*r*zp*(t-s);
  dShape[25][2] =  c3*r*s*t;
}
#endif

// HB (05-21-05): shape fct & derivatives for the 26 nodes wedge element (Iso-parametric formulation)
void 
Penta26ShapeFct(double Shape[26], double dShape[26][3], double m[3])
{
#ifdef USE_NAVY_PENTA26
  //fprintf(stderr," *** WARNING: use NAVY penta26 shape fcts & derivative\n");
  double shape[26];
  double dshape[26][3];

  GMap_cubicLagrangeRegion_wedge(m[0],m[1],m[2],shape);
  GMap_cubicLagrangeDrvRegion_wedge(m[0],m[1],m[2],dshape);
  
  // re-order the shape functions
  int nodeMap[26] = {2,0,1,5,3,4,10,11,6,7,8,9,16,17,12,13,14,15,22,23,18,19,20,21,24,25};
  for(int i=0; i<26; i++){
    Shape[i]     = shape[nodeMap[i]];
    dShape[i][0] = dshape[nodeMap[i]][0];
    dShape[i][1] = dshape[nodeMap[i]][1];
    dShape[i][2] = dshape[nodeMap[i]][2];
  }
#else  
  myPenta26ShapeFct(Shape,m);
  myPenta26dShapeFct(dShape,m);
#endif
}
// HB (05-21-05): shape fct & derivatives for the 26 nodes wedge element(Iso-parametric formulation)
double
Penta26ShapeFct(double Shape[26], double DShape[26][3], double m[3], double X[26], double Y[26], double Z[26])
{
  double dShape[26][3];
  Penta26ShapeFct(Shape, dShape, m);

  return(computePenta26DShapeFct(dShape, X, Y, Z, DShape));
}




