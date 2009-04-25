// ------------------------------------------------------------
// HB - 04-15-05 
// ------------------------------------------------------------

// Std C/C++ lib
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>

double
computeTetra4DShapeFct(double dShape[4][3], double X[4], double Y[4], double Z[4], double (*DShape)[3]=0)
{ 
  double xd1,xd2,xd3;
  double yd1,yd2,yd3;
  double zd1,zd2,zd3;

  xd1 = (X[1]-X[0]); xd2 = (X[2]-X[0]);  xd3 = (X[3]-X[0]);
  yd1 = (Y[1]-Y[0]); yd2 = (Y[2]-Y[0]);  yd3 = (Y[3]-Y[0]);
  zd1 = (Z[1]-Z[0]); zd2 = (Z[2]-Z[0]);  zd3 = (Z[3]-Z[0]);

  //fprintf(stderr," %e  %e  %e \n %e  %e  %e \n %e  %e  %e \n",xd1,xd2,xd3,yd1,yd2,yd3,zd1,zd2,zd3);
                                                                                                                       
  double a11 = yd2*zd3 - yd3*zd2 ; double a12 = yd3*zd1 - yd1*zd3 ; double a13 = yd1*zd2 - yd2*zd1 ;
  double a21 = xd3*zd2 - xd2*zd3 ; double a22 = xd1*zd3 - zd1*xd3 ; double a23 = xd2*zd1 - xd1*zd2 ;
  double a31 = xd2*yd3 - xd3*yd2 ; double a32 = yd1*xd3 - xd1*yd3 ; double a33 = xd1*yd2 - yd1*xd2 ;

   // -------> DETERMINANT OF THE JACOBIAN MATRIX <--- 

  double J = xd1*a11 + yd1*a21 + zd1*a31 ;  
  if ( J == 0.0 ) { 
    fprintf(stderr," *** WARNING: NULL JACOBIAN IN computeTetra4DShapeFct ROUTINE.\n"); 
  }
  //if ( J < 0.0 ) { 
  //  fprintf(stderr," *** WARNING: NEGATIVE JACOBIAN IN computeTetra4DShapeFct ROUTINE.\n"); 
  //  //J *= -1.; // BAD !!! BUT WORKS !!!
  //}
  double cdet = 1.0 / J ;

  // -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REAL" ELEMENT <--- 
 
  DShape[0][0] = cdet * ( a11*dShape[0][0] + a12*dShape[0][1] + a13*dShape[0][2] ) ;
  DShape[1][0] = cdet * ( a11*dShape[1][0] + a12*dShape[1][1] + a13*dShape[1][2] ) ;
  DShape[2][0] = cdet * ( a11*dShape[2][0] + a12*dShape[2][1] + a13*dShape[2][2] ) ;
  DShape[3][0] = cdet * ( a11*dShape[3][0] + a12*dShape[3][1] + a13*dShape[3][2] ) ;

  DShape[0][1] = cdet * ( a21*dShape[0][0] + a22*dShape[0][1] + a23*dShape[0][2] ) ;
  DShape[1][1] = cdet * ( a21*dShape[1][0] + a22*dShape[1][1] + a23*dShape[1][2] ) ;
  DShape[2][1] = cdet * ( a21*dShape[2][0] + a22*dShape[2][1] + a23*dShape[2][2] ) ;
  DShape[3][1] = cdet * ( a21*dShape[3][0] + a22*dShape[3][1] + a23*dShape[3][2] ) ;

  DShape[0][2] = cdet * ( a31*dShape[0][0] + a32*dShape[0][1] + a33*dShape[0][2] ) ;
  DShape[1][2] = cdet * ( a31*dShape[1][0] + a32*dShape[1][1] + a33*dShape[1][2] ) ;
  DShape[2][2] = cdet * ( a31*dShape[2][0] + a32*dShape[2][1] + a33*dShape[2][2] ) ;
  DShape[3][2] = cdet * ( a31*dShape[3][0] + a32*dShape[3][1] + a33*dShape[3][2] ) ; 
  
  return(J); 
}

// HB (04-15-05): shape fct & derivative for the 4 nodes tetrahedral element (Iso-parametric formulation)
// This implementation can be a bit more optimized ...
void
Tetra4ShapeFct(double Shape[4], double dShape[4][3], double m[3])
{
  double r = m[0]; // = x
  double s = m[1]; // = y
  double t = m[2]; // = z
  double u = 1.-r-s-t;

  /* -------> VALUE OF THE SHAPE FCTS IN THE "REFERENCE" ELEMENT <--- */
  
  Shape[0] = u;
  Shape[1] = r;
  Shape[2] = s;
  Shape[3] = t;

  /* -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REFERENCE" ELEMENT <--- */
  
  dShape[0][0] = -1.0 ; dShape[0][1] = -1.0 ; dShape[0][2] = -1.0;
  dShape[1][0] =  1.0 ; dShape[1][1] =  0.0 ; dShape[1][2] =  0.0; 
  dShape[2][0] =  0.0 ; dShape[2][1] =  1.0 ; dShape[2][2] =  0.0;
  dShape[3][0] =  0.0 ; dShape[3][1] =  0.0 ; dShape[3][2] =  1.0;
}

double
Tetra4ShapeFct(double Shape[4], double DShape[4][3], double m[3], double X[4], double Y[4], double Z[4])
{
  double dShape[4][3];
  Tetra4ShapeFct(Shape, dShape, m);

  return(computeTetra4DShapeFct(dShape, X, Y, Z, DShape));
}
