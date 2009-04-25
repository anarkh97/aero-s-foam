// ---------------------------------------------------------------------
// HB - 01-15-06
// ---------------------------------------------------------------------
// Shape fcts & derivatives for 15 nodes wedge element 
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


// HB (01-15-06): for given reference shape fct derivatives dShape & nodes coordinates (X.Y,Z), 
//                compute the jacobian & shape fcts derivative w.r.t to the global coordinate system.
//                assuming an iso-parametric formulation.
double
computePenta15DShapeFct(double dShape[15][3], double X[15], double Y[15], double Z[15], double (*DShape)[3]=0)
{
  double xd1,xd2,xd3;
  double yd1,yd2,yd3;
  double zd1,zd2,zd3;

  xd1 = 0.0 ; xd2 = 0.0 ; xd3 = 0.0 ;
  yd1 = 0.0 ; yd2 = 0.0 ; yd3 = 0.0 ;
  zd1 = 0.0 ; zd2 = 0.0 ; zd3 = 0.0 ;

  for(int i=0; i<15; i+=3){
    xd1 += dShape[i][0]*X[i] + dShape[i+1][0]*X[i+1] + dShape[i+2][0]*X[i+2];
    xd2 += dShape[i][1]*X[i] + dShape[i+1][1]*X[i+1] + dShape[i+2][1]*X[i+2];
    xd3 += dShape[i][2]*X[i] + dShape[i+1][2]*X[i+1] + dShape[i+2][2]*X[i+2];
    
    yd1 += dShape[i][0]*Y[i] + dShape[i+1][0]*Y[i+1] + dShape[i+2][0]*Y[i+2];
    yd2 += dShape[i][1]*Y[i] + dShape[i+1][1]*Y[i+1] + dShape[i+2][1]*Y[i+2];
    yd3 += dShape[i][2]*Y[i] + dShape[i+1][2]*Y[i+1] + dShape[i+2][2]*Y[i+2];
    
    zd1 += dShape[i][0]*Z[i] + dShape[i+1][0]*Z[i+1] + dShape[i+2][0]*Z[i+2];
    zd2 += dShape[i][1]*Z[i] + dShape[i+1][1]*Z[i+1] + dShape[i+2][1]*Z[i+2];
    zd3 += dShape[i][2]*Z[i] + dShape[i+1][2]*Z[i+1] + dShape[i+2][2]*Z[i+2];
  }
  //fprintf(stderr," %e  %e  %e \n %e  %e  %e \n %e  %e  %e \n",xd1,xd2,xd3,yd1,yd2,yd3,zd1,zd2,zd3);
  double a11 = yd2*zd3 - yd3*zd2 ; double a12 = yd3*zd1 - yd1*zd3 ; double a13 = yd1*zd2 - yd2*zd1 ;
  double a21 = xd3*zd2 - xd2*zd3 ; double a22 = xd1*zd3 - zd1*xd3 ; double a23 = xd2*zd1 - xd1*zd2 ;
  double a31 = xd2*yd3 - xd3*yd2 ; double a32 = yd1*xd3 - xd1*yd3 ; double a33 = xd1*yd2 - yd1*xd2 ;

  /* -------> DETERMINANT OF THE JACOBIAN MATRIX <--- */
  double J = xd1*a11 + yd1*a21 + zd1*a31 ;
  if ( J == 0.0 ) {
    fprintf(stderr," *** WARNING: NULL JACOBIAN IN computePenta15DShapeFct.C ROUTINE.\n");
  }
  if ( J < 0.0 ) {
    fprintf(stderr," *** WARNING: NEGATIVE JACOBIAN IN computePenta15DShapeFct.C ROUTINE.\n");
    //J *= -1.; // BAD !!! BUT WORKS !!!
  }
  double cdet = 1.0 / J ;

  /* -------> DERIVATIVE OF THE SHAPE FCTS IN THE "REAL" ELEMENT <--- */
  if(DShape){
    for(int i=0; i<15; i+=3){
      DShape[i  ][0] = cdet * ( a11*dShape[i  ][0] + a12*dShape[i  ][1] + a13*dShape[i  ][2] );
      DShape[i  ][1] = cdet * ( a21*dShape[i  ][0] + a22*dShape[i  ][1] + a23*dShape[i  ][2] );
      DShape[i  ][2] = cdet * ( a31*dShape[i  ][0] + a32*dShape[i  ][1] + a33*dShape[i  ][2] );
     
      DShape[i+1][0] = cdet * ( a11*dShape[i+1][0] + a12*dShape[i+1][1] + a13*dShape[i+1][2] );
      DShape[i+1][1] = cdet * ( a21*dShape[i+1][0] + a22*dShape[i+1][1] + a23*dShape[i+1][2] );
      DShape[i+1][2] = cdet * ( a31*dShape[i+1][0] + a32*dShape[i+1][1] + a33*dShape[i+1][2] );

      DShape[i+2][0] = cdet * ( a11*dShape[i+2][0] + a12*dShape[i+2][1] + a13*dShape[i+2][2] );
      DShape[i+2][1] = cdet * ( a21*dShape[i+2][0] + a22*dShape[i+2][1] + a23*dShape[i+2][2] );
      DShape[i+2][2] = cdet * ( a31*dShape[i+2][0] + a32*dShape[i+2][1] + a33*dShape[i+2][2] );
    } 
  }
  
  return(J);
}

// z in [-1,1]
void Penta15ShapeFct(double Shape[15], double m[3])
{
  double r = m[0];
  double s = m[1];
  double z = m[2];

  double t = 1.-r-s;
  
  double zm  = 0.5*(1.-z);
  double zp  = 0.5*(1.+z);
  double zzm = 1.-z*z;
 
  // Lower & upper wegde's corners  
  Shape[ 0] = zm*t*( 2.*t-z-2. );
  Shape[ 1] = zm*r*( 2.*r-z-2. );
  Shape[ 2] = zm*s*( 2.*s-z-2. );
  Shape[ 3] = zp*t*( 2.*t+z-2. );
  Shape[ 4] = zp*r*( 2.*r+z-2. );
  Shape[ 5] = zp*s*( 2.*s+z-2. );

  // Lower & upper Wegde's triangular edges'nodes
  Shape[ 6] = 4.*zm*t*r;
  Shape[ 7] = 4.*zm*r*s;
  Shape[ 8] = 4.*zm*s*t;
  Shape[ 9] = 4.*zp*t*r;
  Shape[10] = 4.*zp*r*s;
  Shape[11] = 4.*zp*s*t;

  // Wegde vertical edges'nodes
  Shape[12] = zzm*t;
  Shape[13] = zzm*r;
  Shape[14] = zzm*s;
}

// z in [-1,1]
void Penta15dShapeFct(double dShape[15][3], double m[3])
{
  double r = m[0];
  double s = m[1];
  double z = m[2];

  double t = 1.-r-s;

  double zm  = 0.5*(1.-z);
  double zp  = 0.5*(1.+z);
  double zzm = 1.-z*z;

  // Lower & upper wegde's corners
  dShape[ 0][0] = zm*( 2.+z-4.*t ); 
  dShape[ 0][1] = zm*( 2.+z-4.*t );
  dShape[ 0][2] =  t*(  z-t+0.5  );

  dShape[ 1][0] = zm*( 4.*r-z-2. ); 
  dShape[ 1][1] = zm*( 4.*r-z-2. );
  dShape[ 1][2] =  r*(  z-r+0.5  );

  dShape[ 2][0] = zm*( 4.*s-z-2. ); 
  dShape[ 2][1] = zm*( 4.*s-z-2. );
  dShape[ 2][2] =  s*(  z-s+0.5  );

  dShape[ 3][0] = zp*( 2.-z-4.*t );
  dShape[ 3][1] = zp*( 2.-z-4.*t );
  dShape[ 3][2] =  t*(  z+t-0.5  );

  dShape[ 4][0] = zp*( 4.*r+z-2. );
  dShape[ 4][1] = zp*( 4.*r+z-2. );
  dShape[ 4][2] =  r*(  z+r-0.5  );

  dShape[ 5][0] = zp*( 4.*s+z-2. );
  dShape[ 5][1] = zp*( 4.*s+z-2. );
  dShape[ 5][2] =  s*(  z+s-0.5  );

  // Lower & upper Wegde's triangular edges'nodes
  dShape[ 6][0] =  4.*zm*(t-r);
  dShape[ 6][1] = -4.*zm*r    ;
  dShape[ 6][2] = -2.* t*r    ;

  dShape[ 7][0] =  4.*zm*s;
  dShape[ 7][1] =  4.*zm*r;
  dShape[ 7][2] = -2.* r*s;

  dShape[ 8][0] = -4.*zm*s;
  dShape[ 8][1] =  4.*zm*(t-s);
  dShape[ 8][2] = -2.* s*t;

  dShape[ 9][0] =  4.*zp*(t-r);
  dShape[ 9][1] = -4.*zp*r    ;
  dShape[ 9][2] =  2.* t*r    ;

  dShape[10][0] =  4.*zp*s;
  dShape[10][1] =  4.*zp*r;
  dShape[10][2] =  2.* r*s;

  dShape[11][0] = -4.*zp*s;
  dShape[11][1] =  4.*zp*(t-s);
  dShape[11][2] =  2.* s*t;

  // Wegde vertical edges'nodes
  dShape[12][0] = -zzm;
  dShape[12][1] = -zzm;
  dShape[12][2] = -2.*z*t;

  dShape[13][0] =  zzm;
  dShape[13][1] =  0.0;
  dShape[13][2] = -2.*z*r;

  dShape[14][0] =  0.0;
  dShape[14][1] =  zzm;
  dShape[14][2] = -2.*z*s;
}

// HB (01-15-06): shape fct & derivatives for the 15 nodes wedge element (Iso-parametric formulation)
void 
Penta15ShapeFct(double Shape[15], double dShape[15][3], double m[3])
{
  Penta15ShapeFct(Shape,m);
  Penta15dShapeFct(dShape,m);
}

// HB (01-15-06): shape fct & derivatives for the 15 nodes wedge element(Iso-parametric formulation)
double
Penta15ShapeFct(double Shape[15], double DShape[15][3], double m[3], double X[15], double Y[15], double Z[15])
{
  double dShape[15][3];
  Penta15ShapeFct(Shape, dShape, m);

  return(computePenta15DShapeFct(dShape, X, Y, Z, DShape));
}




