#ifdef USE_EIGEN3
#ifndef _MEMBRANEELEMENTTEMPLATE_CPP_
#define _MEMBRANEELEMENTTEMPLATE_CPP_

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <Element.d/Membrane.d/MembraneElementTemplate.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

template<typename doublereal>
void
MembraneElementTemplate<doublereal>
::sands19(doublereal *_xl, doublereal *_yl, doublereal *_zl,
          doublereal e, doublereal nu, doublereal *_h, 
          doublereal *_v, doublereal *_stress, 
          bool strainFlg)
{

/*********************************************************************
*
* THIS SUBROUTINE EVALUATES THE VON MISES STRESS AT THE CENTROID
* OF THE TRIANGLE. IT COMPUTES THE VALUES AT THE TOP AND BOTTOM
* AND PICKS THE MAXIMUN AS OUTPUT
* 
*-------------------------------------------------------------------*
*  CALLED BY : MDERIV.F   
*
*  SUBROUTINES CALLED:
*                      ROTATION
*                      TRIROT
*                      MEMBRA
*                      VONMISES
*-------------------------------------------------------------------*/

//.... REAL ARRAYS

        Eigen::Map<Eigen::Matrix<doublereal,3,1> > xl(_xl), yl(_yl), zl(_zl), h(_h);
        Eigen::Map<Eigen::Matrix<doublereal,7,3> > stress(_stress);
        Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);

//.... DECALRE ALL LOCAL VARIABLES

        int nd,i,j;
        Eigen::Matrix<int,9,1> le; 
        Eigen::Matrix<doublereal,3,3> dm, r1; 
        Eigen::Matrix<doublereal,3,1> xp, yp, zp, xlp, ylp, zlp, v1n, v2n, v3n; 
        Eigen::Matrix<doublereal,18,1> dll; 
        doublereal cb,x21,y21,z21,x32,y32,z32,x13,y13,z13;
        doublereal rx,ry,rz,bx,by,bz,rlr,rlb,bpr,area,ylr,zlr;
        doublereal ycg,xcg,zcg,xlcg,ylcg,zlcg,t;
        Eigen::Matrix<doublereal,18,3> rmem; 
        doublereal rmx,rmy,rmxy,rnx,rny,rnxy,ebar;
        doublereal rmmx,rmmy,rmmxy,rnnx,rnny,rnnxy,sbf;
        Eigen::Matrix<doublereal,3,1> xg,yg,zg;
        Eigen::Matrix<doublereal,6,1> str;

//.... INIALIZE DATA

        le << 0,1,6,7,12,13,5,11,17;
        nd = 18;
        t  = h[0];
        cb = e*t/(1.0-nu*nu);

        rmem.setZero();
        dm.setZero();
        dm(0,0) =    cb;
        dm(0,1) = nu*cb;
        dm(1,0) = nu*cb;
        dm(1,1) =    cb;
        dm(2,2) = (1.0-nu)/2.0*cb;

//.... DIMENSION VARIABLES

         x21 = xl[1] - xl[0];
         y21 = yl[1] - yl[0];
         z21 = zl[1] - zl[0];
         x32 = xl[2] - xl[1];
         y32 = yl[2] - yl[1];
         z32 = zl[2] - zl[1];
         x13 = xl[0] - xl[2];
         y13 = yl[0] - yl[2];
         z13 = zl[0] - zl[2];

//....  TRIANGLE IN SPACE : WE COMPUTE THE LENGTH 
//....  OF ONE SIDE AND THE DISTANCE OF THE
//....  OPPOSING NODE TO THAT SIDE TO COMPUTE THE AREA

      rx   = x21;
      ry   = y21;
      rz   = z21;
      bx   = x32;
      by   = y32;
      bz   = z32;
      rlr  = sqrt( rx*rx + ry*ry + rz*rz );
      rlb  = sqrt( bx*bx + by*by + bz*bz );
      bpr  = abs(rx * bx + ry * by + rz *bz )/rlr;
      area = rlr*(sqrt(rlb*rlb-bpr*bpr))/2.0;

//.... DIRECTION COSINES OF THE LOCAL SYSTEM . 
//.... X' IS DIRCTED PARALLEL TO THE 2-1 SIDE
//.... Z' IS THE EXTERNAL NORMAL (COUNTERCLOCKWISE). 
//.... Y' COMPUTED AS Z'x X'

      xp[0] = x21/rlr;
      xp[1] = y21/rlr;
      xp[2] = z21/rlr;

//.... Z'

      zp[0] = y21 * z32 - z21 * y32;
      zp[1] = z21 * x32 - x21 * z32;
      zp[2] = x21 * y32 - y21 * x32;
      zlr   = sqrt( zp(0)*zp[0] + zp[1]*zp[1] + zp(2)*zp[2] );
      zp[0] = zp[0]/zlr;
      zp[1] = zp[1]/zlr;
      zp[2] = zp[2]/zlr;

//.... Y'

      yp[0] = zp[1] * xp[2] - zp[2] * xp[1];
      yp[1] = zp[2] * xp[0] - zp[0] * xp[2];
      yp[2] = zp[0] * xp[1] - zp[1] * xp[0];
      ylr   = sqrt( yp[0]*yp[0] + yp[1]*yp[1] + yp[2]*yp[2] );
      yp[0] = yp[0]/ylr;
      yp[1] = yp[1]/ylr;
      yp[2] = yp[2]/ylr;

//.... CENTER OF GRAVITY

      xcg = (xl[0] + xl[1] + xl[2])/3.0;
      ycg = (yl[0] + yl[1] + yl[2])/3.0;
      zcg = (zl[0] + zl[1] + zl[2])/3.0;

//.... COMPUTING LOCAL COORDINATES 

      for(i=0; i<3; ++i) {
        xlcg   = xl(i)-xcg;
        ylcg   = yl(i)-ycg;
        zlcg   = zl(i)-zcg;
        xlp[i] = xp[0] * xlcg + xp[1] * ylcg + xp[2] * zlcg;
        ylp[i] = yp[0] * xlcg + yp[1] * ylcg + yp[2] * zlcg;
        zlp[i] = zp[0] * xlcg + zp[1] * ylcg + zp[2] * zlcg;
      }

//.... COMPUTING NODAL ROTATION MATRICES

        v1n[0] = 1.0;
        v1n[1] = 0.0;
        v1n[2] = 0.0;
        v2n[0] = 0.0;
        v2n[1] = 1.0;
        v2n[2] = 0.0;
        v3n[0] = 0.0;
        v3n[1] = 0.0;
        v3n[2] = 1.0;

        rotation(xp.data(),yp.data(),zp.data(),v1n.data(),v2n.data(),v3n.data(),r1.data());

//.... COMPUTE THE VON MISES STRESS

//.... ROTATE NODAL DISPLACEMENTS TO LOCAL SYSTEM

        for(i=0; i<18; ++i) {
          dll[i]= 0.0;
        }

        for(i=0; i<3; ++i) {
          for(j=0; j<3; ++j) {
            dll[i]   =dll[i]   +r1(j,i)*v[j];
            dll[i+3] =dll[i+3] +r1(j,i)*v[j+3];
            dll[i+6] =dll[i+6] +r1(j,i)*v[j+6];
            dll[i+9] =dll[i+9] +r1(j,i)*v[j+9];
            dll[i+12]=dll[i+12]+r1(j,i)*v[j+12];
            dll[i+15]=dll[i+15]+r1(j,i)*v[j+15];
          }
        }

//.... COMPUTE CENTROIDAL MEMBRANE STRAINS

        membra(xlp.data(),ylp.data(),1.5,le.data(),rmem.data());

        rmx  = 0.0;
        rmy  = 0.0;
        rmxy = 0.0;
 
        rnx  = 0.0;
        rny  = 0.0;
        rnxy = 0.0;

        for(j=0; j<18; ++j) {
          rnx  = rnx  + rmem(j,0)*dll[j];
          rny  = rny  + rmem(j,1)*dll[j];
          rnxy = rnxy + rmem(j,2)*dll[j];
        }

// ...  COMPUTE EQUIVALENT STRAINS

       xg[0] = 1.0;
       xg[1] = 0.0;
       xg[2] = 0.0;
       yg[0] = 0.0;
       yg[1] = 1.0;
       yg[2] = 0.0;
       zg[0] = 0.0;
       zg[1] = 0.0;
       zg[2] = 1.0;
       if(strainFlg == true) {
         str(0) = rnx;
         str(1) = rny;
         str(2) = 0.0;
         str(3) = 0.5*rnxy;
         str(4) = -(nu/(1.0-nu))*(rnx + rny);
         str(5) = 0.0;
         transform(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),str.data());
         stress(0,0) = str[0];
         stress(0,1) = str[0];
         stress(0,2) = str[0];
         stress(1,0) = str[1];
         stress(1,1) = str[1];
         stress(1,2) = str[1];
         stress(2,0) = str[2];
         stress(2,1) = str[2];
         stress(2,2) = str[2];
         stress(3,0) = 2.0*str[3];
         stress(3,1) = 2.0*str[3];
         stress(3,2) = 2.0*str[3];
         stress(4,0) = 2.0*str[4];
         stress(4,1) = 2.0*str[4];
         stress(4,2) = 2.0*str[4];
         stress(5,0) = 2.0*str[5];
         stress(5,1) = 2.0*str[5];
         stress(5,2) = 2.0*str[5];

         ebar= ((rnx-rny)*(rnx-rny)+rnx*rnx+rny*rny)/6.0;
         ebar = ebar + 0.25*(rnxy*rnxy);
         ebar = sqrt(3.0*ebar);

         stress(6,0) = ebar;
         stress(6,1) = ebar;
         stress(6,2) = ebar;
         return;
       }

//....  COMPUTE CENTROIDAL STRESS RESULTANTS

       rmmx  = 0.0;
       rmmy  = 0.0;
       rmmxy = 0.0; 

       rnnx  = dm(0,0)*rnx+dm(0,1)*rny+dm(0,2)*rnxy;
       rnny  = dm(1,0)*rnx+dm(1,1)*rny+dm(1,2)*rnxy;
       rnnxy = dm(2,0)*rnx+dm(2,1)*rny+dm(2,2)*rnxy;

// ...  COMPUTE VON MISES STRESS RESULTANT

       vonmis(rmmx,rmmy,rmmxy,rnnx,rnny,rnnxy,t,sbf);

// ... COMPUTE SIGMAXX SIGMAYY AND SIGMAXY

      str[0] = rnnx/t;
      str[1] = rnny/t;
      str[2] = 0.0;
      str[3] = rnnxy/t;
      str[4] = 0.0;
      str[5] = 0.0;

      transform(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),str.data());

      stress(0,0) = str[0];
      stress(0,1) = str[0];
      stress(0,2) = str[0];

      stress(1,0) = str[1];
      stress(1,1) = str[1];
      stress(1,2) = str[1];

      stress(2,0) = str[2];
      stress(2,1) = str[2];
      stress(2,2) = str[2];

      stress(3,0) = str[3];
      stress(3,1) = str[3];
      stress(3,2) = str[3];

      stress(4,0) = str[4];
      stress(4,1) = str[4];
      stress(4,2) = str[4];

      stress(5,0) = str[5];
      stress(5,1) = str[5];
      stress(5,2) = str[5];

      stress(6,0) = sbf; 
      stress(6,1) = sbf;
      stress(6,2) = sbf;
}

template<typename doublereal>
void
MembraneElementTemplate<doublereal>
::rotation(doublereal *_v1o, doublereal *_v2o, doublereal *_v3o,
           doublereal *_v1n, doublereal *_v2n, doublereal *_v3n,
           doublereal *_r)
{
// Declarations
      Eigen::Map<Eigen::Matrix<doublereal,3,1> > v1o(_v1o), v2o(_v2o), v3o(_v3o), v1n(_v1n), v2n(_v2n), v3n(_v3n);
      Eigen::Map<Eigen::Matrix<doublereal,3,3> > r(_r);

// Local Declarations
      int i,j;

// Zero rotation matrix
      r.setZero();

// Compute rotation matrix
      for(i=0; i<3; ++i) {
        r(0,0) = r(0,0) + v1n[i]*v1o[i];
        r(0,1) = r(0,1) + v1n[i]*v2o[i];
        r(0,2) = r(0,2) + v1n[i]*v3o[i];
        r(1,0) = r(1,0) + v2n[i]*v1o[i];
        r(1,1) = r(1,1) + v2n[i]*v2o[i];
        r(1,2) = r(1,2) + v2n[i]*v3o[i];
        r(2,0) = r(2,0) + v3n[i]*v1o[i];
        r(2,1) = r(2,1) + v3n[i]*v2o[i];
        r(2,2) = r(2,2) + v3n[i]*v3o[i];
      }
}

template<typename doublereal>
void
MembraneElementTemplate<doublereal>
::membra(doublereal *_x, doublereal *_y, doublereal alpha,
         int *_le, doublereal *_q)
{
      Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y);
      Eigen::Map<Eigen::Matrix<doublereal,18,3> > q(_q);
      Eigen::Map<Eigen::Matrix<int,9,1> > le(_le);

      Eigen::Matrix<doublereal,9,3> p;
      doublereal area2, c, x21, x32, x13, y21, y32, y13, x12, x23, x31, y12, y23, y31; 

      int i, j, kk, n;

      x21 =      x[1] - x[0];
      x12 =     -x21;
      x32 =      x[2] - x[1];
      x23 =     -x32;
      x13 =      x[0] - x[2];
      x31 =     -x13;
      y21 =      y[1] - y[0];
      y12 =     -y21;
      y32 =      y[2] - y[1];
      y23 =     -y32;
      y13 =      y[0] - y[2];
      y31 =     -y13;
      area2 =    y21*x13 - x21*y13;
      if (area2 <= 0.0) {
        cerr << " ... Error! SM3MB: Zero area\n";
        if (area2 == 0.0)   cerr << " ... Error! SM3MB: Zero area\n";
        exit(-1);
      }
      p(0,0) =   y23;
      p(1,0) =   0.0;
      p(2,0) =   y31;
      p(3,0) =   0.0;
      p(4,0) =   y12;
      p(5,0) =   0.0;
      p(0,1) =   0.0;
      p(1,1) =   x32;
      p(2,1) =   0.0;
      p(3,1) =   x13;
      p(4,1) =   0.0;
      p(5,1) =   x21;
      p(0,2) =   x32;
      p(1,2) =   y23;
      p(2,2) =   x13;
      p(3,2) =   y31;
      p(4,2) =   x21;
      p(5,2) =   y12;
      if (alpha != 0.0) {
        p(6,0) =  y23*(y13-y21)*alpha/6.;
        p(6,1) =  x32*(x31-x12)*alpha/6.;
        p(6,2) =  (x31*y13-x12*y21)*alpha/3.;
        p(7,0) =  y31*(y21-y32)*alpha/6.;
        p(7,1) =  x13*(x12-x23)*alpha/6.;
        p(7,2) =  (x12*y21-x23*y32)*alpha/3.;
        p(8,0) =  y12*(y32-y13)*alpha/6.;
        p(8,1) =  x21*(x23-x31)*alpha/6.;
        p(8,2) =  (x23*y32-x31*y13)*alpha/3.;
        n =        9;
      }
      c =   1.0   /area2;
      for(i=0; i<9; ++i) {
        kk=le[i];
        for(j=0; j<3; ++j) {
          q(kk,j)=p(i,j)*c;
        }
      }
}

template<typename doublereal>
void
MembraneElementTemplate<doublereal>
::transform(doublereal *_xl, doublereal *_yl, doublereal *_zl,
            doublereal *_xg, doublereal *_yg, doublereal *_zg,
            doublereal *_str)
{

/**********************************************************************C
C
C Purpose: to form a transformation matrix from local coordinates to
C          global coordinates
C
C input variables:
C      xl = x local unit vector
C      yl = y local unit vector
C      zl = z local unit vector
C      xg = x global unit vector
C      yg = y global unit vector
C      zg = z global unit vector
C      str = stress/strain 6x1 vector
C            sigmaxx, sigmayy, sigmazz, sigma12, sigma23, sigma13
C
C local variables:
C      l1 = direction cosine between xl and xg
C      l2 = direction cosine between xl and yg
C      l3 = direction cosine between xl and zg
C      m1 = direction cosine between yl and xg
C      m2 = direction cosine between yl and yg
C      m3 = direction cosine between yl and zg
C      n1 = direction cosine between zl and xg
C      n2 = direction cosine between zl and yg
C      n3 = direction cosine between zl and zg
C      t  = transformation matrix from local to global
C
C**********************************************************************/

        Eigen::Map<Eigen::Matrix<doublereal,3,1> > xl(_xl),yl(_yl),zl(_zl), xg(_xg), yg(_yg), zg(_zg);
        Eigen::Map<Eigen::Matrix<doublereal,6,1> > str(_str);

// Local Declarations

        doublereal l1,l2,l3;
        doublereal m1,m2,m3;
        doublereal n1,n2,n3;
        Eigen::Matrix<doublereal,6,1> s;
        Eigen::Matrix<doublereal,6,6> t;

// Copy stress/strain values to a temporary array

        s = str;

// Compute direction cosines
      
        l1 = xg.dot(xl);
        l2 = yg.dot(xl);
        l3 = zg.dot(xl);
        
        m1 = xg.dot(yl);
        m2 = yg.dot(yl);
        m3 = zg.dot(yl);
        
        n1 = xg.dot(zl);
        n2 = yg.dot(zl);
        n3 = zg.dot(zl);

// Construct the 6x6 transformation matrix
       
        t(0,0) = l1*l1;
        t(0,1) = m1*m1;
        t(0,2) = n1*n1;
        t(0,3) = 2.0*l1*m1;
        t(0,4) = 2.0*m1*n1;
        t(0,5) = 2.0*n1*l1;
       
        t(1,0) = l2*l2;
        t(1,1) = m2*m2;
        t(1,2) = n2*n2;
        t(1,3) = 2.0*l2*m2;
        t(1,4) = 2.0*m2*n2;
        t(1,5) = 2.0*n2*l2;
       
        t(2,0) = l3*l3;
        t(2,1) = m3*m3;
        t(2,2) = n3*n3;
        t(2,3) = 2.0*l3*m3;
        t(2,4) = 2.0*m3*n3;
        t(2,5) = 2.0*n3*l3;
       
        t(3,0) = l1*l2;
        t(3,1) = m1*m2;
        t(3,2) = n1*n2;
        t(3,3) = l1*m2 + l2*m1;
        t(3,4) = m1*n2 + m2*n1;
        t(3,5) = n1*l2 + n2*l1;
       
        t(4,0) = l2*l3;
        t(4,1) = m2*m3;
        t(4,2) = n2*n3;
        t(4,3) = l2*m3 + l3*m2;
        t(4,4) = m2*n3 + m3*n2;
        t(4,5) = n2*l3 + n3*l2;
       
        t(5,0) = l3*l1;
        t(5,1) = m3*m1;
        t(5,2) = n3*n1;
        t(5,3) = l3*m1 + l1*m3;
        t(5,4) = m3*n1 + m1*n3;
        t(5,5) = n3*l1 + n1*l3;

// Perform the multiplication {str'} = T{str}
       
        str = t*s;
} 

template<typename doublereal>
void
MembraneElementTemplate<doublereal>
::vonmis(doublereal rmx, doublereal rmy, doublereal rmxy,
         doublereal rnx, doublereal rny, doublereal rnxy,
         doublereal &t,   doublereal &sv)
{
//.... LOCAL VARIABLES
// st = von mises stress in top surface
// sm = von mises stress in median surface
// sb = von mises stress in bottom surface

      doublereal sx,sy,sxy,st,sb,sm,t2,sq3;
      doublereal rnxt,rnyt,rnxyt,rmxt,rmyt,rmxyt;

      t    = abs(t);
      t2   = t*t;

      sq3  = sqrt(3.0);

      rnxt  =  rnx/t;
      rnyt  =  rny/t;
      rnxyt =  rnxy/t;

      rmxt  = 6.0 * rmx / t2;
      rmyt  = 6.0 * rmy / t2;
      rmxyt = 6.0 *rmxy / t2;

// ... COMPUTE VON MISES STRESS IN BOTTOM SURFACE
      sx  =  rnxt -  rmxt;
      sy  =  rnyt -  rmyt;
      sxy = rnxyt - rmxyt;
      compj2(sx,sy,sxy,sb);
      sb = sq3 * sb;

// ... COMPUTE VON MISES STRESS IN MEDIAN SURFACE
      sx  = rnxt;
      sy  = rnyt;
      sxy = rnxyt;
      compj2(sx,sy,sxy,sm);
      sm = sq3 * sm;

// ... COMPUTE VON MISES STRESS IN TOP SURFACE
      sx  =  rnxt +  rmxt;
      sy  =  rnyt +  rmyt;
      sxy = rnxyt + rmxyt;
      compj2(sx,sy,sxy,st);
      st = sq3 *st;

// ... VON MISES STRESS = max(sb,st)
      sv = max(sb,st);
}

template<typename doublereal>
void
MembraneElementTemplate<doublereal>
::compj2(doublereal sx, doublereal sy, doublereal sxy, doublereal &svm)
{

// ... SUBROUTINE TO CALCULATE J2
      doublereal sz,s0,dsx,dsy,dsz,j2;

// ... SET sz = 0 TO REMIND USER OF THIS ASSUMPTION
      sz = 0.0;

// ... COMPUTE AVERAGE HYDROSTATIC STRESS
      s0 = (sx + sy + sz)/3.0;

// ... COMPUTE DEVIATORIC STRESSES
      dsx  = sx - s0;
      dsy  = sy - s0;
      dsz  = sz - s0;

// ... COMPUTE J2

      j2 = 0.5*((dsx*dsx) + (dsy*dsy) + (dsz*dsz)) + (sxy*sxy);

// ... COMPUTE VON MISES STRESS
      svm = sqrt(j2);
}

#endif
#endif
