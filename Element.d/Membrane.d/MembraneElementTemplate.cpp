#ifdef USE_EIGEN3
#ifndef _MEMBRANEELEMENTTEMPLATE_CPP_
#define _MEMBRANEELEMENTTEMPLATE_CPP_

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <Element.d/Membrane.d/MembraneElementTemplate.hpp>
#include <Element.d/Shell.d/ShellElementSemiTemplate.cpp>
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
        Eigen::Matrix<doublereal,6,6> trsM;

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

        this->rotation(xp.data(),yp.data(),zp.data(),v1n.data(),v2n.data(),v3n.data(),r1.data());

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

        this->membra(xlp.data(),ylp.data(),1.5,le.data(),rmem.data());

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
         this->transform(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),str.data(),trsM.data());
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

       this->vonmis(rmmx,rmmy,rmmxy,rnnx,rnny,rnnxy,t,sbf);

// ... COMPUTE SIGMAXX SIGMAYY AND SIGMAXY

      str[0] = rnnx/t;
      str[1] = rnny/t;
      str[2] = 0.0;
      str[3] = rnnxy/t;
      str[4] = 0.0;
      str[5] = 0.0;

      this->transform(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),str.data(),trsM.data());

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
::trimem(int flag, doublereal *_xl, doublereal *_yl, doublereal *_zl,
         doublereal e, doublereal nu, doublereal *_h, doublereal *_rk) 
{
        Eigen::Map<Eigen::Matrix<doublereal,18,18> > rk(_rk);
        Eigen::Map<Eigen::Matrix<doublereal,3,1> > xl(_xl), yl(_yl), zl(_zl), h(_h);

        Eigen::Matrix<doublereal,9,9> sm;
        Eigen::Matrix<doublereal,3,3> dm, r1;
        Eigen::Matrix<doublereal,3,1> xp, yp, zp, xlp, ylp, zlp, v1n, v2n, v3n;

        doublereal cb,x21,y21,z21,x32,y32,z32,x13,z13;
        doublereal rx,ry,rz,bx,by,bz,rlr,rlb,bpr,area,ylr,zlr;
        doublereal ycg,xcg,zcg,xlcg,ylcg,zlcg,esp;
        Eigen::Matrix<int,9,1> le, ptr;

        int i,j;
        char* status = " ";
        le << 0,1,6,7,12,13,5,11,17;
        ptr << 0,1,5,6,7,11,12,13,17;
        v1n << 1.0,0.0,0.0;
        v2n << 0.0,1.0,0.0;
        v3n << 0.0,0.0,1.0;

//.... thickness

        esp = h[0];

//.... elastic matrix

        cb=(e*esp)/(1.0-(nu*nu));
        dm(0,0) = cb;
        dm(0,1) = nu*cb;
        dm(0,2) = 0.0;
        dm(1,0) = dm(0,1);
        dm(1,1) = cb;
        dm(1,2) = 0.0;
        dm(2,0) = 0.0;
        dm(2,1) = 0.0;
        dm(2,2) = ((1.0-nu)/2.0)*cb;

//.... dimension variables

        x21 = xl[1] - xl[0];
        y21 = yl[1] - yl[0];
        z21 = zl[1] - zl[0];
        x32 = xl[2] - xl[1];
        y32 = yl[2] - yl[1];
        z32 = zl[2] - zl[1];
        x13 = xl[0] - xl[2];
        z13 = zl[0] - zl[2];

//.... triangle in space : we compute the length of one 
//.... side and the distance of the
//.... opposing node to that side to compute the area

        rx   = x21;
        ry   = y21;
        rz   = z21;
        bx   = x32;
        by   = y32;
        bz   = z32;

        rlr  = sqrt( rx*rx + ry*ry + rz*rz );
        rlb  = sqrt( bx*bx + by*by + bz*bz );
        bpr  = sqrt((rx * bx + ry * by + rz *bz)*(rx * bx + ry * by + rz *bz))/rlr;
        area = 0.5*rlr*(sqrt(rlb*rlb-bpr*bpr));

//.... direction cosines of the local system . 
//.... X' is dircted parallel to the 2-1 side
//.... Z' is the external normal (counterclockwise). 
//.... Y' computed as Z'x X'

        xp[0] = x21/rlr;
        xp[1] = y21/rlr;
        xp[2] = z21/rlr;

//.... Z'

        zp[0] = y21 * z32 - z21 * y32;
        zp[1] = z21 * x32 - x21 * z32;
        zp[2] = x21 * y32 - y21 * x32;
        zlr   = sqrt( zp[0]*zp[0] + zp[1]*zp[1] + zp[2]*zp[2] );
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

//.... center of gravity

        xcg = (xl[0] + xl[1] + xl[2])/3.0;
        ycg = (yl[0] + yl[1] + yl[2])/3.0;
        zcg = (zl[0] + zl[1] + zl[2])/3.0;

//.... computing local coordinates 

        for(i=0; i<3; ++i) {
          xlcg   = xl[i]-xcg;
          ylcg   = yl[i]-ycg;
          zlcg   = zl[i]-zcg;
          xlp[i] = xp[0] * xlcg + xp[1] * ylcg + xp[2] * zlcg;
          ylp[i] = yp[0] * xlcg + yp[1] * ylcg + yp[2] * zlcg;
          zlp[i] = zp[0] * xlcg + zp[1] * ylcg + zp[2] * zlcg;
        }

//.... cleaning stiffness matrix
        rk.setZero();

//.... forming local basic stiffness for membrane

        this->sm3mb(xlp.data(),ylp.data(),dm.data(),1.5,1.0,le.data(),rk.data(),18,status);

//.... forming local higher order stiffness for membrane

        this->sm3mhe(xlp.data(),ylp.data(),dm.data(),0.32,le.data(),rk.data(),18,status);

// rotate stiffness matrix from local coordinate system to global
// coordinate system in the case of linear FEM. In the case of
// nonlinear FEM with corotational method, do not perform this 
// transformation as the corotational routines expect a stiffness
// matrix in local coordinates.

        if(flag == 1) {
          //.... computing nodal rotation matrices
          this->rotation(xp.data(),yp.data(),zp.data(),v1n.data(),v2n.data(),v3n.data(),r1.data());
          //.... rotate membrane stiffness
          this->trirotation(rk.data(),r1.data());
        }

}

#endif
#endif
