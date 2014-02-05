#ifdef USE_EIGEN3
#ifndef _SHELLELEMENTSEMITEMPLATE_CPP_
#define _SHELLELEMENTSEMITEMPLATE_CPP_

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <Element.d/Shell.d/ShellElementSemiTemplate.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::sands8(doublereal *_xl,doublereal *_yl, doublereal *_zl,
         doublereal e, doublereal nu, doublereal *_h,
         doublereal *_v, doublereal *_stress,
         bool strainFlg, int surface, doublereal thrmStr)
{
/*-----------------------------------------------------------------
C This subroutine evaluates the von mises stress at the centroid
C of the triangle. it computes the values at the top and bottom
C and picks the maximum as output
C for the spacial 18 d.o.f tree node triangle obtained as
C a combination of the aqr bending triangle plus the membrane
C with driling d.o.f. developed by Felippa et al.
C
C MODIFIED: 10-07-97 	K. H. Pierson
C Added stress calculations for: sigmaxx, sigmayy, sigmaxy,
C epsilonxx, epsilonyy, epsilonzz, epsilonxy and an equivalent
C strain similar to the vonmises stress. These stresses and strains
C can be calculated at the top, median or bottom surfaces.
C Stresses and strains are computed locally and then transformed
C to the global coordinate system.
C 
*-------------------------------------------------------------------*
*  CALLED BY : ThreeNodeShell.C 
*
*  SUBROUTINES CALLED:
*                      ROTATION
*                      TRIROT
*                      MEMBRA
*                      MOMEN
*			  TRANSFORMATION
*                      VONMISES
*-------------------------------------------------------------------*/

// t = thickness of triangle

       Eigen::Map<Eigen::Matrix<doublereal,3,1> > xl(_xl), yl(_yl), zl(_zl), h(_h);
       Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);
       Eigen::Map<Eigen::Matrix<doublereal,7,3> > stress(_stress);
   
       doublereal  epsxx,epsyy,epszz,epsxy;
       doublereal  rnxt,rnyt,rnxyt,rmxt,rmyt,rmxyt;
       doublereal  t2,sx,sy,sxy;
       Eigen::Matrix<doublereal,3,3> db,dm,r1;
       Eigen::Matrix<doublereal,18,1> dll;
       Eigen::Matrix<doublereal,3,1> xp,yp,zp,xlp,ylp,zlp,xg,yg,zg;
       Eigen::Matrix<doublereal,6,1> str;
       Eigen::Matrix<doublereal,6,6> tMat;
       Eigen::Matrix<int,9,1> lb,le;
       doublereal  cb,x21,y21,z21,x32,y32,z32,x13,y13,z13;
       doublereal  rx,ry,rz,bx,by,bz,rlr,rlb,bpr,area,ylr,zlr;
       doublereal  ycg,xcg,zcg,xlcg,ylcg,zlcg,f,t;
       Eigen::Matrix<doublereal,18,3> rmom,rmem;
       doublereal  rmx,rmy,rmxy,rnx,rny,rnxy,clr,cqr;
       doublereal  rmmx,rmmy,rmmxy,rnnx,rnny,rnnxy,sbf;
       doublereal  ebar;
       doublereal  factor;
       
       int i,j;
       lb << 2,3,4,8,9,10,14,15,16;
       le << 0,1,6,7,12,13,5,11,17;

       rmom.setZero();
       rmem.setZero();

        f   = 1.0;
        clr = 0.0;
        cqr = 1.0;
        t = h[0];

// set the bending constitutive matrix
        cb=e*(t*t*t)/12.0/(1.0-nu*nu);
        db(0,0) =     cb;
        db(0,1) =  nu*cb;
        db(0,2) = 0.0;
        db(1,0) =  db(0,1);
        db(1,1) =     cb;
        db(1,2) = 0.0;
        db(2,0) = 0.0;
        db(2,1) = 0.0;
        db(2,2) = 0.5*(1.0-nu)*cb;

// set the membrane constitutive matrix
        cb=e*t/(1.0-nu*nu);
        dm(0,0) =     cb;
        dm(0,1) =  nu*cb;
        dm(0,2) = 0.0;
        dm(1,0) =  dm(0,1);
        dm(1,1) =     cb;
        dm(1,2) = 0.0;
        dm(2,0) = 0.0;
        dm(2,1) = 0.0;
        dm(2,2) = 0.5*(1.0-nu)*cb;

// triangular dimension variables
        x21 = xl[1] - xl[0];
        y21 = yl[1] - yl[0];
        z21 = zl[1] - zl[0];
        x32 = xl[2] - xl[1];
        y32 = yl[2] - yl[1];
        z32 = zl[2] - zl[1];
        x13 = xl[0] - xl[2];
        y13 = yl[0] - yl[2];
        z13 = zl[0] - zl[2];
//  triangle in space : we compute the length of one side and the distance of the
//  opposing node to that side to compute the area
       rx   = x21;
       ry   = y21;
       rz   = z21;
       bx   = x32;
       by   = y32;
       bz   = z32;
       rlr  = sqrt( rx*rx + ry*ry + rz*rz );
       rlb  = sqrt( bx*bx + by*by + bz*bz );
       bpr  = abs(rx * bx + ry * by + rz *bz)/rlr;
       area = 0.5*rlr*sqrt(rlb*rlb-bpr*bpr);

// Direction cosines of the local system . X' is directed parallel to the 
// 2-1 side Z' is the external normal (counterclockwise). 
// Y' computed as Z' x X'
       xp[0] = x21/rlr;
       xp[1] = y21/rlr;
       xp[2] = z21/rlr;

// Z' local axis
       zp[0] = y21 * z32 - z21 * y32;
       zp[1] = z21 * x32 - x21 * z32;
       zp[2] = x21 * y32 - y21 * x32;
       zp.normalize();

// Y' local axis
       yp[0] = zp[1] * xp[2] - zp[2] * xp[1];
       yp[1] = zp[2] * xp[0] - zp[0] * xp[2];
       yp[2] = zp[0] * xp[1] - zp[1] * xp[0];
       yp.normalize();

// center of gravity
       xcg = (xl[0] + xl[1] + xl[2])/3.0;
       ycg = (yl[0] + yl[1] + yl[2])/3.0;
       zcg = (zl[0] + zl[1] + zl[2])/3.0;

// computing local coordinates
       for(i=0; i<3; ++i) {
         xlcg   = xl(i) - xcg;
         ylcg   = yl(i) - ycg;
         zlcg   = zl(i) - zcg;
         xlp[i] = xp[0] * xlcg + xp[1] * ylcg + xp[2] * zlcg;
         ylp[i] = yp[0] * xlcg + yp[1] * ylcg + yp[2] * zlcg;
         zlp[i] = zp[0] * xlcg + zp[1] * ylcg + zp[2] * zlcg;
       }

// Set Global axes

       xg[0] = 1.0;
       xg[1] = 0.0;
       xg[2] = 0.0;
       yg[0] = 0.0;
       yg[1] = 1.0;
       yg[2] = 0.0;
       zg[0] = 0.0;
       zg[1] = 0.0;
       zg[2] = 1.0;

// computing nodal rotation matrix
       rotation(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),r1.data());

// compute the von mises stress
// rotate nodal displacements to local system
       dll.setZero();
       for(i=0; i<3; ++i) {
         for(j=0; j<3; ++j) {
           dll[i]    = dll[i]    + r1(j,i)*v[j];
           dll[i+3]  = dll[i+3]  + r1(j,i)*v[j+3];
           dll[i+6]  = dll[i+6]  + r1(j,i)*v[j+6];
           dll[i+9]  = dll[i+9]  + r1(j,i)*v[j+9];
           dll[i+12] = dll[i+12] + r1(j,i)*v[j+12];
           dll[i+15] = dll[i+15] + r1(j,i)*v[j+15];
         }
       }

//      compute centroidal membrane strains
       membra(xlp.data(),ylp.data(),1.5,le.data(),rmem.data());

//      compute centroidal bending strains (curvatures (1/radius))
       momen(xlp.data(),ylp.data(),lb.data(),rmom.data());

//      pjsa 8/11/2011: the matrix rmom returned by momen function is
//      off by a factor of -1/sqrt(area) 
       factor = -1.0 / sqrt( area );

       rmx  = 0.0;
       rmy  = 0.0;
       rmxy = 0.0;
       rnx  = 0.0;
       rny  = 0.0;
       rnxy = 0.0;
       for(j=0; j<18; ++j) {
          rmx  =  rmx + factor*rmom(j,0)*dll[j];
          rmy  =  rmy + factor*rmom(j,1)*dll[j];
          rmxy = rmxy + factor*rmom(j,2)*dll[j];
          rnx  =  rnx + rmem(j,0)*dll[j];
          rny  =  rny + rmem(j,1)*dll[j];
          rnxy = rnxy + rmem(j,2)*dll[j];
       }

//      Subtract off thermal strain portions
       rnx = rnx - thrmStr;
       rny = rny - thrmStr;

//      COMPUTE VON MISES STRAIN RESULTANT

        if(strainFlg == true) {

          t2 = 0.5*t;

          rmx  = t2*rmx;
          rmy  = t2*rmy;
          rmxy = t2*rmxy;

// ... COMPUTE STRAINS AT MEDIAN SURFACE
          if(surface == 2) {
            epsxx = rnx;
            epsyy = rny;
            epsxy = 0.5*rnxy;
// ... COMPUTE STRAINS AT BOTTOM SURFACE
          } else if(surface == 3) {
            epsxx =  rnx -  rmx;
            epsyy =  rny -  rmy;
            epsxy = 0.5*(rnxy - rmxy);
// ... COMPUTE STRAINS AT TOP SURFACE
          } else {
            epsxx =  rnx +  rmx;
            epsyy =  rny +  rmy;
            epsxy = 0.5*(rnxy + rmxy);
          }

          epszz  = -nu/(1.0 - nu)*(epsxx + epsyy);

          str[0] = epsxx;
          str[1] = epsyy;
          str[2] = epszz;
          str[3] = epsxy;
          str[4] = 0.0;
          str[5] = 0.0;

          transform(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),str.data(),tMat.data());

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

          straineq(rmx,rmy,rmxy,rnx,rny,rnxy,nu,surface,ebar);

          stress(6,0) = ebar;
          stress(6,1) = ebar;
          stress(6,2) = ebar; 

          return;
        }

//     compute centroidal stress resultants

//     bending resultants
      rmmx  = db(0,0)*rmx + db(0,1)*rmy + db(0,2)*rmxy;
      rmmy  = db(1,0)*rmx + db(1,1)*rmy + db(1,2)*rmxy;
      rmmxy = db(2,0)*rmx + db(2,1)*rmy + db(2,2)*rmxy;

//     membrane resultants
      rnnx  = dm(0,0)*rnx + dm(0,1)*rny + dm(0,2)*rnxy;
      rnny  = dm(1,0)*rnx + dm(1,1)*rny + dm(1,2)*rnxy;
      rnnxy = dm(2,0)*rnx + dm(2,1)*rny + dm(2,2)*rnxy;

      t = abs(t);
      t2 = t*t;

      rnxt  =  rnnx/t;
      rnyt  =  rnny/t;
      rnxyt = rnnxy/t;

      rmxt  = 6.0* rmmx / t2;
      rmyt  = 6.0* rmmy / t2;
      rmxyt = 6.0*rmmxy / t2;

// COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT TOP SURFACE (DEFAULT)
      sx  =  rnxt +  rmxt;
      sy  =  rnyt +  rmyt;
      sxy = rnxyt + rmxyt;

// COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT MEDIAN SURFACE
      if(surface == 2) {
        sx  = rnxt;
        sy  = rnyt;
        sxy = rnxyt;
      }

// COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT BOTTOM SURFACE
      if(surface == 3) {
        sx  =  rnxt -  rmxt;
        sy  =  rnyt -  rmyt;
        sxy = rnxyt - rmxyt;
      }

      str[0] = sx;
      str[1] = sy;
      str[2] = 0.0; 
      str[3] = sxy;
      str[4] = 0.0;
      str[5] = 0.0;

      transform(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),str.data(),tMat.data());

      stress(0,1) = str[0];
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

//     compute von mises stress resultant

      vonmis(rmmx,rmmy,rmmxy,rnnx,rnny,rnnxy,t,sbf,surface);

      stress(6,0) = sbf;
      stress(6,1) = sbf;
      stress(6,2) = sbf;
}

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::vms8WRTdisp(doublereal *_xl,doublereal *_yl, doublereal *_zl,
              doublereal e, doublereal nu, doublereal *_h,
              doublereal *_v, doublereal *_vmsWRTdisp,
              bool strainFlg, int surface, doublereal thrmStr)
{
/*-----------------------------------------------------------------
C This subroutine evaluates the sensitivity of von mises stress 
C with respect to displacement at the centroid
C of the triangle. it computes the values at the top and bottom
C and picks the maximum as output
C for the spacial 18 d.o.f tree node triangle obtained as
C a combination of the aqr bending triangle plus the membrane
C with driling d.o.f. developed by Felippa et al.
C
C MODIFIED: 10-07-97 	K. H. Pierson
C Added stress calculations for: sigmaxx, sigmayy, sigmaxy,
C epsilonxx, epsilonyy, epsilonzz, epsilonxy and an equivalent
C strain similar to the vonmises stress. These stresses and strains
C can be calculated at the top, median or bottom surfaces.
C Stresses and strains are computed locally and then transformed
C to the global coordinate system.
C 
*-------------------------------------------------------------------*
*  CALLED BY : ThreeNodeShell.C 
*
*  SUBROUTINES CALLED:
*                      ROTATION
*                      TRIROT
*                      MEMBRA
*                      MOMEN
*			  TRANSFORMATION
*                      VONMISES
*-------------------------------------------------------------------*/

// t = thickness of triangle

       if(thrmStr != 0.0) {
         cerr << " ... Error: thermal stress should not be passed in sensitivity computation\n";
         exit(-1);
       }

       Eigen::Map<Eigen::Matrix<doublereal,3,1> > xl(_xl), yl(_yl), zl(_zl), h(_h);
       Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);
       Eigen::Map<Eigen::Matrix<doublereal,3,18> > vmsWRTdisp(_vmsWRTdisp);
   
       doublereal  epsxx,epsyy,epszz,epsxy;
       doublereal  rnxt,rnyt,rnxyt,rmxt,rmyt,rmxyt;
       doublereal  t2,sx,sy,sxy;
       Eigen::Matrix<doublereal,3,3> db,dm,r1;
       Eigen::Matrix<doublereal,18,18> R;
       Eigen::Matrix<doublereal,18,1> dll;
       Eigen::Matrix<doublereal,3,1> xp,yp,zp,xlp,ylp,zlp,xg,yg,zg;
       Eigen::Matrix<doublereal,6,1> str;
       Eigen::Matrix<int,9,1> lb,le;
       doublereal  cb,x21,y21,z21,x32,y32,z32,x13,y13,z13;
       doublereal  rx,ry,rz,bx,by,bz,rlr,rlb,bpr,area,ylr,zlr;
       doublereal  ycg,xcg,zcg,xlcg,ylcg,zlcg,f,t;
       Eigen::Matrix<doublereal,18,3> rmom,rmem;
       Eigen::Matrix<doublereal,6,18> dispToStrain;
       Eigen::Matrix<doublereal,6,6> tMat, Ct, dc;
       doublereal  rmx,rmy,rmxy,rnx,rny,rnxy,clr,cqr;
       doublereal  ebar, sbf;
       doublereal  factor, surfaceconst;
       Eigen::Matrix<doublereal, 6,1> stress_resultants, strain;
//       strain << rnx,rny,rnxy,rmx,rmy,rmxy;
       
       int i,j;
       lb << 2,3,4,8,9,10,14,15,16;
       le << 0,1,6,7,12,13,5,11,17;
       surfaceconst = 1;

       rmom.setZero();       rmem.setZero();
       tMat.setZero();       dispToStrain.setZero();
       Ct.setZero();         dc.setZero();
       db.setZero();         dm.setZero();
       R.setZero();          r1.setZero();
       dll.setZero();        strain.setZero();
       stress_resultants.setZero();

        f   = 1.0;
        clr = 0.0;
        cqr = 1.0;
        t = h[0];

// set the bending constitutive matrix
        cb=e*(t*t*t)/12.0/(1.0-nu*nu);
        db(0,0) =     cb;
        db(0,1) =  nu*cb;
        db(1,0) =  db(0,1);
        db(1,1) =     cb;
        db(2,2) = 0.5*(1.0-nu)*cb;

// set the membrane constitutive matrix
        cb=e*t/(1.0-nu*nu);
        dm(0,0) =     cb;
        dm(0,1) =  nu*cb;
        dm(1,0) =  dm(0,1);
        dm(1,1) =     cb;
        dm(2,2) = 0.5*(1.0-nu)*cb;

// set the global constitutive matrix
        dc << dm, Eigen::Matrix<doublereal,3,3>::Zero(),
              Eigen::Matrix<doublereal,3,3>::Zero(), db;

// triangular dimension variables
        x21 = xl[1] - xl[0];
        y21 = yl[1] - yl[0];
        z21 = zl[1] - zl[0];
        x32 = xl[2] - xl[1];
        y32 = yl[2] - yl[1];
        z32 = zl[2] - zl[1];
        x13 = xl[0] - xl[2];
        y13 = yl[0] - yl[2];
        z13 = zl[0] - zl[2];
//  triangle in space : we compute the length of one side and the distance of the
//  opposing node to that side to compute the area
       rx   = x21;
       ry   = y21;
       rz   = z21;
       bx   = x32;
       by   = y32;
       bz   = z32;
       rlr  = sqrt( rx*rx + ry*ry + rz*rz );
       rlb  = sqrt( bx*bx + by*by + bz*bz );
       bpr  = abs(rx * bx + ry * by + rz *bz)/rlr;
       area = 0.5*rlr*sqrt(rlb*rlb-bpr*bpr);

// Direction cosines of the local system . X' is directed parallel to the 
// 2-1 side Z' is the external normal (counterclockwise). 
// Y' computed as Z' x X'
       xp[0] = x21/rlr;
       xp[1] = y21/rlr;
       xp[2] = z21/rlr;

// Z' local axis
       zp[0] = y21 * z32 - z21 * y32;
       zp[1] = z21 * x32 - x21 * z32;
       zp[2] = x21 * y32 - y21 * x32;
       zp.normalize();

// Y' local axis
       yp[0] = zp[1] * xp[2] - zp[2] * xp[1];
       yp[1] = zp[2] * xp[0] - zp[0] * xp[2];
       yp[2] = zp[0] * xp[1] - zp[1] * xp[0];
       yp.normalize();

// center of gravity
       xcg = (xl[0] + xl[1] + xl[2])/3.0;
       ycg = (yl[0] + yl[1] + yl[2])/3.0;
       zcg = (zl[0] + zl[1] + zl[2])/3.0;

// computing local coordinates
       for(i=0; i<3; ++i) {
         xlcg   = xl(i) - xcg;
         ylcg   = yl(i) - ycg;
         zlcg   = zl(i) - zcg;
         xlp[i] = xp[0] * xlcg + xp[1] * ylcg + xp[2] * zlcg;
         ylp[i] = yp[0] * xlcg + yp[1] * ylcg + yp[2] * zlcg;
         zlp[i] = zp[0] * xlcg + zp[1] * ylcg + zp[2] * zlcg;
       }

// Set Global axes

       xg[0] = 1.0;
       xg[1] = 0.0;
       xg[2] = 0.0;
       yg[0] = 0.0;
       yg[1] = 1.0;
       yg[2] = 0.0;
       zg[0] = 0.0;
       zg[1] = 0.0;
       zg[2] = 1.0;

// computing nodal rotation matrix
       rotation(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),r1.data());
       Eigen::Matrix<doublereal,3,3> zero33;
       zero33.setZero();
       R << r1.transpose()    , zero33, zero33, zero33, zero33, zero33,
            zero33, r1.transpose()    , zero33, zero33, zero33, zero33,
            zero33, zero33, r1.transpose()    , zero33, zero33, zero33,
            zero33, zero33, zero33, r1.transpose()    , zero33, zero33,
            zero33, zero33, zero33, zero33, r1.transpose()    , zero33,
            zero33, zero33, zero33, zero33, zero33, r1.transpose();

// compute the von mises stress
// rotate nodal displacements to local system
       for(i=0; i<3; ++i) {
         for(j=0; j<3; ++j) {
           dll[i]    = dll[i]    + r1(j,i)*v[j];
           dll[i+3]  = dll[i+3]  + r1(j,i)*v[j+3];
           dll[i+6]  = dll[i+6]  + r1(j,i)*v[j+6];
           dll[i+9]  = dll[i+9]  + r1(j,i)*v[j+9];
           dll[i+12] = dll[i+12] + r1(j,i)*v[j+12];
           dll[i+15] = dll[i+15] + r1(j,i)*v[j+15];
         }
       }

//      compute centroidal membrane strains
       membra(xlp.data(),ylp.data(),1.5,le.data(),rmem.data());

//      compute centroidal bending strains (curvatures (1/radius))
       momen(xlp.data(),ylp.data(),lb.data(),rmom.data());

//      pjsa 8/11/2011: the matrix rmom returned by momen function is
//      off by a factor of -1/sqrt(area) 
       factor = -1.0 / sqrt( area );

       dispToStrain << rmem.transpose(), factor*rmom.transpose();

       strain = dispToStrain * dll;

//      Subtract off thermal strain portions
//       strain[0] = strain[0] - thrmStr;
//       strain[1] = strain[1] - thrmStr;

//      COMPUTE VON MISES STRAIN RESULTANT
/*
        if(strainFlg == true) {

          t2 = 0.5*t;

          rmx  = t2*rmx;
          rmy  = t2*rmy;
          rmxy = t2*rmxy;

// ... COMPUTE STRAINS AT MEDIAN SURFACE
          if(surface == 2) {
            epsxx = rnx;
            epsyy = rny;
            epsxy = 0.5*rnxy;
// ... COMPUTE STRAINS AT BOTTOM SURFACE
          } else if(surface == 3) {
            epsxx =  rnx -  rmx;
            epsyy =  rny -  rmy;
            epsxy = 0.5*(rnxy - rmxy);
// ... COMPUTE STRAINS AT TOP SURFACE
          } else {
            epsxx =  rnx +  rmx;
            epsyy =  rny +  rmy;
            epsxy = 0.5*(rnxy + rmxy);
          }

          epszz  = -nu/(1.0 - nu)*(epsxx + epsyy);

          str[0] = epsxx;
          str[1] = epsyy;
          str[2] = epszz;
          str[3] = epsxy;
          str[4] = 0.0;
          str[5] = 0.0;

          transform(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),str.data(),tMat.data());
          straineq(rmx,rmy,rmxy,rnx,rny,rnxy,nu,surface,ebar);

          cerr << "print von Mises strain:\n" << ebar << endl;

          return;
        }
*/

//     compute centroidal stress resultants
      stress_resultants = dc * strain;

      t = abs(t);
      t2 = t*t;

      Ct(0,0) = Ct(1,1) = Ct(2,2) = 1.0/t;
      Ct(3,3) = Ct(4,4) = Ct(5,5) = 6.0/t2; 

      Eigen::Matrix<doublereal,6,1> Sr;
      Sr = Ct*stress_resultants;

// COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT TOP SURFACE (DEFAULT)
      sx  =  Sr[0] +  Sr[3];
      sy  =  Sr[1] +  Sr[4];
      sxy =  Sr[2] +  Sr[5];

// COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT MEDIAN SURFACE
      if(surface == 2) {
        sx  = Sr[0];
        sy  = Sr[1];
        sxy = Sr[2];
        surfaceconst = 0;
      }

// COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT BOTTOM SURFACE
      if(surface == 3) {
        sx  =  Sr[0] -  Sr[3];
        sy  =  Sr[1] -  Sr[4];
        sxy =  Sr[2] -  Sr[5];
        surfaceconst = -1;
      }

      str[0] = sx;
      str[1] = sy;
      str[2] = 0.0; 
      str[3] = sxy;
      str[4] = 0.0;
      str[5] = 0.0;

      Eigen::Matrix<doublereal,6,6> dstressdsr;
      dstressdsr.setZero();
      dstressdsr(0,0) = dstressdsr(1,1) = dstressdsr(3,2) = 1;
      dstressdsr(0,3) = dstressdsr(1,4) = dstressdsr(3,5) = surfaceconst;

      transform(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),str.data(),tMat.data());

//     compute von mises stress resultant

      vonmis(stress_resultants[3],stress_resultants[4],stress_resultants[5],
             stress_resultants[0],stress_resultants[1],stress_resultants[2],t,sbf,surface);
      
      doublereal stress_0dxy = str[0] - str[1];
      doublereal stress_0dyz = str[1] - str[2];
      doublereal stress_0dxz = str[0] - str[2];

      Eigen::Matrix<doublereal,3,6> dvmsdStress;
      dvmsdStress.setZero();
      dvmsdStress(2,0) = dvmsdStress(1,0) = dvmsdStress(0,0) = (2.*str[0]-str[1]-str[2])/(2.*sbf);
      dvmsdStress(2,1) = dvmsdStress(1,1) = dvmsdStress(0,1) = (2.*str[1]-str[0]-str[2])/(2.*sbf);
      dvmsdStress(2,2) = dvmsdStress(1,2) = dvmsdStress(0,2) = (2.*str[2]-str[1]-str[0])/(2.*sbf);
      dvmsdStress(2,3) = dvmsdStress(1,3) = dvmsdStress(0,3) = (3.*str[3])/sbf;
      dvmsdStress(2,4) = dvmsdStress(1,4) = dvmsdStress(0,4) = (3.*str[4])/sbf;
      dvmsdStress(2,5) = dvmsdStress(1,5) = dvmsdStress(0,5) = (3.*str[5])/sbf;

      vmsWRTdisp = dvmsdStress * (tMat * (dstressdsr * (Ct * (dc * (dispToStrain * R)))));

}

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::vms8WRTthic(doublereal *_xl,doublereal *_yl, doublereal *_zl,
              doublereal e, doublereal nu, doublereal *_h,
              doublereal *_v, doublereal *_vmsWRTthic,
              bool strainFlg, int surface, doublereal thrmStr)
{
// t = thickness of triangle

       if(thrmStr != 0.0) {
         cerr << " ... Error: thermal stress should not be passed in sensitivity computation\n";
         exit(-1);
       }

       Eigen::Map<Eigen::Matrix<doublereal,3,1> > xl(_xl), yl(_yl), zl(_zl), h(_h);
       Eigen::Map<Eigen::Matrix<doublereal,18,1> > v(_v);
       Eigen::Map<Eigen::Matrix<doublereal,3,1> > vmsWRTthic(_vmsWRTthic);
   
       doublereal  epsxx,epsyy,epszz,epsxy;
       doublereal  rnxt,rnyt,rnxyt,rmxt,rmyt,rmxyt;
       doublereal  t2,sx,sy,sxy;
       Eigen::Matrix<doublereal,3,3> db,dm,r1;
       Eigen::Matrix<doublereal,3,3> ddbdh,ddmdh;
       Eigen::Matrix<doublereal,18,18> R;
       Eigen::Matrix<doublereal,18,1> dll;
       Eigen::Matrix<doublereal,3,1> xp,yp,zp,xlp,ylp,zlp,xg,yg,zg;
       Eigen::Matrix<doublereal,6,1> str;
       Eigen::Matrix<int,9,1> lb,le;
       doublereal  cb,x21,y21,z21,x32,y32,z32,x13,y13,z13;
       doublereal  rx,ry,rz,bx,by,bz,rlr,rlb,bpr,area,ylr,zlr;
       doublereal  ycg,xcg,zcg,xlcg,ylcg,zlcg,f,t;
       Eigen::Matrix<doublereal,18,3> rmom,rmem;
       Eigen::Matrix<doublereal,6,18> dispToStrain;
       Eigen::Matrix<doublereal,6,6> tMat, Ct, dc, ddcdh, dCtdh;
       doublereal  rmx,rmy,rmxy,rnx,rny,rnxy,clr,cqr;
       doublereal  ebar, sbf;
       doublereal  factor, surfaceconst;
       Eigen::Matrix<doublereal, 6,1> stress_resultants, strain;
//       strain << rnx,rny,rnxy,rmx,rmy,rmxy;
       
       int i,j;
       lb << 2,3,4,8,9,10,14,15,16;
       le << 0,1,6,7,12,13,5,11,17;
       surfaceconst = 1;

       rmom.setZero();       rmem.setZero();
       tMat.setZero();       dispToStrain.setZero();
       Ct.setZero();         dc.setZero();
       dCtdh.setZero();      ddcdh.setZero();
       db.setZero();         dm.setZero();
       ddbdh.setZero();      ddmdh.setZero();
       R.setZero();          r1.setZero();
       dll.setZero();        strain.setZero();
       stress_resultants.setZero();

        f   = 1.0;
        clr = 0.0;
        cqr = 1.0;
        t = h[0];

// set the bending constitutive matrix
        cb=e*(t*t*t)/12.0/(1.0-nu*nu);
        doublereal dcbdh=e*(t*t)/4.0/(1.0-nu*nu);
        db(0,0) =     cb;
        db(0,1) =  nu*cb;
        db(1,0) =  db(0,1);
        db(1,1) =     cb;
        db(2,2) = 0.5*(1.0-nu)*cb;
        ddbdh(0,0) =     dcbdh;
        ddbdh(0,1) =  nu*dcbdh;
        ddbdh(1,0) =  ddbdh(0,1);
        ddbdh(1,1) =     dcbdh;
        ddbdh(2,2) = 0.5*(1.0-nu)*dcbdh;

// set the membrane constitutive matrix
        cb=e*t/(1.0-nu*nu);
        dcbdh=e/(1.0-nu*nu);
        dm(0,0) =     cb;
        dm(0,1) =  nu*cb;
        dm(1,0) =  dm(0,1);
        dm(1,1) =     cb;
        dm(2,2) = 0.5*(1.0-nu)*cb;
        ddmdh(0,0) =     dcbdh;
        ddmdh(0,1) =  nu*dcbdh;
        ddmdh(1,0) =  ddmdh(0,1);
        ddmdh(1,1) =     dcbdh;
        ddmdh(2,2) = 0.5*(1.0-nu)*dcbdh;

// set the global constitutive matrix
        dc << dm, Eigen::Matrix<doublereal,3,3>::Zero(),
              Eigen::Matrix<doublereal,3,3>::Zero(), db;
        ddcdh << ddmdh, Eigen::Matrix<doublereal,3,3>::Zero(),
                 Eigen::Matrix<doublereal,3,3>::Zero(), ddbdh;

// triangular dimension variables
        x21 = xl[1] - xl[0];
        y21 = yl[1] - yl[0];
        z21 = zl[1] - zl[0];
        x32 = xl[2] - xl[1];
        y32 = yl[2] - yl[1];
        z32 = zl[2] - zl[1];
        x13 = xl[0] - xl[2];
        y13 = yl[0] - yl[2];
        z13 = zl[0] - zl[2];
//  triangle in space : we compute the length of one side and the distance of the
//  opposing node to that side to compute the area
       rx   = x21;
       ry   = y21;
       rz   = z21;
       bx   = x32;
       by   = y32;
       bz   = z32;
       rlr  = sqrt( rx*rx + ry*ry + rz*rz );
       rlb  = sqrt( bx*bx + by*by + bz*bz );
       bpr  = abs(rx * bx + ry * by + rz *bz)/rlr;
       area = 0.5*rlr*sqrt(rlb*rlb-bpr*bpr);

// Direction cosines of the local system . X' is directed parallel to the 
// 2-1 side Z' is the external normal (counterclockwise). 
// Y' computed as Z' x X'
       xp[0] = x21/rlr;
       xp[1] = y21/rlr;
       xp[2] = z21/rlr;

// Z' local axis
       zp[0] = y21 * z32 - z21 * y32;
       zp[1] = z21 * x32 - x21 * z32;
       zp[2] = x21 * y32 - y21 * x32;
       zp.normalize();

// Y' local axis
       yp[0] = zp[1] * xp[2] - zp[2] * xp[1];
       yp[1] = zp[2] * xp[0] - zp[0] * xp[2];
       yp[2] = zp[0] * xp[1] - zp[1] * xp[0];
       yp.normalize();

// center of gravity
       xcg = (xl[0] + xl[1] + xl[2])/3.0;
       ycg = (yl[0] + yl[1] + yl[2])/3.0;
       zcg = (zl[0] + zl[1] + zl[2])/3.0;

// computing local coordinates
       for(i=0; i<3; ++i) {
         xlcg   = xl(i) - xcg;
         ylcg   = yl(i) - ycg;
         zlcg   = zl(i) - zcg;
         xlp[i] = xp[0] * xlcg + xp[1] * ylcg + xp[2] * zlcg;
         ylp[i] = yp[0] * xlcg + yp[1] * ylcg + yp[2] * zlcg;
         zlp[i] = zp[0] * xlcg + zp[1] * ylcg + zp[2] * zlcg;
       }

// Set Global axes

       xg[0] = 1.0;
       xg[1] = 0.0;
       xg[2] = 0.0;
       yg[0] = 0.0;
       yg[1] = 1.0;
       yg[2] = 0.0;
       zg[0] = 0.0;
       zg[1] = 0.0;
       zg[2] = 1.0;

// computing nodal rotation matrix
       rotation(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),r1.data());
       Eigen::Matrix<doublereal,3,3> zero33;
       zero33.setZero();
       R << r1.transpose()    , zero33, zero33, zero33, zero33, zero33,
            zero33, r1.transpose()    , zero33, zero33, zero33, zero33,
            zero33, zero33, r1.transpose()    , zero33, zero33, zero33,
            zero33, zero33, zero33, r1.transpose()    , zero33, zero33,
            zero33, zero33, zero33, zero33, r1.transpose()    , zero33,
            zero33, zero33, zero33, zero33, zero33, r1.transpose();

// compute the von mises stress
// rotate nodal displacements to local system
       for(i=0; i<3; ++i) {
         for(j=0; j<3; ++j) {
           dll[i]    = dll[i]    + r1(j,i)*v[j];
           dll[i+3]  = dll[i+3]  + r1(j,i)*v[j+3];
           dll[i+6]  = dll[i+6]  + r1(j,i)*v[j+6];
           dll[i+9]  = dll[i+9]  + r1(j,i)*v[j+9];
           dll[i+12] = dll[i+12] + r1(j,i)*v[j+12];
           dll[i+15] = dll[i+15] + r1(j,i)*v[j+15];
         }
       }

//      compute centroidal membrane strains
       membra(xlp.data(),ylp.data(),1.5,le.data(),rmem.data());

//      compute centroidal bending strains (curvatures (1/radius))
       momen(xlp.data(),ylp.data(),lb.data(),rmom.data());

//      pjsa 8/11/2011: the matrix rmom returned by momen function is
//      off by a factor of -1/sqrt(area) 
       factor = -1.0 / sqrt( area );

       dispToStrain << rmem.transpose(), factor*rmom.transpose();

       strain = dispToStrain * dll;

//      Subtract off thermal strain portions
//       strain[0] = strain[0] - thrmStr;
//       strain[1] = strain[1] - thrmStr;

//      COMPUTE VON MISES STRAIN RESULTANT
/*
        if(strainFlg == true) {

          t2 = 0.5*t;

          rmx  = t2*rmx;
          rmy  = t2*rmy;
          rmxy = t2*rmxy;

// ... COMPUTE STRAINS AT MEDIAN SURFACE
          if(surface == 2) {
            epsxx = rnx;
            epsyy = rny;
            epsxy = 0.5*rnxy;
// ... COMPUTE STRAINS AT BOTTOM SURFACE
          } else if(surface == 3) {
            epsxx =  rnx -  rmx;
            epsyy =  rny -  rmy;
            epsxy = 0.5*(rnxy - rmxy);
// ... COMPUTE STRAINS AT TOP SURFACE
          } else {
            epsxx =  rnx +  rmx;
            epsyy =  rny +  rmy;
            epsxy = 0.5*(rnxy + rmxy);
          }

          epszz  = -nu/(1.0 - nu)*(epsxx + epsyy);

          str[0] = epsxx;
          str[1] = epsyy;
          str[2] = epszz;
          str[3] = epsxy;
          str[4] = 0.0;
          str[5] = 0.0;

          transform(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),str.data(),tMat.data());
          straineq(rmx,rmy,rmxy,rnx,rny,rnxy,nu,surface,ebar);

          cerr << "print von Mises strain:\n" << ebar << endl;

          return;
        }
*/

//     compute centroidal stress resultants
      stress_resultants = dc * strain;

      t = abs(t);
      t2 = t*t;
      doublereal t3 = t*t*t;

      Ct(0,0) = Ct(1,1) = Ct(2,2) = 1.0/t;
      Ct(3,3) = Ct(4,4) = Ct(5,5) = 6.0/t2; 
      dCtdh(0,0) = dCtdh(1,1) = dCtdh(2,2) = -1.0/t2;
      dCtdh(3,3) = dCtdh(4,4) = dCtdh(5,5) = -12.0/t3; 

      Eigen::Matrix<doublereal,6,1> Sr;
      Sr = Ct*stress_resultants;

// COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT TOP SURFACE (DEFAULT)
      sx  =  Sr[0] +  Sr[3];
      sy  =  Sr[1] +  Sr[4];
      sxy =  Sr[2] +  Sr[5];

// COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT MEDIAN SURFACE
      if(surface == 2) {
        sx  = Sr[0];
        sy  = Sr[1];
        sxy = Sr[2];
        surfaceconst = 0;
      }

// COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT BOTTOM SURFACE
      if(surface == 3) {
        sx  =  Sr[0] -  Sr[3];
        sy  =  Sr[1] -  Sr[4];
        sxy =  Sr[2] -  Sr[5];
        surfaceconst = -1;
      }

      str[0] = sx;
      str[1] = sy;
      str[2] = 0.0; 
      str[3] = sxy;
      str[4] = 0.0;
      str[5] = 0.0;

      Eigen::Matrix<doublereal,6,6> dstressdsr;
      dstressdsr.setZero();
      dstressdsr(0,0) = dstressdsr(1,1) = dstressdsr(3,2) = 1;
      dstressdsr(0,3) = dstressdsr(1,4) = dstressdsr(3,5) = surfaceconst;

      transform(xp.data(),yp.data(),zp.data(),xg.data(),yg.data(),zg.data(),str.data(),tMat.data());

//     compute von mises stress resultant

      vonmis(stress_resultants[3],stress_resultants[4],stress_resultants[5],
             stress_resultants[0],stress_resultants[1],stress_resultants[2],t,sbf,surface);
      
      doublereal stress_0dxy = str[0] - str[1];
      doublereal stress_0dyz = str[1] - str[2];
      doublereal stress_0dxz = str[0] - str[2];

      Eigen::Matrix<doublereal,3,6> dvmsdStress;
      dvmsdStress.setZero();
      dvmsdStress(2,0) = dvmsdStress(1,0) = dvmsdStress(0,0) = (2.*str[0]-str[1]-str[2])/(2.*sbf);
      dvmsdStress(2,1) = dvmsdStress(1,1) = dvmsdStress(0,1) = (2.*str[1]-str[0]-str[2])/(2.*sbf);
      dvmsdStress(2,2) = dvmsdStress(1,2) = dvmsdStress(0,2) = (2.*str[2]-str[1]-str[0])/(2.*sbf);
      dvmsdStress(2,3) = dvmsdStress(1,3) = dvmsdStress(0,3) = (3.*str[3])/sbf;
      dvmsdStress(2,4) = dvmsdStress(1,4) = dvmsdStress(0,4) = (3.*str[4])/sbf;
      dvmsdStress(2,5) = dvmsdStress(1,5) = dvmsdStress(0,5) = (3.*str[5])/sbf;

      vmsWRTthic = dvmsdStress * (tMat * (dstressdsr * ((dCtdh * dc + Ct * ddcdh) * (dispToStrain * (R * v)))));

}
template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
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
ShellElementSemiTemplate<doublereal>
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
ShellElementSemiTemplate<doublereal>
::transform(doublereal *_xl, doublereal *_yl, doublereal *_zl,
            doublereal *_xg, doublereal *_yg, doublereal *_zg,
            doublereal *_str, doublereal *_tMat)
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
        Eigen::Map<Eigen::Matrix<doublereal,6,6> > tMat(_tMat);

// Local Declarations

        doublereal l1,l2,l3;
        doublereal m1,m2,m3;
        doublereal n1,n2,n3;
        Eigen::Matrix<doublereal,6,1> s;

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
       
        tMat(0,0) = l1*l1;
        tMat(0,1) = m1*m1;
        tMat(0,2) = n1*n1;
        tMat(0,3) = 2.0*l1*m1;
        tMat(0,4) = 2.0*m1*n1;
        tMat(0,5) = 2.0*n1*l1;
       
        tMat(1,0) = l2*l2;
        tMat(1,1) = m2*m2;
        tMat(1,2) = n2*n2;
        tMat(1,3) = 2.0*l2*m2;
        tMat(1,4) = 2.0*m2*n2;
        tMat(1,5) = 2.0*n2*l2;
       
        tMat(2,0) = l3*l3;
        tMat(2,1) = m3*m3;
        tMat(2,2) = n3*n3;
        tMat(2,3) = 2.0*l3*m3;
        tMat(2,4) = 2.0*m3*n3;
        tMat(2,5) = 2.0*n3*l3;
       
        tMat(3,0) = l1*l2;
        tMat(3,1) = m1*m2;
        tMat(3,2) = n1*n2;
        tMat(3,3) = l1*m2 + l2*m1;
        tMat(3,4) = m1*n2 + m2*n1;
        tMat(3,5) = n1*l2 + n2*l1;
       
        tMat(4,0) = l2*l3;
        tMat(4,1) = m2*m3;
        tMat(4,2) = n2*n3;
        tMat(4,3) = l2*m3 + l3*m2;
        tMat(4,4) = m2*n3 + m3*n2;
        tMat(4,5) = n2*l3 + n3*l2;
       
        tMat(5,0) = l3*l1;
        tMat(5,1) = m3*m1;
        tMat(5,2) = n3*n1;
        tMat(5,3) = l3*m1 + l1*m3;
        tMat(5,4) = m3*n1 + m1*n3;
        tMat(5,5) = n3*l1 + n1*l3;

// Perform the multiplication {str'} = T{str}
       
        str = tMat*s;
} 

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::vonmis(doublereal rmx, doublereal rmy, doublereal rmxy,
         doublereal rnx, doublereal rny, doublereal rnxy,
         doublereal &t,  doublereal &sv, int surface)
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

      if(surface == 3) {
        sv = sb;
        return;
      }

// ... COMPUTE VON MISES STRESS IN MEDIAN SURFACE
      if(surface == 2) {
        sx  = rnxt;
        sy  = rnyt;
        sxy = rnxyt;
        compj2(sx,sy,sxy,sm);
        sv = sq3 * sm;
        return;
      }

// ... COMPUTE VON MISES STRESS IN TOP SURFACE
      sx  =  rnxt +  rmxt;
      sy  =  rnyt +  rmyt;
      sxy = rnxyt + rmxyt;
      compj2(sx,sy,sxy,st);
      st = sq3 *st;

      if(surface == 1) {
        sv = st;
        return;
      }
}

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
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

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::momen(doublereal *_x, doublereal *_y, int *_lb, doublereal *_l)
{
      Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y);
      Eigen::Map<Eigen::Matrix<int,9,1> > lb(_lb);
      Eigen::Map<Eigen::Matrix<doublereal,18,3> > l(_l);

      Eigen::Matrix<doublereal,9,3> lqr;
      doublereal  x0, y0, cab, a1, a2, a3, b1, b2, b3;
      doublereal  x21, x32, x13, y21, y32, y13;
      doublereal  x12, x23, x31, y12, y23, y31;
      doublereal  xl12, xl23, xl31, c12, c23, c31, s12, s23, s31;
      doublereal  cc12, cc23, cc31, ss12, ss23, ss31;
      doublereal  cs12, cs23, cs31;
      doublereal  area2;
      int         i, j,  kk;

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
        cerr << " ... Error! basico: negative area\n";
        if (area2 == 0.0)   cerr << " ... Error! basico: zero area\n";
        exit(-1);
      }
      x0 =      (x[0]+x[1]+x[2])/3.;
      y0 =      (y[0]+y[1]+y[2])/3.;
      cab =      3.0/area2;
      a1 =      -cab*(y[2]-y0);
      a2 =      -cab*(y[0]-y0);
      a3 =      -cab*(y[1]-y0);
      b1 =       cab*(x[2]-x0);
      b2 =       cab*(x[0]-x0);
      b3 =       cab*(x[1]-x0);
      xl12 =    sqrt(x12*x12+y12*y12);
      xl23 =    sqrt(x23*x23+y23*y23);
      xl31 =    sqrt(x31*x31+y31*y31);
      lqr.setZero();

        c12 =       y21/xl12;
        s12 =       x12/xl12;
        c23 =       y32/xl23;
        s23 =       x23/xl23;
        c31 =       y13/xl31;
        s31 =       x31/xl31;
        cc12 =      c12*c12;
        cc23 =      c23*c23;
        cc31 =      c31*c31;
        ss12 =      s12*s12;
        ss23 =      s23*s23;
        ss31 =      s31*s31;
        cs12 =      c12*s12;
        cs23 =      c23*s23;
        cs31 =      c31*s31;
        lqr(0,0) =   cs12 - cs31;
        lqr(0,1) =  -lqr(0,0);
        lqr(0,2) =  (cc31-ss31) - (cc12-ss12);
        lqr(1,0) =  (cc12*x12 + cc31*x31)*.5;
        lqr(1,1) =  (ss12*x12 + ss31*x31)*.5;
        lqr(1,2) =   ss12*y21 + ss31*y13;
        lqr(2,0) = -(cc12*y21 + cc31*y13)*.5;
        lqr(2,1) = -0.5*lqr(1,2);
        lqr(2,2) =  -2.*lqr(1,0);
        lqr(3,0) =  cs23 - cs12;
        lqr(3,1) =  -lqr(3,0);
        lqr(3,2) =  (cc12-ss12) - (cc23-ss23);
        lqr(4,0) =  (cc12*x12 + cc23*x23)*.5;
        lqr(4,1) =  (ss12*x12 + ss23*x23)*.5;
        lqr(4,2) =   ss12*y21 + ss23*y32;
        lqr(5,0) = -(cc12*y21 + cc23*y32)*.5;
        lqr(5,1) = -0.5*lqr(4,2);
        lqr(5,2) =  -2.*lqr(4,0);
        lqr(6,0) =  cs31 - cs23;
        lqr(6,1) =  -lqr(6,0);
        lqr(6,2) =  (cc23-ss23) - (cc31-ss31);
        lqr(7,0) =  (cc23*x23 + cc31*x31)*.5;
        lqr(7,1) =  (ss23*x23 + ss31*x31)*.5;
        lqr(7,2) =   ss23*y32 + ss31*y13;
        lqr(8,0) = -(cc23*y32 + cc31*y13)*.5;
        lqr(8,1) = -0.5*lqr(7,2);
        lqr(8,2) =  -2.*lqr(7,0);

      for(j = 0; j<9; ++j) {
        kk= lb[j];
        l(kk,0) = lqr(j,0)*sqrt(2.0/area2);
        l(kk,1) = lqr(j,1)*sqrt(2.0/area2);
        l(kk,2) = lqr(j,2)*sqrt(2.0/area2);
      }
}

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::straineq(doublereal rmx, doublereal rmy, doublereal rmxy,
           doublereal rnx, doublereal rny, doublereal rnxy,
           doublereal nu,  int surface, doublereal &ebar)
{
/**************************************************************
*       THIS ROUTINE CALCULATES THE EQUIVALENT VON MISES      *
*       STRAIN FOR THE 3 NODE SHELL ELEMENT                   *
***************************************************************
*                                                             *
*       AUTHOR  :       K.H. PIERSON                          *
*       DATE    :       MARCH 1997                            *
*       VERSION :       FEM-C++ 1.00                          *
*                                                             *
***************************************************************
*                                                             *
*	rmx  = centroidal bending curvature (kxxc)            *
*	rmy  = centroidal bending curvature (kyyc)            *
*	rmxy = centroidal bending curvature (kxyc)            *
*	rnx  = centroidal membrane strain   (exxc)            *
*	rny  = centroidal membrane strain   (eyyc)            *
*	rnxy = centroidal membrane strain   (exyc)            *
*	nu   = Poisson's ratio                                *
*	ebar = equivalent strain                              *
*                                                             *
***************************************************************
*                                                             *
*       CALLED BY :  SANDS8.F                                 *
*                                                             *
***************************************************************/

        doublereal ex,ey,ez,exy,etop,ebot,emid;

// ... CALCULATE STRAINS AT MEDIAN SURFACE,
// ... COMPUTE EQUIVALENT STRAIN AT MEDIAN SURFACE,
// ... AND RETURN IF MEDIAN SURFACE IS REQUESTED
        if(surface == 2) {
          ex  = rnx;
          ey  = rny;
          ez  = -nu/(1.0-nu)*(ex + ey);
          exy = 0.5*rnxy;
          equiv(ex,ey,ez,exy,emid);
          ebar = emid;
          return;
        }

// ... CALCULATE STRAINS AT TOP SURFACE
        ex  = rnx + rmx;
        ey  = rny + rmy;
        ez  = -nu/(1.0-nu)*(ex + ey);
        exy = 0.5*(rnxy + rmxy);

// ... COMPUTE EQUIVALENT STRAIN AT TOP SURFACE
        equiv(ex,ey,ez,exy,etop);

// ... RETURN IF TOP SURFACE VALUE IS REQUESTED
        if(surface == 1) {
          ebar = etop;
          return;
        }

// ... CALCULATE STRAINS AT BOTTOM SURFACE
        ex   = rnx - rmx;
        ey   = rny - rmy;
        ez   = -nu/(1.0-nu)*(ex + ey);
        exy  = 0.5*(rnxy - rmxy);

// ... COMPUTE EQUIVALENT STRAIN AT BOTTOM SURFACE
        equiv(ex,ey,ez,exy,ebot);

// ... RETURN IF BOTTOM SURFACE VALUE IS REQUESTED
        if(surface == 3) {
          ebar = ebot;
          return;
        }
}


template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::equiv(doublereal ex, doublereal ey, doublereal ez, doublereal exy, doublereal &eq)
{
// ... LOCAL VARIABLES
        doublereal e0,dex,dey,dez;
        
// ... COMPUTE MEAN HYDROSTATIC STRAIN
        e0 = (ex + ey + ez)/3.0;

// ... COMPUTE DEVIATORIC STRAINS
        dex = ex - e0;
        dey = ey - e0;
        dez = ez - e0;

// ... COMPUTE EQUIVALENT STRAIN
        eq = ((dex*dex + dey*dey + dez*dez)/2.0) + (exy*exy);
        eq = sqrt(3.0 * eq);
}

#endif
#endif
