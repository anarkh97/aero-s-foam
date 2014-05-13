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
*                      TRANSFORMATION
*                      VONMISES
*-------------------------------------------------------------------*/

// t = thickness of triangle
       using std::abs;
       using std::sqrt;

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
       using std::abs;
       using std::sqrt;

       if(thrmStr != 0.0) {
         std::cerr << " ... Error: thermal stress should not be passed in sensitivity computation\n";
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

          std::cerr << "print von Mises strain:\n" << ebar << std::endl;

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
       using std::abs;
       using std::sqrt;

       if(thrmStr != 0.0) {
         std::cerr << " ... Error: thermal stress should not be passed in sensitivity computation\n";
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

          std::cerr << "print von Mises strain:\n" << ebar << std::endl;

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
        std::cerr << " ... Error! SM3MB: Zero area\n";
        if (area2 == 0.0)   std::cerr << " ... Error! SM3MB: Zero area\n";
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
      using std::abs;
      using std::sqrt;
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
      using std::sqrt;
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
      using std::sqrt;
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
        std::cerr << " ... Error! basico: negative area\n";
        if (area2 == 0.0)   std::cerr << " ... Error! basico: zero area\n";
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
        using std::sqrt;
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

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::tria3d(bool flag, doublereal *_xl, doublereal *_yl, doublereal *_zl,
         doublereal e, doublereal nu, doublereal *_h, doublereal *_rk)
{
/*
c-------------------------------------------------------------------*
c This subroutine evaluates the stiffness matrix and mass vector
c for the spatial 18 d.o.f tree node triangle obtained as
c a combination of the aqr bending triangle plus the membrane
c with drilling d.o.f. developed by Felippa et al.
c
c        rkb will be used for the basic stiffness 
c	 rkm will be used for the higher order stiffness
c
c input variables:
c 	xl = x coordinates
c 	yl = y coordinates
c 	zl = z coordinates
c 	e  = elastic modulus
c 	nu = poisson's ratio
c 	h  = thickness
c
c output variables:
c	rk = stiffness matrix
c 
*-------------------------------------------------------------------*

*
*  SUBROUTINES CALLED: BASICO
*                      SM3MB
*                      SMCBH
*                      SM3MHEFF
*                      ROTATION
*                      TRIROT
*-------------------------------------------------------------------*/

        using std::sqrt;

        Eigen::Map<Eigen::Matrix<doublereal,18,18,Eigen::RowMajor> > rk(_rk);
        Eigen::Map<Eigen::Matrix<doublereal,3,1> > xl(_xl), yl(_yl), zl(_zl), h(_h);
        Eigen::Matrix<doublereal,3,1> xp, yp, zp, xlp, ylp, zlp, v1n, v2n, v3n;
        Eigen::Matrix<doublereal,3,3> db, dm, r1;
        Eigen::Matrix<doublereal,18,18> rkm, rkb;
        doublereal cb,x21,y21,z21,x32,y32,z32,x13,y13,z13;
        doublereal rlr,rlb,bpr,area,ylr,zlr;
        doublereal ycg,xcg,zcg,xlcg,ylcg,zlcg,esp;
        const doublereal f(1.0);
        const doublereal clr(0);
        const doublereal cqr(1.0);
        const doublereal alpha(1.5);
        Eigen::Matrix<int,9,1> lb,le;
        int i,j;
        std::string status;
        lb << 2,3,4,8,9,10,14,15,16;
        le << 0,1,6,7,12,13,5,11,17;
        v1n << 1.0,0.0,0.0;
        v2n << 0.0,1.0,0.0;
        v3n << 0.0,0.0,1.0;

// flag = 0 do NOT perform transform from local to global
// flag = 1 perform transformation from local to global

// set thickness
        esp = h[0];

// bending elastic matrix
        cb      = e*esp*esp*esp/12.0/(1.0-(nu*nu));
        db(0,0) =    cb;
        db(0,1) = nu*cb;
        db(0,2) = 0.0;
        db(1,0) = db(0,1);
        db(1,1) = cb;
        db(1,2) = 0.0;
        db(2,0) = 0.0;
        db(2,1) = 0.0;
        db(2,2) = ((1.0-nu)/2.0)*cb;

// membrane elastic matrix
        cb=e*esp/(1.0-(nu*nu));
        dm(0,0) =    cb;
        dm(0,1) = nu*cb;
        dm(0,2) = 0.0;
        dm(1,0) = dm(0,1);
        dm(1,1) = cb;
        dm(1,2) = 0.0;
        dm(2,0) = 0.0;
        dm(2,1) = 0.0;
        dm(2,2) = ((1.0-nu)/2.0)*cb;

// dimension variables
        x21 = xl[1] - xl[0];
        y21 = yl[1] - yl[0];
        z21 = zl[1] - zl[0];
        x32 = xl[2] - xl[1];
        y32 = yl[2] - yl[1];
        z32 = zl[2] - zl[1];
        x13 = xl[0] - xl[2];
        y13 = yl[0] - yl[2];
        z13 = zl[0] - zl[2];
// triangle in space : we compute the length of one side and the distance of the
// opposing node to that side to compute the area
        rlr = sqrt( x21*x21 + y21*y21 + z21*z21 );
        rlb = sqrt( x32*x32 + y32*y32 + z32*z32 );
        bpr = sqrt((x21 * x32 + y21 * y32 + z21 *z32 )*(x21 * x32 + y21 * y32 + z21 *z32 ))/rlr;
        area= rlr*(sqrt(rlb*rlb-bpr*bpr))/2.0;
// direction cosines of the local system . X' is directed parallel 
// to the 2-1 side
// Z' is the external normal (counterclockwise). Y' computed as Z' x X'
        xp[0] = x21/rlr;
        xp[1] = y21/rlr;
        xp[2] = z21/rlr;
// Z'
        zp[0] = y21 * z32 - z21 * y32;
        zp[1] = z21 * x32 - x21 * z32;
        zp[2] = x21 * y32 - y21 * x32;
        zlr   = sqrt( zp[0]*zp[0] + zp[1]*zp[1]+ zp[2]*zp[2] );
        zp[0] = zp[0]/zlr;
        zp[1] = zp[1]/zlr;
        zp[2] = zp[2]/zlr;
// Y'
        yp[0] = zp[1] * xp[2] - zp[2] * xp[1];
        yp[1] = zp[2] * xp[0] - zp[0] * xp[2];
        yp[2] = zp[0] * xp[1] - zp[1] * xp[0];
        ylr   = sqrt( yp[0]*yp[0] + yp[1]*yp[1] + yp[2]*yp[2] );
        yp[0] = yp[0]/ylr;
        yp[1] = yp[1]/ylr;
        yp[2] = yp[2]/ylr;
// compute center of gravity
        xcg = (xl[0] + xl[1] + xl[2])/3.0;
        ycg = (yl[0] + yl[1] + yl[2])/3.0;
        zcg = (zl[0] + zl[1] + zl[2])/3.0;
// compute local coordinates 
        for (i=0; i<3; ++i) {
          xlcg   = xl[i] - xcg;
          ylcg   = yl[i] - ycg;
          zlcg   = zl[i] - zcg;
          xlp(i) = xp[0] * xlcg + xp[1] * ylcg + xp[2] * zlcg;
          ylp(i) = yp[0] * xlcg + yp[1] * ylcg + yp[2] * zlcg;
          zlp(i) = zp[0] * xlcg + zp[1] * ylcg + zp[2] * zlcg;
        }

// zero stiffness matrices
        rk.setZero();
        rkm.setZero();
        rkb.setZero();

// form local basic bending stiffness
        
        basico(xlp.data(),ylp.data(),db.data(),1.0,clr,cqr,lb.data(),rkb.data(),18,status);

// form local basic membrane stiffness
        sm3mb(xlp.data(),ylp.data(),dm.data(),alpha,1.0,le.data(),rkm.data(),18,status);

// form local higher order bending stiffness
        smcbh(xlp.data(),ylp.data(),db.data(),1.0,lb.data(),rkb.data(),18,status);

// form local higher order membrane stiffness
        sm3mhe(xlp.data(),ylp.data(),dm.data(),0.32,le.data(),rkm.data(),18,status);

// add bending stiffness and membrane stiffness
        rk = rkb + rkm;

// rotate stiffness matrix from local coordinate system to global
// coordinate system in the case of linear FEM. In the case of
// nonlinear FEM with corotational method, do not perform this 
// transformation as the corotational routines expect a stiffness
// matrix in local coordinates.
      if(flag) {
//     compute local to global rotation matrix
        rotation(xp.data(),yp.data(),zp.data(),v1n.data(),v2n.data(),v3n.data(),r1.data());
        trirotation(rk.data(),r1.data());
      }
}

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::basico(doublereal *_x, doublereal *_y, doublereal *_db, doublereal f, 
         doublereal clr, doublereal cqr, int *_ls, doublereal *_sm, int m, std::string &status)
{
      using std::sqrt;

      Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y);
      Eigen::Map<Eigen::Matrix<doublereal,3,3> > db(_db);
      Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > sm(_sm,m,m);
      Eigen::Map<Eigen::Matrix<int,9,1> > ls(_ls);
      
//    t y p e   &   d i m e n s i o n
      Eigen::Matrix<doublereal,9,3> llr,lqr,l;
      doublereal         db11, db12, db13, db22, db23, db33;
      doublereal         x0, y0, cab, a1, a2, a3, b1, b2, b3;
      doublereal         x21, x32, x13, y21, y32, y13;
      doublereal         x12, x23, x31, y12, y23, y31;
      doublereal         xl12, xl23, xl31, c12, c23, c31, s12, s23, s31;
      doublereal         cc12, cc23, cc31, ss12, ss23, ss31;
      doublereal         cs12, cs23, cs31, s1, s2, s3;
      doublereal         area2, c;
      int                i, j, ii, jj;

//    l o g i c
      status =   " ";
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
        status = "basico: negative area";
        if (area2 == 0.0)   status = "basico: zero area";
        return;
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

      llr.setZero();
      lqr.setZero();

      if (clr != 0.0) {
        llr(2,0) =   y32*.5;
        llr(5,0) =   y13*.5;
        llr(8,0) =   y21*.5;
        llr(1,1) =   x32*.5;
        llr(4,1) =   x13*.5;
        llr(7,1) =   x21*.5;
        llr(1,2) =  -y32*.5;
        llr(2,2) =  -x32*.5;
        llr(4,2) =  -y13*.5;
        llr(5,2) =  -x13*.5;
        llr(7,2) =  -y21*.5;
        llr(8,2) =  -x21*.5;
      }

      if (cqr != 0.0) {
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
      }

      for(j=0; j<9; ++j) {
        l(j,0) =   clr*llr(j,0) + cqr*lqr(j,0);
        l(j,1) =   clr*llr(j,1) + cqr*lqr(j,1);
        l(j,2) =   clr*llr(j,2) + cqr*lqr(j,2);
      }

      c =        2.0*f/area2;
      db11 =     c*db(0,0);
      db22 =     c*db(1,1);
      db33 =     c*db(2,2);
      db12 =     c*db(0,1);
      db13 =     c*db(0,2);
      db23 =     c*db(1,2);
     
      for(j=0; j<9; ++j) {
        jj =     ls[j];
        s1 =     db11*l(j,0) + db12*l(j,1) + db13*l(j,2);
        s2 =     db12*l(j,0) + db22*l(j,1) + db23*l(j,2);
        s3 =     db13*l(j,0) + db23*l(j,1) + db33*l(j,2);
        for(i=0; i<=j; ++i) {
          ii =      ls[i];
          sm(jj,ii) =  sm(jj,ii) + (s1*l(i,0)+s2*l(i,1)+s3*l(i,2));
          sm(ii,jj) =  sm(jj,ii);
        }
      }
}

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::sm3mb(doublereal *_x, doublereal *_y, doublereal *_dm, 
        doublereal alpha, doublereal f, int *_ls, doublereal *_sm, int m, std::string &status)
{
      Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y);
      Eigen::Map<Eigen::Matrix<doublereal,3,3> > dm(_dm);
      Eigen::Map<Eigen::Matrix<int,9,1> > ls(_ls);
      Eigen::Matrix<doublereal,9,3> p(9,3);
      Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > sm(_sm,m,m);
      doublereal   area2, c;
      doublereal   d11, d12, d13, d22, d23, d33;
      doublereal   x21, x32, x13, y21, y32, y13;
      doublereal   x12, x23, x31, y12, y23, y31;
      doublereal   s1, s2, s3;

      int          i, j, k, l, n;

//                   L O G I C

      status =   " ";
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
        status = "SM3MB: Zero area";
        if (area2 == 0.0)   status = "SM3MB: Zero area";
        return;
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
      n =        6;
      if (alpha != 0.0) {
        doublereal coef1  = alpha/6.0;
        doublereal coef2  = alpha/3.0;
        p(6,0) =  y23*(y13-y21)*coef1;
        p(6,1) =  x32*(x31-x12)*coef1;
        p(6,2) =  (x31*y13-x12*y21)*coef2;
        p(7,0) =  y31*(y21-y32)*coef1;
        p(7,1) =  x13*(x12-x23)*coef1;
        p(7,2) =  (x12*y21-x23*y32)*coef2;
        p(8,0) =  y12*(y32-y13)*coef1;
        p(8,1) =  x21*(x23-x31)*coef1;
        p(8,2) =  (x23*y32-x31*y13)*coef2;
        n = 9;
      }
      c =       0.5*f/area2;
      d11 =     c * dm(0,0);
      d22 =     c * dm(1,1);
      d33 =     c * dm(2,2);
      d12 =     c * dm(0,1);
      d13 =     c * dm(0,2);
      d23 =     c * dm(1,2);
      for(j = 0; j < n; ++j) {
        l =      ls[j];
        s1 =     d11*p(j,0) + d12*p(j,1) + d13*p(j,2);
        s2 =     d12*p(j,0) + d22*p(j,1) + d23*p(j,2);
        s3 =     d13*p(j,0) + d23*p(j,1) + d33*p(j,2);
        for(i = 0; i <= j; ++i) {
          k =      ls[i];
          sm(k,l) =  sm(k,l) + (s1*p(i,0) + s2*p(i,1) + s3*p(i,2));
          sm(l,k) =  sm(k,l);
        }
      }
      return;
}

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::smcbh(doublereal *_x, doublereal *_y, doublereal *_db, 
        doublereal f, int *_ls, doublereal *_sm, int m, std::string &status)
{
      using std::sqrt;

      Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y);
      Eigen::Map<Eigen::Matrix<int,9,1> > ls(_ls);
      Eigen::Map<Eigen::Matrix<doublereal,3,3> > db(_db);
      Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > sm(_sm,m,m);
      Eigen::Matrix<doublereal,3,3> pg,sq,sds;
      Eigen::Matrix<doublereal,6,6> rsd;
      Eigen::Matrix<doublereal,6,9> q;
      Eigen::Matrix<doublereal,3,6> rm;
      doublereal l1, l2, l3;
      doublereal x0, y0, x1, x2, x3, y1, y2, y3;
      doublereal x21, x32, x13, y21, y32, y13, area, area2;
      doublereal l21, l32, l13, bl2, al2, bl3, al3, bl1, al1;
      doublereal cc, d11, d22, d33, d12, d13, d23;
      doublereal s1, s2, s3, s4, s5, s6;
      int    i, j, k, l;
      pg << 0,0.5,0.5,0.5,0.,0.5,0.5,0.5,0.;

      status =   " ";
      rm.setZero();
      q.setZero();
      rsd.setZero();
      sds.setZero();
      x0=(x[0]+x[1]+x[2])/3.;
      y0=(y[0]+y[1]+y[2])/3.;
      x1= x[0] -x0;
      x2= x[1] -x0;
      x3= x[2] -x0;
      y1= y[0] -y0;
      y2= y[1] -y0;
      y3= y[2] -y0;
      x21 = x2 - x1;
      x32 = x3 - x2;
      x13 = x1 - x3;
      y21 = y2 - y1;
      y32 = y3 - y2;
      y13 = y1 - y3;
      area2 = y21*x13 - x21*y13;
      if (area2 <= 0.0) {
        status = "nega_area";
        if (area2 == 0.0)   status = "zero_area";
        return;
      }
//          side lenghts
      l21 =  sqrt( x21*x21+y21*y21 );
      l32 =  sqrt( x32*x32+y32*y32 );
      l13 =  sqrt( x13*x13+y13*y13 );
//          side proyections
      bl2=((x3-x1)*(x2-x1)+(y3-y1)*(y2-y1))/(l21*l21);
      al2=1.0-bl2;
      bl3=((x1-x2)*(x3-x2)+(y1-y2)*(y3-y2))/(l32*l32);
      al3=1.0-bl3;
      bl1=((x2-x3)*(x1-x3)+(y2-y3)*(y1-y3))/(l13*l13);
      al1=1.0-bl1;
//          inverse of the matrix relating inside curvatures
//          xx,yy,xy with boundary curvatures
      cc=area2*area2*area2;
      sq(0,0)=( -x21*y21*y32*y32 + x32*y21*y21*y32)/cc;
      sq(0,1)=(  x13*y13*y32*y32 - x32*y13*y13*y32)/cc;
      sq(0,2)=(  x21*y21*y13*y13 - x13*y21*y21*y13)/cc;
      sq(1,0)=(  x21*x32*x32*y21 - x21*x21*x32*y32)/cc;
      sq(1,1)=( -x13*x32*x32*y13 + x13*x13*x32*y32)/cc;
      sq(1,2)=( -x21*x13*x13*y21 + x21*x21*x13*y13)/cc;
      sq(2,0)=(  x21*x21*y32*y32  - x32*x32*y21*y21)/cc;
      sq(2,1)=( -x13*x13*y32*y32  + x32*x32*y13*y13)/cc;
      sq(2,2)=( -x21*x21*y13*y13  + x13*x13*y21*y21)/cc;
      d11 =      db(0,0);
      d22 =      db(1,1);
      d33 =      db(2,2);
      d12 =      db(0,1);
      d13 =      db(0,2);
      d23 =      db(1,2);
      area =     0.5*area2;
      for(j=0; j<3; j++) {
        s1=   d11*sq(0,j) + d12*sq(1,j) + d13*sq(2,j);
        s2=   d12*sq(0,j) + d22*sq(1,j) + d23*sq(2,j);
        s3=   d13*sq(0,j) + d23*sq(1,j) + d33*sq(2,j);
        for(i=0; i<=j; ++i) {
          sds(i,j)= sds(i,j)+ (s1*sq(0,i)+s2*sq(1,i)+s3*sq(2,i));
          sds(j,i)= sds(i,j);
        }
      }
      for(j=0; j<3; ++j)
        for(i=0; i<3; ++i)
          sds(i,j)=sds(i,j)*area/3.0*f;
//    matrix q
      q(0,0)=  6.0;
      q(0,1)= -2.0*y13;
      q(0,2)=  2.0*x13;
      q(0,6)= -6.0;
      q(0,7)= -4.0*y13;
      q(0,8)=  4.0*x13;
      q(1,0)= -6.0;
      q(1,1)=  4.0*y13;
      q(1,2)= -4.0*x13;
      q(1,6)=  6.0;
      q(1,7)=  2.0*y13;
      q(1,8)= -2.0*x13;
      q(2,0)= -6.0;
      q(2,1)= -4.0*y21;
      q(2,2)=  4.0*x21;
      q(2,3)=  6.0;
      q(2,4)= -2.0*y21;
      q(2,5)=  2.0*x21;
      q(3,0)=  6.0;
      q(3,1)=  2.0*y21;
      q(3,2)= -2.0*x21;
      q(3,3)= -6.0;
      q(3,4)=  4.0*y21;
      q(3,5)= -4.0*x21;
      q(4,3)= -6.0;
      q(4,4)= -4.0*y32;
      q(4,5)=  4.0*x32;
      q(4,6)=  6.0;
      q(4,7)= -2.0*y32;
      q(4,8)=  2.0*x32;
      q(5,3)=  6.0;
      q(5,4)=  2.0*y32;
      q(5,5)= -2.0*x32;
      q(5,6)= -6.0;
      q(5,7)=  4.0*y32;
      q(5,8)= -4.0*x32;
//    numerical integration
      for(k=0; k<3; ++k) {
        l1=pg(k,0);
        l2=pg(k,1);
        l3=pg(k,2);

//    compute rm in the integration point
        rm(0,0)= l3 +al1 *l2 -(1.+al1)/3.;
        rm(0,1)= l1 +bl1 *l2 -(1.+bl1)/3.;
        rm(1,2)= l1 +al2 *l3 -(1.+al2)/3.;
        rm(1,3)= l2 +bl2 *l3 -(1.+bl2)/3.;
        rm(2,4)= l2 +al3 *l1 -(1.+al3)/3.;
        rm(2,5)= l3 +bl3 *l1 -(1.+bl3)/3.;
        for (j=0; j<6; ++j) {
          s1=sds(0,0)*rm(0,j)+sds(0,1)*rm(1,j)+sds(0,2)*rm(2,j);
          s2=sds(1,0)*rm(0,j)+sds(1,1)*rm(1,j)+sds(1,2)*rm(2,j);
          s3=sds(2,0)*rm(0,j)+sds(2,1)*rm(1,j)+sds(2,2)*rm(2,j);
          for( i=0; i<=j; ++i) {
            rsd(i,j)=rsd(i,j)+(s1*rm(0,i)+s2*rm(1,i)+s3*rm(2,i));
            rsd(j,i)=rsd(i,j);
          }
        }
      }
      for (j=0; j<9; ++j) {
        k =ls[j];
        s1=rsd(0,0)*q(0,j)+rsd(0,1)*q(1,j)+rsd(0,2)*q(2,j)+rsd(0,3)*q(3,j)+rsd(0,4)*q(4,j)+rsd(0,5)*q(5,j);
        s2=rsd(1,0)*q(0,j)+rsd(1,1)*q(1,j)+rsd(1,2)*q(2,j)+rsd(1,3)*q(3,j)+rsd(1,4)*q(4,j)+rsd(1,5)*q(5,j);
        s3=rsd(2,0)*q(0,j)+rsd(2,1)*q(1,j)+rsd(2,2)*q(2,j)+rsd(2,3)*q(3,j)+rsd(2,4)*q(4,j)+rsd(2,5)*q(5,j);
        s4=rsd(3,0)*q(0,j)+rsd(3,1)*q(1,j)+rsd(3,2)*q(2,j)+rsd(3,3)*q(3,j)+rsd(3,4)*q(4,j)+rsd(3,5)*q(5,j);
        s5=rsd(4,0)*q(0,j)+rsd(4,1)*q(1,j)+rsd(4,2)*q(2,j)+rsd(4,3)*q(3,j)+rsd(4,4)*q(4,j)+rsd(4,5)*q(5,j);
        s6=rsd(5,0)*q(0,j)+rsd(5,1)*q(1,j)+rsd(5,2)*q(2,j)+rsd(5,3)*q(3,j)+rsd(5,4)*q(4,j)+rsd(5,5)*q(5,j);
        for (i=0; i<=j; ++i) {
          l = ls[i];
          sm(l,k)=sm(l,k)+(s1*q(0,i)+s2*q(1,i)+s3*q(2,i)+s4*q(3,i)+s5*q(4,i)+s6*q(5,i));
          sm(k,l)=sm(l,k);
        }
       }
}

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::sm3mhe(doublereal *_x, doublereal *_y, doublereal *_dm, doublereal f, int *_ls, doublereal *_sm, int m, std::string &status)
{
      Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y);
      Eigen::Map<Eigen::Matrix<int,9,1> > ls(_ls);
      Eigen::Map<Eigen::Matrix<doublereal,3,3> > dm(_dm);
      Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > sm(_sm,m,m);

      doublereal  x0,y0, x10,x20,x30, y10,y20,y30;
      doublereal  x12, x21, x23, x32, x31, x13;
      doublereal  y12, y21, y23, y32, y31, y13;
      doublereal  aa12,aa23,aa31,ss12,ss23,ss31,ss1,ss2,ss3;
      doublereal  caa12,caa23,caa31, sum;
      doublereal  ca,cax10,cax20,cax30,cay10,cay20,cay30;
      doublereal  area, area2, kfac;
      Eigen::Matrix<doublereal,6,6> kqh;
      Eigen::Matrix<doublereal,6,3> hmt, hqt;
      Eigen::Matrix<doublereal,3,3> kth;
      Eigen::Matrix<doublereal,3,1> s;
      Eigen::Matrix<doublereal,6,1> w, xyij;
      doublereal  e11,e22,e33,e12,e13,e23;
      int         i,j,k,l;

//                   L O G I C

      status =   " ";
      if (f == 0.0) return;
      x12 =      x[0] - x[1];
      x21 =     -x12;
      x23 =      x[1] - x[2];
      x32 =     -x23;
      x31 =      x[2] - x[0];
      x13 =     -x31;
      y12 =      y[0] - y[1];
      y21 =     -y12;
      y23 =      y[1] - y[2];
      y32 =     -y23;
      y31 =      y[2] - y[0];
      y13 =     -y31;
      area2 =    x21*y31-x31*y21;
      if (area2 <= 0.0) {
        status = "SM3MBE: Negative area";
        if (area2 == 0.0)   status = "SM3MBE: Zero area";
        return;
      }
      area  =    0.5*area2;
      x0 =       (x[0]+x[1]+x[2])/3.0;
      y0 =       (y[0]+y[1]+y[2])/3.0;
      x10 =      x[0] - x0;
      x20 =      x[1] - x0;
      x30 =      x[2] - x0;
      y10 =      y[0] - y0;
      y20 =      y[1] - y0;
      y30 =      y[2] - y0;
      aa12 =     2.25*(x30*x30+y30*y30);
      aa23 =     2.25*(x10*x10+y10*y10);
      aa31 =     2.25*(x20*x20+y20*y20);
      caa12 =    15./(32.*aa12);
      caa23 =    15./(32.*aa23);
      caa31 =    15./(32.*aa31);
      ss12 =     x12*x12+y12*y12;
      ss23 =     x23*x23+y23*y23;
      ss31 =     x31*x31+y31*y31;
      ss1 =      0.25*(ss12-ss31);
      ss2 =      0.25*(ss23-ss12);
      ss3 =      0.25*(ss31-ss23);
      cay10 =    0.1875*y10;
      cay20 =    0.1875*y20;
      cay30 =    0.1875*y30;
      cax10 =    0.1875*x10;
      cax20 =    0.1875*x20;
      cax30 =    0.1875*x30;
      hmt(0,0) = caa12*((-ss3+0.6*aa12)*y30+area*x30);
      hmt(0,1) =  3.*cay30 - hmt(0,0);
      hmt(0,2) = cay30;
      hmt(1,0) = cay10;
      hmt(1,1) = caa23*((-ss1+0.6*aa23)*y10+area*x10);
      hmt(1,2) =  3.*cay10 - hmt(1,1);
      hmt(2,0) = caa31*((ss2+0.6*aa31)*y20-area*x20);
      hmt(2,1) = cay20;
      hmt(2,2) =  3.*cay20 - hmt(2,0);
      hmt(3,0) = caa12*((ss3-0.6*aa12)*x30+area*y30);
      hmt(3,1) = -3.*cax30 - hmt(3,0);
      hmt(3,2) = -cax30;
      hmt(4,0) = -cax10;
      hmt(4,1) = caa23*((ss1-0.6*aa23)*x10+area*y10);
      hmt(4,2) = -3.*cax10 - hmt(4,1);
      hmt(5,0) = caa31*((-ss2-0.6*aa31)*x20-area*y20);
      hmt(5,1) = -cax20;
      hmt(5,2) = -3.*cax20 - hmt(5,0);
      for(j = 0; j<3; ++j) {
        sum =    (2./9.)*(hmt(0,j)+hmt(1,j)+hmt(2,j));
        hqt(0,j) =  sum - (4./3.)*hmt(0,j);
        hqt(1,j) =  sum - (4./3.)*hmt(1,j);
        hqt(2,j) =  sum - (4./3.)*hmt(2,j);
        sum =    (2./9.)*(hmt(3,j)+hmt(4,j)+hmt(5,j));
        hqt(3,j) =  sum - (4./3.)*hmt(3,j);
        hqt(4,j) =  sum - (4./3.)*hmt(4,j);
        hqt(5,j) =  sum - (4./3.)*hmt(5,j);
      }
      kfac =     1.5*f/area2;
      e11 =      kfac * dm(0,0);
      e22 =      kfac * dm(1,1);
      e33 =      kfac * dm(2,2);
      e12 =      kfac * dm(0,1);
      e13 =      kfac * dm(0,2);
      e23 =      kfac * dm(1,2);
      kqh(0,0) = 2*(e11*y30*y30-2*e13*x30*y30+e33*x30*x30);
      kqh(0,1) = ((e13*x10-e11*y10)*y30+(e13*y10-e33*x10)*x30);
      kqh(0,2) = ((e13*x20-e11*y20)*y30+(e13*y20-e33*x20)*x30);
      kqh(0,3) = 2*(e13*y30*y30-(e33+e12)*x30*y30+e23*x30*x30);
      kqh(0,4) = ((e12*x10-e13*y10)*y30+(e33*y10-e23*x10)*x30);
      kqh(0,5) = ((e12*x20-e13*y20)*y30+(e33*y20-e23*x20)*x30);
      kqh(1,0) = kqh(0,1);
      kqh(1,1) = 2*(e11*y10*y10-2*e13*x10*y10+e33*x10*x10);
      kqh(1,2) = ((e13*x10-e11*y10)*y20+(e13*y10-e33*x10)*x20);
      kqh(1,3) = ((e33*x10-e13*y10)*y30+(e12*y10-e23*x10)*x30);
      kqh(1,4) = 2*(e13*y10*y10-(e33+e12)*x10*y10+e23*x10*x10);
      kqh(1,5) = ((e33*x10-e13*y10)*y20+(e12*y10-e23*x10)*x20);
      kqh(2,0) = kqh(0,2);
      kqh(2,1) = kqh(1,2);
      kqh(2,2) = 2*(e11*y20*y20-2*e13*x20*y20+e33*x20*x20);
      kqh(2,3) = ((e33*x20-e13*y20)*y30+(e12*y20-e23*x20)*x30);
      kqh(2,4) = ((e12*x10-e13*y10)*y20+(e33*y10-e23*x10)*x20);
      kqh(2,5) = 2*(e13*y20*y20-(e33+e12)*x20*y20+e23*x20*x20);
      kqh(3,0) = kqh(0,3);
      kqh(3,1) = kqh(1,3);
      kqh(3,2) = kqh(2,3);
      kqh(3,3) = 2*(e33*y30*y30-2*e23*x30*y30+e22*x30*x30);
      kqh(3,4) = ((e23*x10-e33*y10)*y30+(e23*y10-e22*x10)*x30);
      kqh(3,5) = ((e23*x20-e33*y20)*y30+(e23*y20-e22*x20)*x30);
      kqh(4,0) = kqh(0,4);
      kqh(4,1) = kqh(1,4);
      kqh(4,2) = kqh(2,4);
      kqh(4,3) = kqh(3,4);
      kqh(4,4) = 2*(e33*y10*y10-2*e23*x10*y10+e22*x10*x10);
      kqh(4,5) = ((e23*x10-e33*y10)*y20+(e23*y10-e22*x10)*x20);
      kqh(5,0) = kqh(0,5);
      kqh(5,1) = kqh(1,5);
      kqh(5,2) = kqh(2,5);
      kqh(5,3) = kqh(3,5);
      kqh(5,4) = kqh(4,5);
      kqh(5,5) = 2*(e33*y20*y20-2*e23*x20*y20+e22*x20*x20);
      kth(0,0) =  0.0;
      kth(0,1) =  0.0;
      kth(1,1) =  0.0;
      kth(0,2) =  0.0;
      kth(1,2) =  0.0;
      kth(2,2) =  0.0;
      for(j = 0; j<3; ++j) {
        for(i = 0; i<6; ++i) {
          w[i] =    kqh(i,0)*hqt(0,j) + kqh(i,1)*hqt(1,j)
                   + kqh(i,2)*hqt(2,j) + kqh(i,3)*hqt(3,j)
                   + kqh(i,4)*hqt(4,j) + kqh(i,5)*hqt(5,j);
        }
        for(i = 0; i<=j; ++i) {
          kth(i,j) = kth(i,j) + hqt(0,i)*w[0] + hqt(1,i)*w[1]
                              + hqt(2,i)*w[2] + hqt(3,i)*w[3]
                              + hqt(4,i)*w[4] + hqt(5,i)*w[5];
          kth(j,i) = kth(i,j);
        }
      }
      s[0] =   kth(0,0) + kth(0,1) + kth(0,2);
      s[1] =   kth(1,0) + kth(1,1) + kth(1,2);
      s[2] =   kth(2,0) + kth(2,1) + kth(2,2);
      ca =       0.25/area;
      xyij[0] =  ca*x32;
      xyij[1] =  ca*y32;
      xyij[2] =  ca*x13;
      xyij[3] =  ca*y13;
      xyij[4] =  ca*x21;
      xyij[5] =  ca*y21;
      for(j = 0; j < 9; ++j) {
        l =    ls[j];
        for(i = 0; i<3; ++i) {
          if (j < 6) {
             w[i] =  s[i]*xyij[j];
          } else {
             w[i] =  kth(i,j-6);
          }
        }
        sum =    w[0] + w[1] + w[2];
        for(i = 0; i<=j; ++i) {
          k =      ls[i];
          if (i < 6) {
             sm(k,l) =  sm(k,l) + sum*xyij(i);
          } else {
             sm(k,l) =  sm(k,l) + w(i-6);
          }
          sm(l,k) =  sm(k,l);
        }
      }
}

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::trirotation(doublereal *_sm, doublereal *_r1)
{
      Eigen::Map<Eigen::Matrix<doublereal,3,3> > r1(_r1);
      Eigen::Map<Eigen::Matrix<doublereal,18,18> > sm(_sm);

      Eigen::Matrix<doublereal,18,1> prod1;
      int i, j, k;

// premultiplication by rotation matrix

      for(j=0; j<18; ++j) {
        for(k=0; k<3; ++k) {
          prod1[k]    = 0.0;
          prod1[k+3]  = 0.0;
          prod1[k+6]  = 0.0;
          prod1[k+9]  = 0.0;
          prod1[k+12] = 0.0;
          prod1[k+15] = 0.0;
          for(i=0; i<3; ++i) {
            prod1[k]    = prod1[k]    + r1(k,i)*sm(i   ,j);
            prod1[k+3]  = prod1[k+3]  + r1(k,i)*sm(i+3 ,j);
            prod1[k+6]  = prod1[k+6]  + r1(k,i)*sm(i+6 ,j);
            prod1[k+9]  = prod1[k+9]  + r1(k,i)*sm(i+9 ,j);
            prod1[k+12] = prod1[k+12] + r1(k,i)*sm(i+12,j);
            prod1[k+15] = prod1[k+15] + r1(k,i)*sm(i+15,j);
          }
        }

        for(k=0; k<18; ++k) {
          sm(k,j) = prod1[k];
        }
      }

// postmultiplication by rotation matrix transposed

      for(j=0; j<18; ++j) {
        for(k=0; k<3; ++k) {
          prod1[k]   =0.0;
          prod1[k+3] =0.0;
          prod1[k+6] =0.0;
          prod1[k+9] =0.0;
          prod1[k+12]=0.0;
          prod1[k+15]=0.0;
          for(i=0; i<3; ++i) {
            prod1[k]   =prod1[k]   +sm(j,i)*r1(k,i);
            prod1[k+3] =prod1[k+3] +sm(j,i+3)*r1(k,i);
            prod1[k+6] =prod1[k+6] +sm(j,i+6)*r1(k,i);
            prod1[k+9] =prod1[k+9] +sm(j,i+9)*r1(k,i);
            prod1[k+12]=prod1[k+12]+sm(j,i+12)*r1(k,i);
            prod1[k+15]=prod1[k+15]+sm(j,i+15)*r1(k,i);
          }
        }

        for(k=0; k<18; ++k) sm(j,k)=prod1[k];
      }
}

template<typename doublereal>
void
ShellElementSemiTemplate<doublereal>
::tria3dthickness(bool flag, doublereal *_xl, doublereal *_yl, doublereal *_zl,
                  doublereal e, doublereal nu, doublereal *_h, doublereal *_drkdh)
{
        using std::sqrt;

        Eigen::Map<Eigen::Matrix<doublereal,18,18,Eigen::RowMajor> > drkdh(_drkdh);
        Eigen::Map<Eigen::Matrix<doublereal,3,1> > xl(_xl), yl(_yl), zl(_zl), h(_h);
        Eigen::Matrix<doublereal,3,1> xp, yp, zp, xlp, ylp, zlp, v1n, v2n, v3n;
        Eigen::Matrix<doublereal,3,3> db, dm, r1;
        Eigen::Matrix<doublereal,18,18> drkmdh, drkbdh;
        doublereal cb,x21,y21,z21,x32,y32,z32,x13,y13,z13;
        doublereal rlr,rlb,bpr,area,ylr,zlr;
        doublereal ycg,xcg,zcg,xlcg,ylcg,zlcg,esp;
        const doublereal f(1.0);
        const doublereal clr(0);
        const doublereal cqr(1.0);
        const doublereal alpha(1.5);
        Eigen::Matrix<int,9,1> lb,le;
        int i,j;
        std::string status;
        lb << 2,3,4,8,9,10,14,15,16;
        le << 0,1,6,7,12,13,5,11,17;
        v1n << 1.0,0.0,0.0;
        v2n << 0.0,1.0,0.0;
        v3n << 0.0,0.0,1.0;

// flag = 0 do NOT perform transform from local to global
// flag = 1 perform transformation from local to global

// set thickness
        esp = h[0];

// bending elastic matrix
        cb      = e*esp*esp/4.0/(1.0-(nu*nu));
        db(0,0) =    cb;
        db(0,1) = nu*cb;
        db(0,2) = 0.0;
        db(1,0) = db(0,1);
        db(1,1) = cb;
        db(1,2) = 0.0;
        db(2,0) = 0.0;
        db(2,1) = 0.0;
        db(2,2) = ((1.0-nu)/2.0)*cb;

// membrane elastic matrix
        cb=e/(1.0-(nu*nu));
        dm(0,0) =    cb;
        dm(0,1) = nu*cb;
        dm(0,2) = 0.0;
        dm(1,0) = dm(0,1);
        dm(1,1) = cb;
        dm(1,2) = 0.0;
        dm(2,0) = 0.0;
        dm(2,1) = 0.0;
        dm(2,2) = ((1.0-nu)/2.0)*cb;

// dimension variables
        x21 = xl[1] - xl[0];
        y21 = yl[1] - yl[0];
        z21 = zl[1] - zl[0];
        x32 = xl[2] - xl[1];
        y32 = yl[2] - yl[1];
        z32 = zl[2] - zl[1];
        x13 = xl[0] - xl[2];
        y13 = yl[0] - yl[2];
        z13 = zl[0] - zl[2];
// triangle in space : we compute the length of one side and the distance of the
// opposing node to that side to compute the area
        rlr = sqrt( x21*x21 + y21*y21 + z21*z21 );
        rlb = sqrt( x32*x32 + y32*y32 + z32*z32 );
        bpr = sqrt((x21 * x32 + y21 * y32 + z21 *z32 )*(x21 * x32 + y21 * y32 + z21 *z32 ))/rlr;
        area= rlr*(sqrt(rlb*rlb-bpr*bpr))/2.0;
// direction cosines of the local system . X' is directed parallel 
// to the 2-1 side
// Z' is the external normal (counterclockwise). Y' computed as Z' x X'
        xp[0] = x21/rlr;
        xp[1] = y21/rlr;
        xp[2] = z21/rlr;
// Z'
        zp[0] = y21 * z32 - z21 * y32;
        zp[1] = z21 * x32 - x21 * z32;
        zp[2] = x21 * y32 - y21 * x32;
        zlr   = sqrt( zp[0]*zp[0] + zp[1]*zp[1]+ zp[2]*zp[2] );
        zp[0] = zp[0]/zlr;
        zp[1] = zp[1]/zlr;
        zp[2] = zp[2]/zlr;
// Y'
        yp[0] = zp[1] * xp[2] - zp[2] * xp[1];
        yp[1] = zp[2] * xp[0] - zp[0] * xp[2];
        yp[2] = zp[0] * xp[1] - zp[1] * xp[0];
        ylr   = sqrt( yp[0]*yp[0] + yp[1]*yp[1] + yp[2]*yp[2] );
        yp[0] = yp[0]/ylr;
        yp[1] = yp[1]/ylr;
        yp[2] = yp[2]/ylr;
// compute center of gravity
        xcg = (xl[0] + xl[1] + xl[2])/3.0;
        ycg = (yl[0] + yl[1] + yl[2])/3.0;
        zcg = (zl[0] + zl[1] + zl[2])/3.0;
// compute local coordinates 
        for (i=0; i<3; ++i) {
          xlcg   = xl[i] - xcg;
          ylcg   = yl[i] - ycg;
          zlcg   = zl[i] - zcg;
          xlp(i) = xp[0] * xlcg + xp[1] * ylcg + xp[2] * zlcg;
          ylp(i) = yp[0] * xlcg + yp[1] * ylcg + yp[2] * zlcg;
          zlp(i) = zp[0] * xlcg + zp[1] * ylcg + zp[2] * zlcg;
        }

// zero stiffness matrices
        drkdh.setZero();
        drkmdh.setZero();
        drkbdh.setZero();

// form local basic bending stiffness
        
        basico(xlp.data(),ylp.data(),db.data(),1.0,clr,cqr,lb.data(),drkbdh.data(),18,status);

// form local basic membrane stiffness
        sm3mb(xlp.data(),ylp.data(),dm.data(),alpha,1.0,le.data(),drkmdh.data(),18,status);

// form local higher order bending stiffness
        smcbh(xlp.data(),ylp.data(),db.data(),1.0,lb.data(),drkbdh.data(),18,status);

// form local higher order membrane stiffness
        sm3mhe(xlp.data(),ylp.data(),dm.data(),0.32,le.data(),drkmdh.data(),18,status);

// add bending stiffness and membrane stiffness
        drkdh = drkbdh + drkmdh;

// rotate stiffness matrix from local coordinate system to global
// coordinate system in the case of linear FEM. In the case of
// nonlinear FEM with corotational method, do not perform this 
// transformation as the corotational routines expect a stiffness
// matrix in local coordinates.
      if(flag) {
//     compute local to global rotation matrix
        rotation(xp.data(),yp.data(),zp.data(),v1n.data(),v2n.data(),v3n.data(),r1.data());
        trirotation(drkdh.data(),r1.data());
      }
}

#endif
#endif
