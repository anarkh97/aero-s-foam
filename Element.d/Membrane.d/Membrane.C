#include <cstdio>
#include <cmath>

#include <Element.d/Membrane.d/Membrane.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Corotational.d/utilities.h>
#include <Element.d/Membrane.d/MembraneStressWRTThicknessSensitivity.h>
#include <Element.d/Membrane.d/MembraneStressWRTDisplacementSensitivity.h>
#include <Element.d/Membrane.d/MembraneStiffnessWRTThicknessSensitivity.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/Vector.h>
#include <Corotational.d/Shell3Corotator.h>

extern int verboseFlag;
extern "C"      {
void _FORTRAN(trimem)(int&, double* ,double* ,double* ,double& , double& ,
                      double* ,double* );
void _FORTRAN(sands19)(double*,double*,double*,double&,double&,double*,double*,
        double*, const int&, const int&, const int&, const int&, const int&);
}

Membrane::Membrane(int* nodeNumbers)
{
	nn[0] = nodeNumbers[0];
	nn[1] = nodeNumbers[1];
	nn[2] = nodeNumbers[2];
}

Element *
Membrane::clone()
{
	return new Membrane(*this);
}

void
Membrane::renum(int *renumberingTable)
{
	nn[0] = renumberingTable[nn[0]];
	nn[1] = renumberingTable[nn[1]];
	nn[2] = renumberingTable[nn[2]];
}

void
Membrane::renum(EleRenumMap& renumberingTable)
{
	nn[0] = renumberingTable[nn[0]];
	nn[1] = renumberingTable[nn[1]];
	nn[2] = renumberingTable[nn[2]];
}

void
Membrane::getVonMises(Vector& stress,Vector& weight,CoordSet &cs, Vector& elDisp, 
                      int strInd,int,double*,double ylayer, double zlayer, int avgnum)
{
	weight = 1.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        h[0] = h[1] = h[2] = prop->eh;

        double elStress[3][7];

        int  msize = 1;
        int maxstr = 7;
        int maxgus = 3;
        int    elm = 1; 

	// SET STRAIN FLAG IF USER WANTS STRAIN OUTPUT 
	int strainFlg;
	if(strInd > 6) strainFlg = 1; 

#ifdef USE_EIGEN3
       sands19(x,y,z,prop->E,prop->nu,h,elDisp.data(),(double*)elStress, strainFlg);
#else
       _FORTRAN(sands19)(x,y,z,prop->E,prop->nu,h,elDisp.data(),(double*)elStress,
                         msize,maxstr,maxgus,elm,strainFlg);
#endif

	if(strInd < 7) {
          stress[0] = elStress[0][strInd];
          stress[1] = elStress[1][strInd];
          stress[2] = elStress[2][strInd];
	}
	else {
          stress[0] = elStress[0][strInd-7];
          stress[1] = elStress[1][strInd-7];
          stress[2] = elStress[2][strInd-7];
	}

}

void 
Membrane::getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                             int senMethod, double *ndTemps, int avgnum, double ylayer, double zlayer)
{
#ifdef USE_EIGEN3
   if(strInd != 6) {
     std::cerr << " ... Error: strInd must be 6 in TwoNodeTruss::getVonMisesDisplacementSensitivity\n";
     exit(-1);
   }
   if(dStdDisp.numRow() != 3 || dStdDisp.numCol() !=18) {
     std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
     exit(-1);
   }
   if(ndTemps != 0) {
     std::cerr << " ... Error: thermal stress should not be passed in sensitivity computation\n";
     exit(-1);
   }
  weight = 1;
  // scalar parameters
  Eigen::Array<double,12,1> dconst;

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);

  dconst[0] = nd1.x; dconst[1] = nd2.x; dconst[2] = nd3.x; // x coordinates
  dconst[3] = nd1.y; dconst[4] = nd2.y; dconst[5] = nd3.y; // y coordinates
  dconst[6] = nd1.z; dconst[7] = nd2.z; dconst[8] = nd3.z; // z coordinates
  dconst[9] = prop->E;     // E
  dconst[10] = prop->nu;   // nu
  dconst[11] = prop->eh;   // thickness

  Eigen::Array<double,3,1> globalx = dconst.segment<3>(0).cast<double>();
  Eigen::Array<double,3,1> globaly = dconst.segment<3>(3).cast<double>();
  Eigen::Array<double,3,1> globalz = dconst.segment<3>(6).cast<double>();
   
  // integer parameters
  Eigen::Array<int,1,1> iconst;
  iconst[0] = surface; // surface
  // inputs
  Eigen::Matrix<double,18,1> q = Eigen::Map<Eigen::Matrix<double,18,1> >(elDisp.data()).segment(0,18); //displacements

  // Jacobian evaluation
  Eigen::Matrix<double,3,18> dStressdDisp;
  Eigen::Matrix<double,7,3> stress;
  if(verboseFlag) std::cerr << "senMethod is " << senMethod << std::endl;
 
  if(avgnum == 0 || avgnum == 1) { // NODALFULL or ELEMENTAL
    if(senMethod == 0) { // analytic
/*      dStressdDisp.setZero();
      double h[3]; 
      h[0] = h[1] = h[2] = prop->eh;
      vms8WRTdisp(globalx.data(), globaly.data(), globalz.data(),
                  prop->E, prop->nu, h, q.data(), 
                  dStressdDisp.data(), 0, surface, 0);   
      dStdDisp.copy(dStressdDisp.data());
      if(verboseFlag) std::cerr << "dStressdDisp(analytic) =\n" << dStressdDisp << std::endl;
*/
    }

    if(senMethod == 1) { // automatic differentiation
      Simo::Jacobian<double,MembraneStressWRTDisplacementSensitivity> dSdu(dconst,iconst);
      dStressdDisp = dSdu(q, 0);
      dStdDisp.copy(dStressdDisp.data());
      if(verboseFlag) std::cerr << "dStressdDisp(AD) =\n" << dStressdDisp << std::endl;
    }
 

    if(senMethod == 2) { // finite difference
      // finite difference
      dStressdDisp.setZero();
      MembraneStressWRTDisplacementSensitivity<double> foo(dconst,iconst);
      Eigen::Matrix<double,18,1> qp, qm;
      double h(1e-6);
      Eigen::Matrix<double,3,1> S = foo(q,0);
//      cout << "displacement = " << q.transpose() << endl;
      for(int i=0; i<18; ++i) {
        qp = q;             qm = q;
        qp[i] += h;         qm[i] -= h;
//        if(q[i] == 0) { qp[i] = h;   qm[i] = -h; }
//        else { qp[i] = q[i]*(1 + h);   qm[i] = q[i]*(1 - h); }
        Eigen::Matrix<double,3,1> Sp = foo(qp, 0);
        Eigen::Matrix<double,3,1> Sm = foo(qm, 0);
        Eigen::Matrix<double,3,1> fd = (Sp - Sm)/(2*h);
//        if(i==2) {
//          Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, " "); 
//          cout << Sp.transpose().format(HeavyFmt) << endl;
//          cout << Sm.transpose().format(HeavyFmt) << endl;
//          cout << fd.transpose().format(HeavyFmt) << endl;
//        }
        for(int j=0; j<3; ++j) {
          dStressdDisp(j,i) = fd[j];
        }
      }
      dStdDisp.copy(dStressdDisp.data());
      if(verboseFlag) std::cerr << "dStressdDisp(FD) =\n" << dStressdDisp << std::endl;
    }
  } else dStdDisp.zero(); // NODALPARTIAL or GAUSS or any others
#else
  std::cerr << " ... Error! Membrane::getVonMisesDisplacementSensitivity needs Eigen library.\n";
  exit(-1);
#endif
}

void 
Membrane::getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                          int senMethod, double *, int avgnum, double ylayer, double zlayer)
{
#ifdef USE_EIGEN3
   if(strInd != 6) {
     std::cerr << " ... Error: strInd must be 6 in TwoNodeTruss::getVonMisesThicknessSensitivity\n";
     exit(-1);
   }
   if(dStdThick.size() !=3) {
     std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
     exit(-1);
   }
  weight = 1;
  // scalar parameters
  Eigen::Array<double,29,1> dconst;

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);

  dconst[0] = nd1.x; dconst[1] = nd2.x; dconst[2] = nd3.x; // x coordinates
  dconst[3] = nd1.y; dconst[4] = nd2.y; dconst[5] = nd3.y; // y coordinates
  dconst[6] = nd1.z; dconst[7] = nd2.z; dconst[8] = nd3.z; // z coordinates
  for(int i=0; i<18; ++i) {
    dconst[i+9] = elDisp[i];
  }
  dconst[27] = prop->E;     // E
  dconst[28] = prop->nu;   // nu

  Eigen::Array<double,3,1> globalx = dconst.segment<3>(0).cast<double>();
  Eigen::Array<double,3,1> globaly = dconst.segment<3>(3).cast<double>();
  Eigen::Array<double,3,1> globalz = dconst.segment<3>(6).cast<double>();
   
  // integer parameters
  Eigen::Array<int,1,1> iconst;
  iconst[0] = surface; // surface
  // inputs
  Eigen::Matrix<double,1,1> q;
  q[0] = prop->eh;

  // Jacobian evaluation
  Eigen::Matrix<double,3,1> dStressdThic;
  Eigen::Matrix<double,7,3> stress;
  if(verboseFlag) std::cerr << "senMethod is " << senMethod << std::endl;
 
  if(avgnum == 0 || avgnum == 1) { // NODALFULL or ELEMENTAL
    if(senMethod == 0) { // analytic
/*      dStressdThic.setZero();
      double h[3]; 
      h[0] = h[1] = h[2] = q[0];
      vms8WRTthic(globalx.data(), globaly.data(), globalz.data(),
                  prop->E, prop->nu, h, elDisp.data(), 
                  dStressdThic.data(), 0, surface, 0);   
      dStdThick.copy(dStressdThic.data());
      if(verboseFlag) std::cerr << "dStressdThic(analytic) =\n" << dStressdThic << std::endl;
*/
    }

    if(senMethod == 1) { // automatic differentiation
      Simo::Jacobian<double,MembraneStressWRTThicknessSensitivity> dSdu(dconst,iconst);
      dStressdThic = dSdu(q, 0);
      dStdThick.copy(dStressdThic.data());
      if(verboseFlag) std::cerr << "dStressdThic(AD) =\n" << dStressdThic << std::endl;
    }
 

    if(senMethod == 2) { // finite difference
      // finite difference
      dStressdThic.setZero();
      MembraneStressWRTThicknessSensitivity<double> foo(dconst,iconst);
      Eigen::Matrix<double,1,1> qp, qm;
      double h(1e-6);
      Eigen::Matrix<double,3,1> S = foo(q,0);
//      cout << "displacement = " << q.transpose() << endl;
      qp = q;             qm = q;
      qp[0] += h;         qm[0] -= h;
//        if(q[i] == 0) { qp[i] = h;   qm[i] = -h; }
//        else { qp[i] = q[i]*(1 + h);   qm[i] = q[i]*(1 - h); }
      Eigen::Matrix<double,3,1> Sp = foo(qp, 0);
      Eigen::Matrix<double,3,1> Sm = foo(qm, 0);
      Eigen::Matrix<double,3,1> fd = (Sp - Sm)/(2*h);
//        if(i==2) {
//          Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, " "); 
//          cout << Sp.transpose().format(HeavyFmt) << endl;
//          cout << Sm.transpose().format(HeavyFmt) << endl;
//          cout << fd.transpose().format(HeavyFmt) << endl;
//        }
      for(int j=0; j<3; ++j) {
        dStressdThic(j,0) = fd[j];
      }
      dStdThick.copy(dStressdThic.data());
      if(verboseFlag) std::cerr << "dStressdThic(FD) =\n" << dStressdThic << std::endl;
    }
  } else dStdThick.zero(); // NODALPARTIAL or GAUSS or any others
#else
  std::cerr << " ... Error! Membrane::getVonMisesThicknessSensitivity needs Eigen library\n";
  exit(-1);
#endif
}

void
Membrane::getAllStress(FullM& stress,Vector& weight,CoordSet &cs, Vector& elDisp, int strInd,int,double*)
{
	weight = 1.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        h[0] = h[1] = h[2] = prop->eh;

        double elStress[3][7];

        int  msize = 1;
        int maxstr = 7;
        int maxgus = 3;
        int    elm = 1; 

       _FORTRAN(sands19)(x,y,z,prop->E,prop->nu,h,elDisp.data(),(double*)elStress,
                         msize,maxstr,maxgus,elm,strInd);

// Store all Stress or all Strain as defined by strInd
        int i,j;
        for (i=0; i<3; ++i) {
          for (j=0; j<6; ++j) {
            stress[i][j] = elStress[i][j];
          }
        }

// Get Element Principals
        double svec[6], pvec[3];
        for (j=0; j<6; ++j) {
          svec[j] = stress[0][j];
        }
// Convert Engineering to Tensor Strains
        if(strInd != 0) {
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
        }
        pstress(svec,pvec);
        for (i=0; i<3; ++i) {
          for (j=0; j<3; ++j) {
            stress[i][j+6] = pvec[j];
          }
        }
}

double
Membrane::getMass(CoordSet& cs)
{
        if (prop == NULL) return 0.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double r1[3], r2[3], r3[3], v1[3], v2[3], v3[3];

        r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
        r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
        r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

        v1[0] = r3[0] - r1[0];
        v1[1] = r3[1] - r1[1];
        v1[2] = r3[2] - r1[2];

        v2[0] = r2[0] - r1[0];
        v2[1] = r2[1] - r1[1];
        v2[2] = r2[2] - r1[2];

        v3[0] =  v1[1]*v2[2] - v2[1]*v1[2];
        v3[1] = -v1[0]*v2[2] + v2[0]*v1[2];
        v3[2] =  v1[0]*v2[1] - v2[0]*v1[1];

        double area = 0.5*sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
        double mass = area*prop->rho*prop->eh;

        return mass;

}

double
Membrane::getMassSensitivityWRTthickness(CoordSet& cs)
{
        if (prop == NULL) return 0.0;

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double r1[3], r2[3], r3[3], v1[3], v2[3], v3[3];

        r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
        r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
        r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

        v1[0] = r3[0] - r1[0];
        v1[1] = r3[1] - r1[1];
        v1[2] = r3[2] - r1[2];

        v2[0] = r2[0] - r1[0];
        v2[1] = r2[1] - r1[1];
        v2[2] = r2[2] - r1[2];

        v3[0] =  v1[1]*v2[2] - v2[1]*v1[2];
        v3[1] = -v1[0]*v2[2] + v2[0]*v1[2];
        v3[2] =  v1[0]*v2[1] - v2[0]*v1[1];

        double area = 0.5*sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
        double massWRTthic = area*prop->rho;

        return massWRTthic;
}

double
Membrane::weight(CoordSet& cs, double *gravityAcceleration)
{
  if (prop == NULL) {
    return 0.0;
  }

  double _mass = getMass(cs);
  double gravAccNorm = sqrt(gravityAcceleration[0]*gravityAcceleration[0] + 
                            gravityAcceleration[1]*gravityAcceleration[1] +
                            gravityAcceleration[2]*gravityAcceleration[2]);
  return _mass*gravAccNorm;
}

double
Membrane::weightDerivativeWRTthickness(CoordSet& cs, double *gravityAcceleration, int senMethod)
{
  if (prop == NULL) {
    return 0.0;
  }

  if(senMethod == 0) {
    double _weight = weight(cs, gravityAcceleration);
    double thick = prop->eh;
    return _weight/thick;
  } else {
    fprintf(stderr," ... Error: Membrane::weightDerivativeWRTthickness for automatic differentiation and finite difference is not implemented\n");
    exit(-1);
  }
}

void
Membrane::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                          Vector& gravityForce, int gravflg, GeomState *geomState)
{
        double mass = getMass(cs);
        double massPerNode = mass/3.0;
        double fx, fy, fz;

        // Lumped
        if(gravflg != 2) {

          fx = massPerNode*gravityAcceleration[0];
          fy = massPerNode*gravityAcceleration[1];
          fz = massPerNode*gravityAcceleration[2];

        }
        // Consistent
        else {
          int i;
          Node &nd1 = cs.getNode(nn[0]);
          Node &nd2 = cs.getNode(nn[1]);
          Node &nd3 = cs.getNode(nn[2]);
          double x[3], y[3], z[3], localg[3];
          double T1[3],T2[3],T3[3];

          // Set the coordinates
          x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
          x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
          x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

          // Local X-axis
          T1[0] = x[1] - x[0];
          T1[1] = y[1] - y[0];
          T1[2] = z[1] - z[0];
          normalize( T1 );
          // 2nd Vector In Plane
          T2[0] = x[2] - x[0];
          T2[1] = y[2] - y[0];
          T2[2] = z[2] - z[0];
          normalize( T2 );
          // Local Z-axis as cross product of x-axis and in-plane vector
          crossprod( T1, T2, T3 );
          normalize( T3 );
          // Local Y-axis as cross product of x-axis and z-axis
          crossprod( T3, T1, T2 );
          normalize( T2 );

          for(i=0; i<3; ++i)
            localg[i] = 0.0;

          for(i=0; i<3; ++i) {
            localg[0] += T1[i]*gravityAcceleration[i];
            localg[1] += T2[i]*gravityAcceleration[i];
            localg[2] += T3[i]*gravityAcceleration[i];
          }
          double localf[3];
          localf[0] = massPerNode*localg[0];
          localf[1] = massPerNode*localg[1];
          localf[2] = massPerNode*localg[2];

          fx = (T1[0]*localf[0]) + (T2[0]*localf[1]) + (T3[0]*localf[2]);
          fy = (T1[1]*localf[0]) + (T2[1]*localf[1]) + (T3[1]*localf[2]);
          fz = (T1[2]*localf[0]) + (T2[2]*localf[1]) + (T3[2]*localf[2]);
        }

        gravityForce[0] =  fx;
        gravityForce[1] =  fy;
        gravityForce[2] =  fz;
        gravityForce[3] =  fx;
        gravityForce[4] =  fy;
        gravityForce[5] =  fz;
        gravityForce[6] =  fx;
        gravityForce[7] =  fy;
        gravityForce[8] =  fz;
}

void
Membrane::getGravityForceSensitivityWRTthickness(CoordSet& cs, double *gravityAcceleration,
                                                 Vector& gravityForceSensitivity, int gravflg, GeomState *geomState)
{
        double mass = getMass(cs);
        double massPerNodePerThick = mass/(3.0*prop->eh);
        double fx, fy, fz;

        // Lumped
        if(gravflg != 2) {

          fx = massPerNodePerThick*gravityAcceleration[0];
          fy = massPerNodePerThick*gravityAcceleration[1];
          fz = massPerNodePerThick*gravityAcceleration[2];

        }
        // Consistent
        else {
          int i;
          Node &nd1 = cs.getNode(nn[0]);
          Node &nd2 = cs.getNode(nn[1]);
          Node &nd3 = cs.getNode(nn[2]);
          double x[3], y[3], z[3], localg[3];
          double T1[3],T2[3],T3[3];

          // Set the coordinates
          x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
          x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
          x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

          // Local X-axis
          T1[0] = x[1] - x[0];
          T1[1] = y[1] - y[0];
          T1[2] = z[1] - z[0];
          normalize( T1 );
          // 2nd Vector In Plane
          T2[0] = x[2] - x[0];
          T2[1] = y[2] - y[0];
          T2[2] = z[2] - z[0];
          normalize( T2 );
          // Local Z-axis as cross product of x-axis and in-plane vector
          crossprod( T1, T2, T3 );
          normalize( T3 );
          // Local Y-axis as cross product of x-axis and z-axis
          crossprod( T3, T1, T2 );
          normalize( T2 );

          for(i=0; i<3; ++i)
            localg[i] = 0.0;

          for(i=0; i<3; ++i) {
            localg[0] += T1[i]*gravityAcceleration[i];
            localg[1] += T2[i]*gravityAcceleration[i];
            localg[2] += T3[i]*gravityAcceleration[i];
          }
          double localf[3];
          localf[0] = massPerNodePerThick*localg[0];
          localf[1] = massPerNodePerThick*localg[1];
          localf[2] = massPerNodePerThick*localg[2];

          fx = (T1[0]*localf[0]) + (T2[0]*localf[1]) + (T3[0]*localf[2]);
          fy = (T1[1]*localf[0]) + (T2[1]*localf[1]) + (T3[1]*localf[2]);
          fz = (T1[2]*localf[0]) + (T2[2]*localf[1]) + (T3[2]*localf[2]);
        }

        gravityForceSensitivity[0] =  fx;
        gravityForceSensitivity[1] =  fy;
        gravityForceSensitivity[2] =  fz;
        gravityForceSensitivity[3] =  fx;
        gravityForceSensitivity[4] =  fy;
        gravityForceSensitivity[5] =  fz;
        gravityForceSensitivity[6] =  fx;
        gravityForceSensitivity[7] =  fy;
        gravityForceSensitivity[8] =  fz;
}

FullSquareMatrix
Membrane::massMatrix(CoordSet &cs,double *mel,int cmflg)
{
        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        double x21   = x[1] - x[0];
        double y21   = y[1] - y[0];
        double z21   = z[1] - z[0];
        double x32   = x[2] - x[1];
        double y32   = y[2] - y[1];
        double z32   = z[2] - z[1];
        double x13   = x[0] - x[2];
        double y13   = y[0] - y[2];
        double z13   = z[0] - z[2];

        double rl[3];
        rl[0] = sqrt(x21*x21 + y21*y21 + z21*z21);
        rl[1] = sqrt(x32*x32 + y32*y32 + z32*z32);
        rl[2] = sqrt(x13*x13 + y13*y13 + z13*z13);

        double rmas = getMass(cs)/3.0; 

        FullSquareMatrix ret(9,mel);

        ret.zero();

        ret[0][0] = rmas;
        ret[1][1] = rmas;
        ret[2][2] = (rl[0]*rl[0]+rl[2]*rl[2])/420.0*rmas;
        ret[3][3] = rmas;
        ret[4][4] = rmas;
        ret[5][5] = (rl[1]*rl[1]+rl[0]*rl[0])/420.0*rmas;
        ret[6][6] = rmas;
        ret[7][7] = rmas;
        ret[8][8] = (rl[2]*rl[2]+rl[1]*rl[1])/420.0*rmas;

        return ret;
}

FullSquareMatrix
Membrane::stiffness(CoordSet &cs, double *d, int flg)
{
	Node &nd1 = cs.getNode(nn[0]);
	Node &nd2 = cs.getNode(nn[1]);
	Node &nd3 = cs.getNode(nn[2]);

	double x[3], y[3], z[3], h[3];

        // Set the coordinates
	x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
	x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
	x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        // Set the thickness
	h[0] = h[1] = h[2] = prop->eh;

        if(h[0] <= 0.0)
          fprintf(stderr,"ERROR: Zero shell thickness (ThreeNodeShell.C) %d %d %d\n",
                nn[0], nn[1], nn[2]);

#ifdef USE_EIGEN3
        trimem(flg, x, y, z, prop->E, prop->nu, h, d);
#else
        _FORTRAN(trimem)(flg, x, y, z, prop->E, prop->nu, h, (double *)d);
#endif

        FullSquareMatrix ret(18,d);

        //std::cerr << "here in Membrane::stiffness\n"; ret.print();

        return ret;
}

void 
Membrane::getStiffnessThicknessSensitivity(CoordSet &cs, FullSquareMatrix &dStiffdThick, int flg, int senMethod)
{
#ifdef USE_EIGEN3
  if(dStiffdThick.dim() != 18) {
     std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
     exit(-1); 
  }

  // scalar parameters
  Eigen::Array<double,11,1> dconst;

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  Node &nd3 = cs.getNode(nn[2]);

  dconst[0] = nd1.x; dconst[1] = nd2.x; dconst[2] = nd3.x; // x coordinates
  dconst[3] = nd1.y; dconst[4] = nd2.y; dconst[5] = nd3.y; // y coordinates
  dconst[6] = nd1.z; dconst[7] = nd2.z; dconst[8] = nd3.z; // z coordinates
  dconst[9] = prop->E; // E
  dconst[10] = prop->nu;   // nu
  // integer parameters
  Eigen::Array<int,1,1> iconst;
  iconst[0] = flg; 
  // inputs
  Eigen::Matrix<double,1,1> q;
  q[0] = prop->eh;   // value of thickness at which jacobian is to be evaluated

  Eigen::Matrix<double,18,18> dStiffnessdThick;
  if(senMethod == 0) { // analytic
/*
    double x[3], y[3], z[3], h[3];
    x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
    x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
    x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
    h[0] = h[1] = h[2] = prop->eh;
   
    andesstfWRTthick(glNum+1, dStiffnessdThick.data(), prop->nu, x, y, z, type, flg);
    if(verboseFlag) std::cerr << "dStiffnessdThick(analytic) =\n" << dStiffnessdThick << std::endl;
*/
  }

  if(senMethod == 1) { // automatic differentiation
    Simo::FirstPartialSpaceDerivatives<double, MembraneStiffnessWRTThicknessSensitivity> dSdh(dconst,iconst); 
    Eigen::Array<Eigen::Matrix<double,18,18>,1,1> dStifdThick = dSdh(q, 0);
    if(verboseFlag) std::cerr << "dStifdThick(AD) =\n" << dStifdThick[0] << std::endl;
    dStiffnessdThick = dStifdThick[0];
  }

  if(senMethod == 2) { // finite difference
    MembraneStiffnessWRTThicknessSensitivity<double> foo(dconst,iconst);
    Eigen::Matrix<double,1,1> qp, qm;
    double h(1e-6);
    qp[0] = q[0] + h;   qm[0] = q[0] - h;
    Eigen::Matrix<double,18,18> Sp = foo(qp, 0);
    Eigen::Matrix<double,18,18> Sm = foo(qm, 0);
    dStiffnessdThick = (Sp-Sm)/(2*h);
    if(verboseFlag) std::cerr << "dStiffnessdThick(FD) =\n" << dStiffnessdThick << std::endl;
  }

  dStiffdThick.copy(dStiffnessdThick.data());
#endif
}

int
Membrane::numNodes()
{
 	return 3;
}

int*
Membrane::nodes(int *p)
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
Membrane::numDofs()
{
 	return 18;
}

int*
Membrane::dofs(DofSetArray &dsa, int *p)
{
 	if(p == 0) p = new int[18];
        dsa.number(nn[0], DofSet::XYZdisp | DofSet::XYZrot, p  );
        dsa.number(nn[1], DofSet::XYZdisp | DofSet::XYZrot, p+6);
        dsa.number(nn[2], DofSet::XYZdisp | DofSet::XYZrot, p+12);
	return p;
}

void
Membrane::markDofs(DofSetArray &dsa)
{
        dsa.mark(nn, 3, DofSet::XYZdisp | DofSet::XYZrot );
}

int
Membrane::getTopNumber()
{
  return 119;//4;
}

Corotator *
Membrane::getCorotator(CoordSet &cs, double *kel, int fitAlg, int)
{
 int flag = 0; // signals stiffness routine to keep local matrix 
 FullSquareMatrix myStiff = stiffness(cs, kel, flag);
 return new Shell3Corotator(nn[0], nn[1], nn[2], myStiff, fitAlg);
}

