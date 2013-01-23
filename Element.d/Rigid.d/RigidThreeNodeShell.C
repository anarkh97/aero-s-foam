#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/RigidThreeNodeShell.h>
#include <Element.d/Rigid.d/RigidBeam.h>
#include <Corotational.d/utilities.h>

extern "C"      {
void _FORTRAN(mass8)(double* ,double* ,double* ,double* , double& ,
                     double* ,const int&, double* ,double*,const int&,
                     double&,const int&);
}

RigidThreeNodeShell::RigidThreeNodeShell(int *_nn)
 : SuperElement(true)
{
  nnodes = 3;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 2;
  subElems = new Element * [nSubElems];
  for(int i = 0; i < nSubElems; ++i) {
    int indices[2] = { i+1, 0 };
    subElems[i] = new RigidBeam(indices);
  }
}

FullSquareMatrix
RigidThreeNodeShell::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
        if(prop == NULL || prop->rho == 0 || prop->eh == 0) {
           FullSquareMatrix ret(numDofs(),mel);
           ret.zero();
           return ret;
        }

        Node &nd1 = cs.getNode(nn[0]);
        Node &nd2 = cs.getNode(nn[1]);
        Node &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];
        double *gravityAcceleration = NULL, *grvfor = NULL, totmas = 0.0;

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        h[0] = h[1] = h[2] = prop->eh;

        int grvflg = 0, masflg = 0;

        double *elmass = new double[324];
        _FORTRAN(mass8)(x,y,z,h,prop->rho,elmass,18,
                        gravityAcceleration,grvfor,grvflg,totmas,masflg);

        FullSquareMatrix ret(numDofs(),mel);
        ret.zero();

        for(int i=0; i<18; ++i)
          for(int j=0; j<18; ++j) ret[i][j] = elmass[18*i+j];
        delete [] elmass;

        return ret;
}

double
RigidThreeNodeShell::getMass(CoordSet& cs)
{
        if (prop == NULL || prop->rho == 0 || prop->eh == 0) return 0.0;

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

        crossprod(v1, v2, v3);

        double area = 0.5*sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
        double mass = area*prop->rho*prop->eh;

        return mass;

}

void
RigidThreeNodeShell::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                                     Vector& gravityForce, int gravflg, GeomState *geomState)
{
        if(prop == NULL || prop->rho == 0 || prop->eh == 0) {
          gravityForce.zero(); 
          return;
        }

        double mass = getMass(cs);
        double massPerNode = mass/3.0;

        double fx = massPerNode*gravityAcceleration[0];
        double fy = massPerNode*gravityAcceleration[1];
        double fz = massPerNode*gravityAcceleration[2];
        double mx[3],my[3],mz[3];
        int i;
        for(i = 0; i < 3; ++i) {
          mx[i] = 0.0;
          my[i] = 0.0;
          mz[i] = 0.0;
        }

        // Lumped with no fixed-end moments
        if (gravflg == 0) {

        }
        // Consistent or lumped with fixed end moments.  Compute treating shell as 3 beams.
        else {

          Node &nd1 = cs.getNode(nn[0]);
          Node &nd2 = cs.getNode(nn[1]);
          Node &nd3 = cs.getNode(nn[2]);
          double x[3], y[3], z[3];
          x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
          x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
          x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

          double T1[3],T2[3],T3[3];
          // Vector 1 from Node 1->2
          T1[0] = x[1] - x[0];
          T1[1] = y[1] - y[0];
          T1[2] = z[1] - z[0];
          normalize( T1 );
          // Vector 2 from Node 1->3
          T2[0] = x[2] - x[0];
          T2[1] = y[2] - y[0];
          T2[2] = z[2] - z[0];
          normalize( T2 );
          // Local Z-axis as cross between V1 and V2
          crossprod( T1, T2, T3 );
          normalize( T3);

          int beam, beamnode[3][2];
          beamnode[0][0] = 0;
          beamnode[0][1] = 1;
          beamnode[1][0] = 0;
          beamnode[1][1] = 2;
          beamnode[2][0] = 1;
          beamnode[2][1] = 2;

          for(beam=0; beam<3; ++beam) {
            double length, dx, dy, dz, localg[3];
            int n1, n2;
            n1 = beamnode[beam][0];
            n2 = beamnode[beam][1];
            dx = x[n2] - x[n1];
            dy = y[n2] - y[n1];
            dz = z[n2] - z[n1];
            length = sqrt(dx*dx + dy*dy + dz*dz);
            // Local X-axis from Node 1->2
            T1[0] = x[n2] - x[n1];
            T1[1] = y[n2] - y[n1];
            T1[2] = z[n2] - z[n1];
            normalize( T1 );
            // Local Y-axis as cross between Z and X
            crossprod( T3, T1, T2 );
            normalize( T2);

            for(i = 0; i < 3; ++i)
              localg[i] = 0.0;
            for(i = 0; i < 3; ++i) {
              localg[0] += T1[i]*gravityAcceleration[i];
              localg[1] += T2[i]*gravityAcceleration[i];
              localg[2] += T3[i]*gravityAcceleration[i];
            }
            double lmy,lmz;
            if (gravflg == 2) { // consistent
              lmy = -massPerNode*localg[2]*length/12.0;
              lmz = massPerNode*localg[1]*length/12.0;
            }
            else { // lumped with fixed-end moments
              lmy = -massPerNode*localg[2]*length/16.0;
              lmz = massPerNode*localg[1]*length/16.0;
            }
            mx[n1] += ((T2[0]*lmy) + (T3[0]*lmz));
            my[n1] += ((T2[1]*lmy) + (T3[1]*lmz));
            mz[n1] += ((T2[2]*lmy) + (T3[2]*lmz));
            mx[n2] -= ((T2[0]*lmy) + (T3[0]*lmz));
            my[n2] -= ((T2[1]*lmy) + (T3[1]*lmz));
            mz[n2] -= ((T2[2]*lmy) + (T3[2]*lmz));
          }
        }

        gravityForce[0]  = fx;
        gravityForce[1]  = fy;
        gravityForce[2]  = fz;
        gravityForce[3]  = mx[0];
        gravityForce[4]  = my[0];
        gravityForce[5]  = mz[0];
        gravityForce[6]  = fx;
        gravityForce[7]  = fy;
        gravityForce[8]  = fz;
        gravityForce[9]  = mx[1];
        gravityForce[10] = my[1];
        gravityForce[11] = mz[1];
        gravityForce[12] = fx;
        gravityForce[13] = fy;
        gravityForce[14] = fz;
        gravityForce[15] = mx[2];
        gravityForce[16] = my[2];
        gravityForce[17] = mz[2];

        for(int i=18; i<numDofs(); ++i) gravityForce[i] = 0;
}
#endif
