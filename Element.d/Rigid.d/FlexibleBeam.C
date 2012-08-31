#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/FlexibleBeam.h>
#include <Element.d/Joint.d/LinearTranslationalSpring.h>
#include <Element.d/Joint.d/NonlinearTranslationalSpring2.h>
#include <Element.d/Joint.d/NonlinearTorsionalSpring.h>
#include <Corotational.d/utilities.h>

extern "C" {
 void _FORTRAN(e3dmas)(double&, double*,
                       double&, double*, double*, double*,
                       double*, double*, const int&, double&, const int&);
 void _FORTRAN(mass6)(double&, double&, double&, double&, double*,
                      double&, double*, double*, double*,
                      double*, double*, const int&, double&,
                      const int&, double*);
}

FlexibleBeam::FlexibleBeam(int* _nn)
 : SuperElement(true),
   elemframe(0)
{
  nnodes = 2;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 6; 
  subElems = new Element * [nSubElems];
  int indices[2] = { 0, 1 };

  subElems[0] = new LinearTranslationalSpring(indices);
  subElems[1] = new NonlinearTranslationalSpring2(indices, 1);
  subElems[2] = new NonlinearTranslationalSpring2(indices, 2);
  subElems[3] = new NonlinearTorsionalSpring(indices, 2, 1);
  subElems[4] = new NonlinearTorsionalSpring(indices, 1, 0);
  subElems[5] = new NonlinearTorsionalSpring(indices, 2, 0);
}

void
FlexibleBeam::setProp(StructProp *p, bool myProp)
{
  SuperElement::setProp(p, myProp);

  double G = p->E/(2*(1+p->nu));
  subElems[0]->getProperty()->penalty = p->E*p->A/l0;
  subElems[1]->getProperty()->penalty = 3*p->E*p->Iyy/(l0*l0*l0);
  subElems[2]->getProperty()->penalty = 3*p->E*p->Izz/(l0*l0*l0);
  subElems[3]->getProperty()->penalty = G*p->Ixx/l0;
  subElems[4]->getProperty()->penalty = p->E*p->Iyy/l0;
  subElems[5]->getProperty()->penalty = p->E*p->Izz/l0;
}

void
FlexibleBeam::buildFrame(CoordSet& cs)
{
  getLength(cs,l0);
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);
  if(elemframe != 0) {
    EFrame &theFrame = *elemframe;
    theFrame[0][0] = nd2.x-nd1.x;
    theFrame[0][1] = nd2.y-nd1.y;
    theFrame[0][2] = nd2.z-nd1.z;
    normalize(theFrame[0]);
    crossprod(theFrame[0],theFrame[1],theFrame[2]);
    normalize(theFrame[2]);
    crossprod(theFrame[2],theFrame[0],theFrame[1]);
  }
  else {
    c0[0][0] = nd2.x-nd1.x;
    c0[0][1] = nd2.y-nd1.y;
    c0[0][2] = nd2.z-nd1.z;
    normalize(c0[0]);
    double N1 = sqrt( c0[0][0]*c0[0][0] + c0[0][1]*c0[0][1] );
    double N2 = sqrt( c0[0][0]*c0[0][0] + c0[0][2]*c0[0][2] );

    if (N1 > N2) {
      c0[1][0] = -c0[0][1]/N1;
      c0[1][1] = c0[0][0]/N1;
      c0[1][2] = 0.0;
    }
    else {
      c0[1][0] = c0[0][2]/N2;
      c0[1][1] = 0.0;
      c0[1][2] = -c0[0][0]/N2;
    }

    c0[2][0] = c0[0][1] * c0[1][2] - c0[0][2] * c0[1][1];
    c0[2][1] = c0[0][2] * c0[1][0] - c0[0][0] * c0[1][2];
    c0[2][2] = c0[0][0] * c0[1][1] - c0[0][1] * c0[1][0];

    elemframe = &c0;
  }
  SuperElement::setFrame(elemframe);
  SuperElement::buildFrame(cs);
}

FullSquareMatrix
FlexibleBeam::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  if(prop == NULL || prop->rho == 0 || prop->A == 0) {
    FullSquareMatrix ret(numDofs(),mel);
    ret.zero();
    return ret;
  }

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double x[2], y[2], z[2];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

  double *gravityAcceleration = 0, *grvfor = 0, totmas = 0.0;

  int grvflg = 0, masflg = 0;

  double elmass[144];
  // Lumped Mass Matrix
  if (cmflg == 0) {
    _FORTRAN(e3dmas)(prop->rho,elmass,prop->A,
                     x,y,z,gravityAcceleration,grvfor,grvflg,totmas,masflg);

  // Consistent (Full) Mass Matrix
  } else {
    _FORTRAN(mass6)(prop->rho,prop->Ixx, prop->Iyy, prop->Izz,
                    elmass,prop->A,x,y,z,gravityAcceleration,
                    grvfor,grvflg,totmas,masflg,(double *) *elemframe);
  }

  FullSquareMatrix ret(numDofs(), mel);
  ret.zero();

  for(int i=0; i<12; ++i)
    for(int j=0; j<12; ++j) ret[i][j] = elmass[12*i+j];

  return ret;
}

double
FlexibleBeam::getMass(CoordSet& cs)
{
  if(prop == NULL || prop->rho == 0 || prop->A == 0)
    return 0;
  else 
    return l0*(prop->rho)*(prop->A);
}

void
FlexibleBeam::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                           Vector& gravityForce, int gravflg, GeomState *geomState)
{
        if(prop == NULL || prop->rho == 0 || prop->A == 0) {
          gravityForce.zero();
          return;
        }

        double massPerNode = 0.5*getMass(cs);

        double t0n[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
        double length;

        if (geomState) {
            updTransMatrix(cs, geomState, t0n, length);

        }  else  {
           getLength(cs, length);

           for(int i=0; i<3; ++i) {
              for(int j=0; j<3; ++j) {
                  t0n[i][j] = (*elemframe)[i][j] ;

              }
           }
        }

        int i;
        double localg[3];

        for(i=0; i<3; ++i)
          localg[i] = 0.0;

        for(i=0; i<3; ++i) {
          localg[0] += t0n[0][i]*gravityAcceleration[i];
          localg[1] += t0n[1][i]*gravityAcceleration[i];
          localg[2] += t0n[2][i]*gravityAcceleration[i];
        }
        double localf[3], localm[3];
        double globalf[3], globalm[3];
        localf[0] =  massPerNode*localg[0];
        localf[1] =  massPerNode*localg[1];
        localf[2] =  massPerNode*localg[2];
        if (gravflg == 2) { // consistent
          localm[0] =  0.0;
          localm[1] = -massPerNode*localg[2]*length/6.0;
          localm[2] =  massPerNode*localg[1]*length/6.0;
        }
        else if (gravflg == 1) { // lumped with fixed-end moments
          localm[0] =  0.0;
          localm[1] = -massPerNode*localg[2]*length/8.0;
          localm[2] =  massPerNode*localg[1]*length/8.0;
        }
        else localm[0] = localm[1] = localm[2] = 0.0; // lumped without fixed-end moments

        for(i=0; i<3; ++i) {
          globalf[i] = (t0n[0][i]*localf[0]) + (t0n[1][i]*localf[1]) + (t0n[2][i]*localf[2]);
          globalm[i] = (t0n[1][i]*localm[1]) + (t0n[2][i]*localm[2]);
        }

        gravityForce[0]  =  globalf[0];
        gravityForce[1]  =  globalf[1];
        gravityForce[2]  =  globalf[2];
        gravityForce[3]  =  globalm[0];
        gravityForce[4]  =  globalm[1];
        gravityForce[5]  =  globalm[2];
        gravityForce[6]  =  globalf[0];
        gravityForce[7]  =  globalf[1];
        gravityForce[8]  =  globalf[2];
        gravityForce[9]  = -globalm[0];
        gravityForce[10] = -globalm[1];
        gravityForce[11] = -globalm[2];

        for(int i=12; i<numDofs(); ++i) gravityForce[i] = 0;
}

void
FlexibleBeam::updTransMatrix(CoordSet& cs, GeomState *geomState, double t0n[3][3], double &length)
{
// Returns t0n[3][3] and length

   double  xn[2][3];

   double  zVecL[2][3];
   double (* rot[2])[3][3];

   // Get Nodes current coordinates
   NodeState &ns1 = (*geomState)[nn[0]];
   NodeState &ns2 = (*geomState)[nn[1]];

   xn[0][0]  = ns1.x; // x coordinate of node state 1
   xn[0][1]  = ns1.y; // y coordinate of node state 1
   xn[0][2]  = ns1.z; // z coordinate of node state 1

   xn[1][0]  = ns2.x; // x coordinate of node state 2
   xn[1][1]  = ns2.y; // y coordinate of node state 2
   xn[1][2]  = ns2.z; // z coordinate of node state 2

   double dx = xn[1][0] - xn[0][0];
   double dy = xn[1][1] - xn[0][1];
   double dz = xn[1][2] - xn[0][2];

   length = sqrt(dx*dx + dy*dy + dz*dz);

   rot[0]    = &(ns1.R); // rotation tensor of node state 1
   rot[1]    = &(ns2.R); // rotation tensor of node state 2

// Compute nodal rotated Z-axis in global coordinate system

   int i, nod;
   for(nod=0; nod<2; ++nod ) {
      for(i=0; i<3; ++i ) {
         zVecL[nod][i] = (*rot[nod])[i][0]*(*elemframe)[2][0]
                        +(*rot[nod])[i][1]*(*elemframe)[2][1]
                        +(*rot[nod])[i][2]*(*elemframe)[2][2];
      }
   }
/* Fitalg 1: Z-axis from node 1 */
   // We are setting fit Alg. to 2 to average z vectors.
   int fitAlg = 2;

   if (fitAlg == 1) {
      t0n[2][0] = zVecL[0][0];
      t0n[2][1] = zVecL[0][1];
      t0n[2][2] = zVecL[0][2];
   }

/* Fitalg .ne. 1: Z-axis as sum of nodal z-axis */
   else {
      t0n[2][0] = zVecL[0][0] + zVecL[1][0];
      t0n[2][1] = zVecL[0][1] + zVecL[1][1];
      t0n[2][2] = zVecL[0][2] + zVecL[1][2];
   }


   t0n[0][0]  = xn[1][0] - xn[0][0];
   t0n[0][1]  = xn[1][1] - xn[0][1];
   t0n[0][2]  = xn[1][2] - xn[0][2];

/* X-axis along element in Cn */
   normalize( t0n[0] );

/* Y-axis as cross product between z and x */
    crossprod( t0n[2], t0n[0], t0n[1]);
    normalize( t0n[1] );

/* Z-axis as cross product between x and y */
    crossprod( t0n[0], t0n[1], t0n[2]);
    normalize( t0n[2] );
}

void
FlexibleBeam::getLength(CoordSet& cs, double &length)
{
// Returns length of element

     Node &nd1 = cs.getNode(nn[0]);
     Node &nd2 = cs.getNode(nn[1]);

     double x[2], y[2], z[2];

     x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
     x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

     double dx = x[1] - x[0];
     double dy = y[1] - y[0];
     double dz = z[1] - z[0];

     length = sqrt(dx*dx + dy*dy + dz*dz);

}
#endif
