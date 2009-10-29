#include <Element.d/Beam.d/RigidMpcBeam.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>
#include <Corotational.d/utilities.h>
#define FOLLOWER_FORCE

RigidMpcBeam::RigidMpcBeam(int* _nn)
 : RigidMpcElement(2, 6, DofSet::XYZdisp | DofSet::XYZrot, 6, _nn)
{
  /* do nothing */
}

void 
RigidMpcBeam::computeMPCs(CoordSet& cs)
{
  for(int i = 0; i < nMpcs; ++i)
    mpcs[i] = new LMPCons(0, 0.0);

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double dx = nd2.x - nd1.x;
  double dy = nd2.y - nd1.y;
  double dz = nd2.z - nd1.z;

  double length = sqrt( dx*dx + dy*dy + dz*dz );

  if(length == 0.0) {
    cerr << " *** WARNING: Rigid beam has zero length, nodes: " << nn[0]+1 << " " << nn[1]+1 << endl;
    for(int i = 0; i < 6; ++i) {
      mpcs[i]->addterm(new LMPCTerm(nn[0], i, 1.0));
      mpcs[i]->addterm(new LMPCTerm(nn[1], i, -1.0));
    }
  }
  else {
    double c1[3], c2[3], c3[3];
    c1[0] = dx/length;
    c1[1] = dy/length;
    c1[2] = dz/length;

    // translation in x, y, z (1 constraint equation)
    for(int i = 0; i < 3; ++i) {
      mpcs[0]->addterm(new LMPCTerm(nn[0], i, c1[i]));
      mpcs[0]->addterm(new LMPCTerm(nn[1], i, -c1[i]));
    }

    // rotation in x, y, z (3 constraint equations)
    for(int i = 0; i < 3; ++i) {
      mpcs[i+1]->addterm(new LMPCTerm(nn[0], i+3, 1.0));
      mpcs[i+1]->addterm(new LMPCTerm(nn[1], i+3, -1.0));
    }

    // relative rotation in y and z (2 constraint equations)
    double N1 = sqrt( c1[0]*c1[0] + c1[1]*c1[1] );
    double N2 = sqrt( c1[0]*c1[0] + c1[2]*c1[2] );

    if (N1 > N2) {
      c2[0] = -c1[1]/N1;
      c2[1] = c1[0]/N1;
      c2[2] = 0.0;
    }
    else {
      c2[0] = c1[2]/N2;
      c2[1] = 0.0;
      c2[2] = -c1[0]/N2;
    }

    mpcs[4]->addterm(new LMPCTerm(nn[0], 0, c2[0]));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 1, c2[1]));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 2, c2[2]));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 3, -c2[1] * dz + c2[2] * dy));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 4, c2[0] * dz - c2[2] * dx));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 5, -c2[0] * dy + c2[1] * dx));
    mpcs[4]->addterm(new LMPCTerm(nn[1], 0, -c2[0]));
    mpcs[4]->addterm(new LMPCTerm(nn[1], 1, -c2[1]));
    mpcs[4]->addterm(new LMPCTerm(nn[1], 2, -c2[2]));

    c3[0] = c1[1] * c2[2] - c1[2] * c2[1];
    c3[1] = c1[2] * c2[0] - c1[0] * c2[2];
    c3[2] = c1[0] * c2[1] - c1[1] * c2[0];

    mpcs[5]->addterm(new LMPCTerm(nn[0], 0, c3[0]));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 1, c3[1]));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 2, c3[2]));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 3, -c3[1] * dz + c3[2] * dy));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 4, c3[0] * dz - c3[2] * dx));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 5, -c3[0] * dy + c3[1] * dx));
    mpcs[5]->addterm(new LMPCTerm(nn[1], 0, -c3[0]));
    mpcs[5]->addterm(new LMPCTerm(nn[1], 1, -c3[1]));
    mpcs[5]->addterm(new LMPCTerm(nn[1], 2, -c3[2]));
  }
}

int 
RigidMpcBeam::getTopNumber() 
{ 
  return 106; 
}

int
RigidMpcBeam::buildFrame(CoordSet& cs)
{
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
  return 0;
} 


void
RigidMpcBeam::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                   GeomState* geomState, int cflg)
{
  double t0n[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  double normal[3], normal2[3];
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double length;

  if (geomState) {
    updTransMatrix(cs, geomState, t0n, length);
#ifdef FOLLOWER_FORCE
    normal[0] = t0n[1][0];
    normal[1] = t0n[1][1];
    normal[2] = t0n[1][2];
    normal2[0] = t0n[2][0];
    normal2[1] = t0n[2][1];
    normal2[2] = t0n[2][2];
#else
    normal[0] = (*elemframe)[1][0];
    normal[1] = (*elemframe)[1][1];
    normal[2] = (*elemframe)[1][2];
    normal2[0] = (*elemframe)[2][0];
    normal2[1] = (*elemframe)[2][1];
    normal2[2] = (*elemframe)[2][2];
#endif
  }
  else {
    // Obtain normal to beam (second vector in element frame)
    normal[0] = (*elemframe)[1][0];
    normal[1] = (*elemframe)[1][1];
    normal[2] = (*elemframe)[1][2];
    normal2[0] = (*elemframe)[2][0];
    normal2[1] = (*elemframe)[2][1];
    normal2[2] = (*elemframe)[2][2];
    getLength(cs, length);
  }
  double pressureForce = 0.5*pressure*length;
  px = pressureForce*normal[0];
  py = pressureForce*normal[1];
  pz = pressureForce*normal[2];

  // Consistent
  double localMz = pressureForce*length/6.0;
  double mx = localMz*normal2[0];
  double my = localMz*normal2[1];
  double mz = localMz*normal2[2];
  elPressureForce[0]  = px;
  elPressureForce[1]  = py;
  elPressureForce[2]  = pz;
  elPressureForce[3]  = mx;
  elPressureForce[4]  = my;
  elPressureForce[5]  = mz;
  elPressureForce[6]  = px;
  elPressureForce[7]  = py;
  elPressureForce[8]  = pz;
  elPressureForce[9]  = -mx;
  elPressureForce[10] = -my;
  elPressureForce[11] = -mz;
  if(mode == 1) for(int i = 12; i < 18; ++i) elPressureForce[i] = 0.0;
}

void
RigidMpcBeam::updTransMatrix(CoordSet& cs, GeomState* geomState, double t0n[3][3], double& length)
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
RigidMpcBeam::getLength(CoordSet& cs, double& length)
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

