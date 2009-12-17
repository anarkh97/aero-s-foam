#include <Element.d/Beam.d/RigidMpcBeam.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>
#include <Corotational.d/utilities.h>
#define FOLLOWER_FORCE

RigidMpcBeam::RigidMpcBeam(int* _nn)
 : RigidMpcElement(2, 6, DofSet::XYZdisp | DofSet::XYZrot, 6, _nn)
{
  elemframe = 0;
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
    c0[0][0] = dx/length;
    c0[0][1] = dy/length;
    c0[0][2] = dz/length;

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

    // translation in x, y, z (1 constraint equation)
    for(int i = 0; i < 3; ++i)
      mpcs[0]->addterm(new LMPCTerm(nn[0], i, c0[0][i]));
    for(int i = 0; i < 3; ++i)
      mpcs[0]->addterm(new LMPCTerm(nn[1], i, -c0[0][i]));

    // rotation in x, y, z (3 constraint equations)
    for(int i = 0; i < 3; ++i) {
      for(int j = 0; j < 3; ++j) {
        double coef = (i == j) ? 1.0 : 0.0;
        mpcs[i+1]->addterm(new LMPCTerm(nn[0], 3+j, coef));
      }
      for(int j = 0; j < 3; ++j) {
        double coef = (i == j) ? -1.0 : 0.0;
        mpcs[i+1]->addterm(new LMPCTerm(nn[1], 3+j, coef));
      }
    }

    mpcs[4]->addterm(new LMPCTerm(nn[0], 0, c0[1][0]));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 1, c0[1][1]));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 2, c0[1][2]));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 3, -c0[1][1] * dz + c0[1][2] * dy));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 4, c0[1][0] * dz - c0[1][2] * dx));
    mpcs[4]->addterm(new LMPCTerm(nn[0], 5, -c0[1][0] * dy + c0[1][1] * dx));
    mpcs[4]->addterm(new LMPCTerm(nn[1], 0, -c0[1][0]));
    mpcs[4]->addterm(new LMPCTerm(nn[1], 1, -c0[1][1]));
    mpcs[4]->addterm(new LMPCTerm(nn[1], 2, -c0[1][2]));

    mpcs[5]->addterm(new LMPCTerm(nn[0], 0, c0[2][0]));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 1, c0[2][1]));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 2, c0[2][2]));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 3, -c0[2][1] * dz + c0[2][2] * dy));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 4, c0[2][0] * dz - c0[2][2] * dx));
    mpcs[5]->addterm(new LMPCTerm(nn[0], 5, -c0[2][0] * dy + c0[2][1] * dx));
    mpcs[5]->addterm(new LMPCTerm(nn[1], 0, -c0[2][0]));
    mpcs[5]->addterm(new LMPCTerm(nn[1], 1, -c0[2][1]));
    mpcs[5]->addterm(new LMPCTerm(nn[1], 2, -c0[2][2]));
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

#include "exp-map.h"
void 
RigidMpcBeam::updateLMPCs(GeomState& gState, CoordSet& cs)
{
  // nodes' original coordinates
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double dx0 = nd2.x - nd1.x;
  double dy0 = nd2.y - nd1.y;
  double dz0 = nd2.z - nd1.z;
  double l0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);

  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  double dx = ns2.x - ns1.x;
  double dy = ns2.y - ns1.y;
  double dz = ns2.z - ns1.z;
  double l = sqrt(dx*dx + dy*dy + dz*dz);

  if(l0 == 0.0) {
    /* linear */
  }
  else {
    double d[3] = { dx, dy, dz };

    // rotated cframes
    double c1[3][3], c2[3][3];
    mat_mult_mat(ns1.R, c0, c1, 1);
    mat_mult_mat(ns2.R, c0, c2, 1);

    // partial derivatives of constraint functions wrt x, y, z 
    for(int i = 0; i < 3; ++i) {
      mpcs[0]->terms[0+i].coef.r_value = -d[i]/l;
      mpcs[0]->terms[3+i].coef.r_value = d[i]/l;
      mpcs[4]->terms[0+i].coef.r_value = -c1[2][i];
      mpcs[4]->terms[6+i].coef.r_value = c1[2][i];
      mpcs[5]->terms[0+i].coef.r_value = -c1[1][i];
      mpcs[5]->terms[6+i].coef.r_value = c1[1][i];
    }

    // rotation parameters (thetax, thetay, thetaz)
    double r1[3], r2[3];
    mat_to_vec(ns1.R, r1);
    mat_to_vec(ns2.R, r2);
    //cerr << "theta1 = " << sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]) << ", theta2 = " << sqrt(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]) << endl;

    double dRdvi1[3][3], dRdvi2[3][3], d1[3][3], d2[3][3];
    for(int k=0; k<3; ++k) {
      // partial derivatives of rotation matrices wrt kth rotation parameters
      Partial_R_Partial_EM3(r1, k, dRdvi1);
      Partial_R_Partial_EM3(r2, k, dRdvi2);

      // partial derivatives of constraint functions wrt kth rotation parameters
      mat_mult_mat(dRdvi1, c0, d1, 1);
      mat_mult_mat(dRdvi2, c0, d2, 1);
      mpcs[1]->terms[0+k].coef.r_value = d1[2][0]*c2[1][0] + d1[2][1]*c2[1][1] + d1[2][2]*c2[1][2];
      mpcs[1]->terms[3+k].coef.r_value = c1[2][0]*d2[1][0] + c1[2][1]*d2[1][1] + c1[2][2]*d2[1][2];
      mpcs[2]->terms[0+k].coef.r_value = d1[2][0]*c2[0][0] + d1[2][1]*c2[0][1] + d1[2][2]*c2[0][2];
      mpcs[2]->terms[3+k].coef.r_value = c1[2][0]*d2[0][0] + c1[2][1]*d2[0][1] + c1[2][2]*d2[0][2];
      mpcs[3]->terms[0+k].coef.r_value = d1[1][0]*c2[0][0] + d1[1][1]*c2[0][1] + d1[1][2]*c2[0][2];
      mpcs[3]->terms[3+k].coef.r_value = c1[1][0]*d2[0][0] + c1[1][1]*d2[0][1] + c1[1][2]*d2[0][2];
      mpcs[4]->terms[3+k].coef.r_value = d1[2][0]*d[0] + d1[2][1]*d[1] + d1[2][2]*d[2];
      mpcs[5]->terms[3+k].coef.r_value = d1[1][0]*d[0] + d1[1][1]*d[1] + d1[1][2]*d[2];
    }

    // values of constraint functions
    mpcs[0]->rhs.r_value = l - l0;
    mpcs[1]->rhs.r_value = c1[2][0]*c2[1][0] + c1[2][1]*c2[1][1] + c1[2][2]*c2[1][2];
    mpcs[2]->rhs.r_value = c1[2][0]*c2[0][0] + c1[2][1]*c2[0][1] + c1[2][2]*c2[0][2];
    mpcs[3]->rhs.r_value = c1[1][0]*c2[0][0] + c1[1][1]*c2[0][1] + c1[1][2]*c2[0][2];
    mpcs[4]->rhs.r_value = c1[2][0]*d[0] + c1[2][1]*d[1] + c1[2][2]*d[2];
    mpcs[5]->rhs.r_value = c1[1][0]*d[0] + c1[1][1]*d[1] + c1[1][2]*d[2];
/*
    if(first) {
      for(int i=0; i<6; ++i) { first_rhs[i] = mpcs[i]->rhs.r_value; mpcs[i]->rhs.r_value = 0.0; }
    //  cerr << "first\n";
      first = false;
    } 
    else
      for(int i=0; i<6; ++i) mpcs[i]->rhs.r_value -= first_rhs[i];
    //cerr << "rhs = "; for(int i=0; i<6; ++i) cerr << mpcs[i]->rhs.r_value << " "; cerr << endl;
*/
  }

  //for(int i = 0; i<6; ++i) mpcs[i]->print();
}

void
RigidMpcBeam::getJacobian(GeomState& gState, CoordSet& cs, int k, FullSquareMatrix& J)
{
  J.zero();
  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  double lx = ns1.x-ns2.x;
  double ly = ns1.y-ns2.y;
  double lz = ns1.z-ns2.z;
  double l = sqrt(lx*lx + ly*ly + lz*lz);

  double dx = ns2.x - ns1.x;
  double dy = ns2.y - ns1.y;
  double dz = ns2.z - ns1.z;
  double d[3] = { dx, dy, dz };

  // rotated cframes
  double c1[3][3], c2[3][3];
  mat_mult_mat(ns1.R, c0, c1, 1);
  mat_mult_mat(ns2.R, c0, c2, 1);

  // rotation parameters (thetax, thetay, thetaz)
  double r1[3], r2[3];
  mat_to_vec(ns1.R, r1);
  mat_to_vec(ns2.R, r2);

  switch(k) { 

    case 0 : {

// [ 1/l - lx^2/l^3, -lx*ly/l^3,     -lx*lz/l^3,     0, 0, 0, lx^2/l^3 - 1/l, lx*ly/l^3,      lx*lz/l^3,      0, 0, 0]
// [-lx*ly/l^3,      1/l - ly^2/l^3, -ly*lz/l^3,     0, 0, 0, lx*ly/l^3,      ly^2/l^3 - 1/l, ly*lz/l^3,      0, 0, 0]
// [-lx*lz/l^3,      -ly*lz/l^3,     1/l - lz^2/l^3, 0, 0, 0, lx*lz/l^3,      ly*lz/l^3,      lz^2/l^3 - 1/l, 0, 0, 0]
// [ 0,              0,              0,              0, 0, 0, 0,              0,              0,              0, 0, 0]
// [ 0,              0,              0,              0, 0, 0, 0,              0,              0,              0, 0, 0]
// [ 0,              0,              0,              0, 0, 0, 0,              0,              0,              0, 0, 0]
// [ lx^2/l^3 - 1/l, lx*ly/l^3,      lx*lz/l^3,      0, 0, 0, 1/l - lx^2/l^3, -lx*ly/l^3,     -lx*lz/l^3,     0, 0, 0]
// [ lx*ly/l^3,      ly^2/l^3 - 1/l, ly*lz/l^3,      0, 0, 0, -lx*ly/l^3,     1/l - ly^2/l^3, -ly*lz/l^3,     0, 0, 0]
// [ lx*lz/l^3,      ly*lz/l^3,      lz^2/l^3 - 1/l, 0, 0, 0, -lx*lz/l^3,     -ly*lz/l^3,     1/l - lz^2/l^3, 0, 0, 0]
// [ 0,              0,              0,              0, 0, 0, 0,              0,              0,              0, 0, 0]
// [ 0,              0,              0,              0, 0, 0, 0,              0,              0,              0, 0, 0]
// [ 0,              0,              0,              0, 0, 0, 0,              0,              0,              0, 0, 0]

      double l2 = l*l;
      double l3 = l*l*l;
      J[0][0] = J[6][6] = (l2-lx*lx)/l3;
      J[1][1] = J[7][7] = (l2-ly*ly)/l3;
      J[2][2] = J[8][8] = (l2-lz*lz)/l3;

      J[0][1] = J[1][0] = J[6][7] = J[7][6] = -lx*ly/l3;
      J[0][2] = J[2][0] = J[6][8] = J[8][6] = -lx*lz/l3;
      J[1][2] = J[2][1] = J[7][8] = J[8][7] = -ly*lz/l3;

      J[0][6] = J[6][0] = -J[0][0];
      J[1][7] = J[7][1] = -J[1][1];
      J[2][8] = J[8][2] = -J[2][2];

      J[0][7] = J[7][0] = J[1][6] = J[6][1] = -J[0][1];
      J[0][8] = J[8][0] = J[2][6] = J[6][2] = -J[0][2];
      J[1][8] = J[8][1] = J[2][7] = J[7][2] = -J[1][2];

    } break;

    case 1 : {
      double d2Rdvidvj1[3][3], d2Rdvidvj2[3][3], dRdvi1[3][3], dRdvj2[3][3], d1[3][3], d2[3][3];
      for(int i=0; i<3; ++i) {
        // second partial derivatives of rotation matrices wrt rotation parameters
        for(int j=i; j<3; ++j) {
          Second_Partial_R_Partial_EM3(r1, i, j, d2Rdvidvj1);
          Second_Partial_R_Partial_EM3(r2, i, j, d2Rdvidvj2);

          mat_mult_mat(d2Rdvidvj1, c0, d1, 1);
          mat_mult_mat(d2Rdvidvj2, c0, d2, 1);

          J[3+i][3+j] = J[3+j][3+i] = d1[2][0]*c2[1][0] + d1[2][1]*c2[1][1] + d1[2][2]*c2[1][2];
          J[9+i][9+j] = J[9+j][9+i] = c1[2][0]*d2[1][0] + c1[2][1]*d2[1][1] + c1[2][2]*d2[1][2];
        }
      }
      for(int i=0; i<3; ++i) {
        Partial_R_Partial_EM3(r1, i, dRdvi1);
        mat_mult_mat(dRdvi1, c0, d1, 1);
        for(int j=0; j<3; ++j) {
          Partial_R_Partial_EM3(r2, j, dRdvj2);
          mat_mult_mat(dRdvj2, c0, d2, 1);

          J[3+i][9+j] = J[9+j][3+i] = d1[2][0]*d2[1][0] + d1[2][1]*d2[1][1] + d1[2][2]*d2[1][2];
        }
      }
    } break;

    case 2 : {
      double d2Rdvidvj1[3][3], d2Rdvidvj2[3][3], dRdvi1[3][3], dRdvj2[3][3], d1[3][3], d2[3][3];
      for(int i=0; i<3; ++i) {
        // second partial derivatives of rotation matrices wrt rotation parameters
        for(int j=i; j<3; ++j) {
          Second_Partial_R_Partial_EM3(r1, i, j, d2Rdvidvj1);
          Second_Partial_R_Partial_EM3(r2, i, j, d2Rdvidvj2);

          mat_mult_mat(d2Rdvidvj1, c0, d1, 1);
          mat_mult_mat(d2Rdvidvj2, c0, d2, 1);

          J[3+i][3+j] = J[3+j][3+i] = d1[2][0]*c2[0][0] + d1[2][1]*c2[0][1] + d1[2][2]*c2[0][2];
          J[9+i][9+j] = J[9+j][9+i] = c1[2][0]*d2[0][0] + c1[2][1]*d2[0][1] + c1[2][2]*d2[0][2];
        }
      }
      for(int i=0; i<3; ++i) {
        Partial_R_Partial_EM3(r1, i, dRdvi1);
        mat_mult_mat(dRdvi1, c0, d1, 1);
        for(int j=0; j<3; ++j) {
          Partial_R_Partial_EM3(r2, j, dRdvj2);
          mat_mult_mat(dRdvj2, c0, d2, 1);

          J[3+i][9+j] = J[9+j][3+i] = d1[2][0]*d2[0][0] + d1[2][1]*d2[0][1] + d1[2][2]*d2[0][2];
        }
      }
    } break;

    case 3 : {
      double d2Rdvidvj1[3][3], d2Rdvidvj2[3][3], dRdvi1[3][3], dRdvj2[3][3], d1[3][3], d2[3][3];
      for(int i=0; i<3; ++i) {
        // second partial derivatives of rotation matrices wrt rotation parameters
        for(int j=i; j<3; ++j) {
          Second_Partial_R_Partial_EM3(r1, i, j, d2Rdvidvj1);
          Second_Partial_R_Partial_EM3(r2, i, j, d2Rdvidvj2);

          mat_mult_mat(d2Rdvidvj1, c0, d1, 1);
          mat_mult_mat(d2Rdvidvj2, c0, d2, 1);

          J[3+i][3+j] = J[3+j][3+i] = d1[1][0]*c2[0][0] + d1[1][1]*c2[0][1] + d1[1][2]*c2[0][2];
          J[9+i][9+j] = J[9+j][9+i] = c1[1][0]*d2[0][0] + c1[1][1]*d2[0][1] + c1[1][2]*d2[0][2];
        }
      }
      for(int i=0; i<3; ++i) {
        Partial_R_Partial_EM3(r1, i, dRdvi1);
        mat_mult_mat(dRdvi1, c0, d1, 1);
        for(int j=0; j<3; ++j) {
          Partial_R_Partial_EM3(r2, j, dRdvj2);
          mat_mult_mat(dRdvj2, c0, d2, 1);

          J[3+i][9+j] = J[9+j][3+i] = d1[1][0]*d2[0][0] + d1[1][1]*d2[0][1] + d1[1][2]*d2[0][2];
        }
      }
    } break;

    case 4 : {
      double d2Rdvidvj1[3][3], dRdvi1[3][3], d1[3][3];
      for(int i=0; i<3; ++i) {
        // second partial derivatives of rotation matrices wrt rotation parameters
        for(int j=i; j<3; ++j) {
          Second_Partial_R_Partial_EM3(r1, i, j, d2Rdvidvj1);

          mat_mult_mat(d2Rdvidvj1, c0, d1, 1);

          J[3+i][3+j] = J[3+j][3+i] = d1[2][0]*d[0] + d1[2][1]*d[1] + d1[2][2]*d[2];
        }
      }
      for(int i=0; i<3; ++i) {
        Partial_R_Partial_EM3(r1, i, dRdvi1);
        mat_mult_mat(dRdvi1, c0, d1, 1);
        for(int j=0; j<3; ++j) {
          J[3+i][j] = J[j][3+i] = -d1[2][j];
          J[3+i][6+j] = J[6+j][3+i] = d1[2][j];
        }
      }
    } break;

    case 5 : {
      double d2Rdvidvj1[3][3], dRdvi1[3][3], d1[3][3];
      for(int i=0; i<3; ++i) {
        // second partial derivatives of rotation matrices wrt rotation parameters
        for(int j=i; j<3; ++j) {
          Second_Partial_R_Partial_EM3(r1, i, j, d2Rdvidvj1);

          mat_mult_mat(d2Rdvidvj1, c0, d1, 1);

          J[3+i][3+j] = J[3+j][3+i] = d1[1][0]*d[0] + d1[1][1]*d[1] + d1[1][2]*d[2];
        }
      }
      for(int i=0; i<3; ++i) {
        Partial_R_Partial_EM3(r1, i, dRdvi1);
        mat_mult_mat(dRdvi1, c0, d1, 1);
        for(int j=0; j<3; ++j) {
          J[3+i][j] = J[j][3+i] = -d1[1][j];
          J[3+i][6+j] = J[6+j][3+i] = d1[1][j];
        }
      }
    } break;

  }
  //cerr << "k = " << k << ", J = \n"; J.print();
}

