#include <Element.d/Beam.d/RigidMpcBeam.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>
#include <Corotational.d/utilities.h>
#define FOLLOWER_FORCE

#define PJSA_DEBUG

RigidMpcBeam::RigidMpcBeam(int* _nn)
 : RigidMpcElement(2, 6, DofSet::XYZdisp | DofSet::XYZrot, 6, _nn)
{
  elemframe = 0;
}

void 
RigidMpcBeam::computeMPCs(CoordSet& cs)
{
#ifdef PJSA_DEBUG
  for(int i = 0; i < nMpcs; ++i)
    mpcs[i] = new LMPCons(0, 0.0);

  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double lx0 = nd2.x - nd1.x;
  double ly0 = nd2.y - nd1.y;
  double lz0 = nd2.z - nd1.z;

  for(int i = 0; i < 3; ++i)
    mpcs[i]->addterm(new LMPCTerm(nn[0], i, -1.0));
  mpcs[0]->addterm(new LMPCTerm(nn[0], 3, 0.0));
  mpcs[0]->addterm(new LMPCTerm(nn[0], 4, -lz0));
  mpcs[0]->addterm(new LMPCTerm(nn[0], 5, ly0));
  mpcs[1]->addterm(new LMPCTerm(nn[0], 3, lz0));
  mpcs[1]->addterm(new LMPCTerm(nn[0], 4, 0.0));
  mpcs[1]->addterm(new LMPCTerm(nn[0], 5, -lx0));
  mpcs[2]->addterm(new LMPCTerm(nn[0], 3, -ly0));
  mpcs[2]->addterm(new LMPCTerm(nn[0], 4, lx0)); 
  mpcs[2]->addterm(new LMPCTerm(nn[0], 5, 0.0));
  for(int i = 0; i < 3; ++i)
      mpcs[i]->addterm(new LMPCTerm(nn[1], i, 1.0));

  for(int i = 0; i < 3; ++i) {
    mpcs[i+3]->addterm(new LMPCTerm(nn[0], i+3, 1.0));
    mpcs[i+3]->addterm(new LMPCTerm(nn[1], i+3, -1.0));
  }
  //for(int i = 0; i<6; ++i) {
  //  cerr << "lmpc " << i+1 << " has coefs "; for(int j=0; j<mpcs[i]->nterms; ++j) cerr << mpcs[i]->terms[j].coef.r_value << " "; cerr << endl;
  //}
#else
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
    for(int i = 0; i < 3; ++i)
      mpcs[0]->addterm(new LMPCTerm(nn[0], i, c1[i]));
    for(int i = 0; i < 3; ++i)
      mpcs[0]->addterm(new LMPCTerm(nn[1], i, -c1[i]));

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
#endif
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
#ifdef PJSA_DEBUG
  // nodes' original coordinates
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double lx0 = nd2.x - nd1.x;
  double ly0 = nd2.y - nd1.y;
  double lz0 = nd2.z - nd1.z;
  double l0 = sqrt(lx0*lx0+ly0*ly0+lz0*lz0);

  // nodes' current state
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];
  double lx = ns2.x - ns1.x;
  double ly = ns2.y - ns1.y;
  double lz = ns2.z - ns1.z;
  double l = sqrt(lx*lx+ly*ly+lz*lz);
  double r1[3], r2[3];
  mat_to_vec(ns1.R, r1);
  mat_to_vec(ns2.R, r2);
  double thetax1 = r1[0], thetay1 = r1[1], thetaz1 = r1[2];
  double thetax12 = thetax1*thetax1, thetay12 = thetay1*thetay1, thetaz12 = thetaz1*thetaz1;

/*
  mpcs[0]->terms[1].coef.r_value = - ly0*(thetay1/2 + (thetax1*thetaz1)/3) - lz0*(thetaz1/2 - (thetax1*thetay1)/3);
  mpcs[0]->terms[2].coef.r_value = lx0*thetay1 - ly0*(thetax1/2 + (thetay1*thetaz1)/3) + lz0*(thetax12/6 + thetay12/2 + thetaz12/6 - 1); //-lz0*exp(r1[1]);
  mpcs[0]->terms[3].coef.r_value = lx0*thetaz1 - lz0*(thetax1/2 - (thetay1*thetaz1)/3) - ly0*(thetax12/6 + thetay12/6 + thetaz12/2 - 1); //ly0*exp(r1[2]);
  mpcs[1]->terms[1].coef.r_value = ly0*thetax1 - lx0*(thetay1/2 - (thetax1*thetaz1)/3) - lz0*(thetax12/2 + thetay12/6 + thetaz12/6 - 1); //lz0*exp(r1[0]);
  mpcs[1]->terms[2].coef.r_value = - lx0*(thetax1/2 - (thetay1*thetaz1)/3) - lz0*(thetaz1/2 + (thetax1*thetay1)/3);
  mpcs[1]->terms[3].coef.r_value = ly0*thetaz1 - lz0*(thetay1/2 + (thetax1*thetaz1)/3) + lx0*(thetax12/6 + thetay12/6 + thetaz12/2 - 1); //-lx0*exp(r1[2]);
  mpcs[2]->terms[1].coef.r_value = lz0*thetax1 - lx0*(thetaz1/2 + (thetax1*thetay1)/3) + ly0*(thetax12/2 + thetay12/6 + thetaz12/6 - 1); //-ly0*exp(r1[0]);
  mpcs[2]->terms[2].coef.r_value = lz0*thetay1 - ly0*(thetaz1/2 - (thetax1*thetay1)/3) - lx0*(thetax12/6 + thetay12/2 + thetaz12/6 - 1); //lx0*exp(r1[1]);
  mpcs[2]->terms[3].coef.r_value = - lx0*(thetax1/2 + (thetay1*thetaz1)/3) - ly0*(thetay1/2 - (thetax1*thetaz1)/3);
*/

  mpcs[0]->rhs.r_value = ns2.x - ns1.x - ns1.R[0][0]*lx0 - ns1.R[0][1]*ly0 - ns1.R[0][2]*lz0;
  mpcs[1]->rhs.r_value = ns2.y - ns1.y - ns1.R[1][0]*lx0 - ns1.R[1][1]*ly0 - ns1.R[1][2]*lz0;
  mpcs[2]->rhs.r_value = ns2.z - ns1.z - ns1.R[2][0]*lx0 - ns1.R[2][1]*ly0 - ns1.R[2][2]*lz0;
  mpcs[3]->rhs.r_value = (r1[0] - r2[0]);
  mpcs[4]->rhs.r_value = (r1[1] - r2[1]);
  mpcs[5]->rhs.r_value = (r1[2] - r2[2]);

/*
  cerr << "R1 = " << ns1.R[0][0] << " " << ns1.R[0][1] << " " << ns1.R[0][2] << endl;
  cerr << "     " << ns1.R[1][0] << " " << ns1.R[1][1] << " " << ns1.R[1][2] << endl;
  cerr << "     " << ns1.R[2][0] << " " << ns1.R[2][1] << " " << ns1.R[2][2] << endl;
  cerr << "R2 = " << ns2.R[0][0] << " " << ns2.R[0][1] << " " << ns2.R[0][2] << endl;
  cerr << "     " << ns2.R[1][0] << " " << ns2.R[1][1] << " " << ns2.R[1][2] << endl;
  cerr << "     " << ns2.R[2][0] << " " << ns2.R[2][1] << " " << ns2.R[2][2] << endl;
  cerr << "r1 = " << r1[0] << " " << r1[1] << " " << r1[2] << endl;
  cerr << "r2 = " << r2[0] << " " << r2[1] << " " << r2[2] << endl;
  for(int i=0; i<6; ++i) cerr << "mpc " << i << " rhs = " << mpcs[i]->rhs.r_value << endl;
  cerr << "l - l0 = " << l - l0 << endl;
*/

/*
  double R[4][4];
  EM3_To_R(r1, R);
  cerr << "R  = " << R[0][0] << " " << R[0][1] << " " << R[0][2] << " " << R[0][3] << endl;
  cerr << "     " << R[1][0] << " " << R[1][1] << " " << R[1][2] << " " << R[1][3] << endl;
  cerr << "     " << R[2][0] << " " << R[2][1] << " " << R[2][2] << " " << R[1][3] << endl;
  cerr << "     " << R[3][0] << " " << R[3][1] << " " << R[3][2] << " " << R[3][3] << endl;
*/
  int rep;
  double dRdvi[4][4];
  for(int i=0; i<3; ++i) {
    rep = Partial_R_Partial_EM3(r1, i, dRdvi);
    //if(rep) cerr << "dynamic reparametrization activated\n";
    //cerr << "i = " << i << ", rep = " << rep << endl;
    //cerr << "dRdvi  = " << dRdvi[0][0] << " " << dRdvi[0][1] << " " << dRdvi[0][2] << " " << dRdvi[0][3] << endl;
    //cerr << "         " << dRdvi[1][0] << " " << dRdvi[1][1] << " " << dRdvi[1][2] << " " << dRdvi[1][3] << endl;
    //cerr << "         " << dRdvi[2][0] << " " << dRdvi[2][1] << " " << dRdvi[2][2] << " " << dRdvi[1][3] << endl;
    //cerr << "         " << dRdvi[3][0] << " " << dRdvi[3][1] << " " << dRdvi[3][2] << " " << dRdvi[3][3] << endl;
    for(int j=0; j<3; ++j)
      mpcs[j]->terms[i+1].coef.r_value = -(dRdvi[j][0]*lx0+dRdvi[j][1]*ly0+dRdvi[j][2]*lz0);
  }
  //for(int i = 0; i<6; ++i) {
  //  cerr << "lmpc " << i+1 << " has coefs "; for(int j=0; j<mpcs[i]->nterms; ++j) cerr << mpcs[i]->terms[j].coef.r_value << " "; cerr << endl;
  //}
/*
  mpcs[0]->rhs.r_value = ns2.x - ns1.x - ns1.R[0][0]*lx0 - ns1.R[0][1]*ly0 - ns1.R[0][2]*lz0;
  mpcs[1]->rhs.r_value = ns2.y - ns1.y - ns1.R[1][0]*lx0 - ns1.R[1][1]*ly0 - ns1.R[1][2]*lz0;
  mpcs[2]->rhs.r_value = ns2.z - ns1.z - ns1.R[2][0]*lx0 - ns1.R[2][1]*ly0 - ns1.R[2][2]*lz0;
  mpcs[3]->rhs.r_value = (r1[0] - r2[0]);
  mpcs[4]->rhs.r_value = (r1[1] - r2[1]);
  mpcs[5]->rhs.r_value = (r1[2] - r2[2]); 
*/
#else
  // nodes' original coordinates
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double dx0 = nd2.x - nd1.x;
  double dy0 = nd2.y - nd1.y;
  double dz0 = nd2.z - nd1.z;
  double l0 = sqrt(dx*dx + dy*dy + dz*dz);

  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  double dx = ns2.x - ns1.x;
  double dy = ns2.y - ns1.y;
  double dz = ns2.z - ns1.z;
  double l = sqrt(lx*lx + ly*ly + lz*lz);

  if(l0 == 0.0) {
    /* revisit this case */
  }
  else {
    double c1[3], c2[3], c3[3];
    c1[0] = dx/l;
    c1[1] = dy/l;
    c1[2] = dz/l;

    // translation in x, y, z (1 constraint equation)
    for(int i = 0; i < 3; ++i) {
      mpcs[0]->terms[i].coef.r_value = c1[i];
      mpcs[0]->terms[3+i].coef.r_value = -c1[i] ;
    }
    mpcs[0]->rhs.r_value = l - l0;

    // rotation in x, y, z (3 constraint equations)
    mpcs[1]->rhs.r_value = thetax1-thetax2;
    mpcs[2]->rhs.r_value = thetay1-thetay2;
    mpcs[3]->rhs.r_value = thetaz1-thetaz2;

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

    mpcs[4]->terms[0].coef.r_value = c2[0];
    mpcs[4]->terms[1].coef.r_value = c2[1];
    mpcs[4]->terms[2].coef.r_value = c2[2];
    mpcs[4]->terms[3].coef.r_value = -c2[1] * dz + c2[2] * dy;
    mpcs[4]->terms[4].coef.r_value = c2[0] * dz - c2[2] * dx;
    mpcs[4]->terms[5].coef.r_value = -c2[0] * dy + c2[1] * dx;
    mpcs[4]->terms[6].coef.r_value = -c2[0];
    mpcs[4]->terms[7].coef.r_value = -c2[1];
    mpcs[4]->terms[8].coef.r_value = -c2[2];
    mpcs[4]->rhs.r_value = xxxx;

    c3[0] = c1[1] * c2[2] - c1[2] * c2[1];
    c3[1] = c1[2] * c2[0] - c1[0] * c2[2];
    c3[2] = c1[0] * c2[1] - c1[1] * c2[0];

    mpcs[5]->terms[0].coef.r_value = c3[0];
    mpcs[5]->terms[1].coef.r_value = c3[1];
    mpcs[5]->terms[2].coef.r_value = c3[2];
    mpcs[5]->terms[3].coef.r_value = -c3[1] * dz + c3[2] * dy;
    mpcs[5]->terms[4].coef.r_value = c3[0] * dz - c3[2] * dx;
    mpcs[5]->terms[5].coef.r_value = -c3[0] * dy + c3[1] * dx;
    mpcs[5]->terms[6].coef.r_value = -c3[0];
    mpcs[5]->terms[7].coef.r_value = -c3[1];
    mpcs[5]->terms[8].coef.r_value = -c3[2];
    mpcs[5]->rhs.r_value = yyyy;
  }
#endif
}

void
RigidMpcBeam::getJacobian(GeomState& gState, CoordSet& cs, int i, FullSquareMatrix& J)
{
#ifdef PJSA_DEBUG
  // nodes' original coordinates
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double lx0 = nd2.x - nd1.x;
  double ly0 = nd2.y - nd1.y;
  double lz0 = nd2.z - nd1.z;
    
  // nodes' current state
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];
  double lx = ns2.x - ns1.x;
  double ly = ns2.y - ns1.y;
  double lz = ns2.z - ns1.z;
  double r1[3];
  mat_to_vec(ns1.R, r1);
  double thetax1 = r1[0], thetay1 = r1[1], thetaz1 = r1[2];
  double thetax12 = thetax1*thetax1, thetay12 = thetay1*thetay1, thetaz12 = thetaz1*thetaz1;
  double thetax13 = thetax12*thetax1, thetay13 = thetay12*thetay1, thetaz13 = thetaz12*thetaz1;

  J.zero();
/*
  switch(i) {
    case 0 : {
      J[3][3] = (lz0*thetay1)/3 - (ly0*thetaz1)/3;
      J[4][4] = lx0 - (ly0*thetaz1)/3 + lz0*thetay1;
      J[5][5] = lx0 - ly0*thetaz1 + (lz0*thetay1)/3;
      J[3][4] = J[4][3] = (lz0*thetax1)/3 - ly0/2;
      J[3][5] = J[5][3] = - lz0/2 - (ly0*thetax1)/3;
      J[4][5] = J[5][4] = (lz0*thetaz1)/3 - (ly0*thetay1)/3;
    } break;
    case 1 : {
      J[3][3] = ly0 + (lx0*thetaz1)/3 - lz0*thetax1;
      J[4][4] = (lx0*thetaz1)/3 - (lz0*thetax1)/3;
      J[5][5] = ly0 + lx0*thetaz1 - (lz0*thetax1)/3;
      J[3][4] = J[4][3] = - lx0/2 - (lz0*thetay1)/3;
      J[3][5] = J[5][3] = (lx0*thetax1)/3 - (lz0*thetaz1)/3;
      J[4][5] = J[5][4] = (lx0*thetay1)/3 - lz0/2;
    } break;
    case 2 : {
      J[3][3] = lz0 - (lx0*thetay1)/3 + ly0*thetax1;
      J[4][4] = lz0 - lx0*thetay1 + (ly0*thetax1)/3;
      J[5][5] = (ly0*thetax1)/3 - (lx0*thetay1)/3;
      J[3][4] = J[4][3] = (ly0*thetay1)/3 - (lx0*thetax1)/3;
      J[3][5] = J[5][3] = (ly0*thetaz1)/3 - lx0/2;
      J[4][5] = J[5][4] = - ly0/2 - (lx0*thetaz1)/3;
    } break;
    default : {
    } break;
  }
*/

/*
  switch(i) {
    case 0 : {
      J[3][3] = lz0*(thetay1/3.0 + (thetax1*thetaz1)/4.0) - ly0*(thetaz1/3.0 - (thetax1*thetay1)/4.0) - lx0*(thetay12/12.0 + thetaz12/12.0);
      J[4][4] = lz0*(thetay1 + (thetax1*thetaz1)/12.0) - lx0*(thetax12/12.0 + thetay12/2 + thetaz12/6.0 - 1.0) - ly0*(thetaz1/3.0 - (thetax1*thetay1)/4.0);
      J[5][5] = lz0*(thetay1/3.0 + (thetax1*thetaz1)/4.0) - lx0*(thetax12/12.0 + thetay12/6.0 + thetaz12/2 - 1.0) - ly0*(thetaz1 - (thetax1*thetay1)/12.0);
      J[4][3] = J[3][4] = lz0*(thetax1/3.0 + (thetay1*thetaz1)/12.0) + ly0*(thetax12/8.0 + thetay12/8.0 + thetaz12/24 - 0.5) - (lx0*thetax1*thetay1)/6.0;
      J[5][3] = J[3][5] = lz0*(thetax12/8.0 + thetay12/24 + thetaz12/8.0 - 0.5) - ly0*(thetax1/3.0 - (thetay1*thetaz1)/12.0) - (lx0*thetax1*thetaz1)/6.0;
      J[5][4] = J[4][5] = lz0*(thetaz1/3.0 + (thetax1*thetay1)/12.0) - ly0*(thetay1/3.0 - (thetax1*thetaz1)/12.0) - (lx0*thetay1*thetaz1)/3.0;
    } break;
    case 1 : {
      J[3][3] = lx0*(thetaz1/3.0 + (thetax1*thetay1)/4.0) - ly0*(thetax12/2 + thetay12/12.0 + thetaz12/6.0 - 1.0) - lz0*(thetax1 - (thetay1*thetaz1)/12.0);
      J[4][4] = lx0*(thetaz1/3.0 + (thetax1*thetay1)/4.0) - lz0*(thetax1/3.0 - (thetay1*thetaz1)/4.0) - ly0*(thetax12/12.0 + thetaz12/12.0);
      J[5][5] = lx0*(thetaz1 + (thetax1*thetay1)/12.0) - ly0*(thetax12/6.0 + thetay12/12.0 + thetaz12/2 - 1.0) - lz0*(thetax1/3.0 - (thetay1*thetaz1)/4.0);
      J[4][3] = J[3][4] = lx0*(thetax12/8.0 + thetay12/8.0 + thetaz12/24 - 0.5) - lz0*(thetay1/3.0 - (thetax1*thetaz1)/12.0) - (ly0*thetax1*thetay1)/6.0;
      J[5][3] = J[3][5] = lx0*(thetax1/3.0 + (thetay1*thetaz1)/12.0) - lz0*(thetaz1/3.0 - (thetax1*thetay1)/12.0) - (ly0*thetax1*thetaz1)/3.0;
      J[5][4] = J[4][5] = lx0*(thetay1/3.0 + (thetax1*thetaz1)/12.0) + lz0*(thetax12/24 + thetay12/8.0 + thetaz12/8.0 - 0.5) - (ly0*thetay1*thetaz1)/6.0;
    } break;
    case 2 : {
      J[3][3] = ly0*(thetax1 + (thetay1*thetaz1)/12.0) - lz0*(thetax12/2 + thetay12/6.0 + thetaz12/12.0 - 1.0) - lx0*(thetay1/3.0 - (thetax1*thetaz1)/4.0);
      J[4][4] = ly0*(thetax1/3.0 + (thetay1*thetaz1)/4.0) - lz0*(thetax12/6.0 + thetay12/2 + thetaz12/12.0 - 1.0) - lx0*(thetay1 - (thetax1*thetaz1)/12.0);
      J[5][5] = ly0*(thetax1/3.0 + (thetay1*thetaz1)/4.0) - lx0*(thetay1/3.0 - (thetax1*thetaz1)/4.0) - lz0*(thetax12/12.0 + thetay12/12.0);
      J[4][3] = J[3][4] = ly0*(thetay1/3.0 + (thetax1*thetaz1)/12.0) - lx0*(thetax1/3.0 - (thetay1*thetaz1)/12.0) - (lz0*thetax1*thetay1)/3.0;
      J[5][3] = J[3][5] = ly0*(thetaz1/3.0 + (thetax1*thetay1)/12.0) + lx0*(thetax12/8.0 + thetay12/24 + thetaz12/8.0 - 0.5) - (lz0*thetax1*thetaz1)/6.0;
      J[5][4] = J[4][5] = ly0*(thetax12/24 + thetay12/8.0 + thetaz12/8.0 - 0.5) - lx0*(thetaz1/3.0 - (thetax1*thetay1)/12.0) - (lz0*thetay1*thetaz1)/6.0;
    } break;
    default : {
    } break;
  }
*/


switch(i) {
case 0 : {
double H[3][3] = { 
{ ly0*((thetax1*thetay1)/4 - thetaz1/3 + (thetax12*thetaz1)/60 + (thetaz1*(12*thetax12 + 2*thetax12 + 4*thetax12))/120) - lx0*(thetax12/12 + thetax12/12) + lz0*(thetay1/3 + (thetax1*thetaz1)/4 - (thetay1*thetax12)/60 - (thetay1*(12*thetax12 + 4*thetax12 + 2*thetax12))/120), 
  ly0*(thetax12/8 + (thetax1*thetay1*thetaz1)/15 + thetax12/8 + thetax12/24 - 0.5) - lz0*((thetax1*(3*thetax12 + thetax12 + thetax12))/120 - thetax1/3 - (thetay1*thetaz1)/12 + (thetax1*thetax12)/12 + (thetax1*thetax12)/40 + (thetax1*(thetax12 + thetax12))/120) + lx0*((thetax12*thetaz1)/120 - (thetax1*thetay1)/6 - (thetaz1*(thetax12 + 3*thetax12 + thetax12))/120 + (thetax12*thetaz1)/60 + (thetaz1*(thetax12 + thetax12))/120), 
 ly0*((thetax1*(3*thetax12 + thetax12 + thetax12))/120 - thetax1/3 + (thetay1*thetaz1)/12 + (thetax1*thetax12)/40 + (thetax1*thetax12)/12 + (thetax1*(thetax12 + thetax12))/120) + lz0*(thetax12/8 - (thetax1*thetay1*thetaz1)/15 + thetax12/24 + thetax12/8 - 0.5) - lx0*((thetax1*thetaz1)/6 - (thetay1*(thetax12 + thetax12 + 3*thetax12))/120 + (thetax12*thetay1)/120 + (thetay1*thetax12)/60 + (thetay1*(thetax12 + thetax12))/120) },
{ ly0*(thetax12/8 + (thetax1*thetay1*thetaz1)/15 + thetax12/8 + thetax12/24 - 0.5) - lz0*((thetax1*(3*thetax12 + thetax12 + thetax12))/120 - thetax1/3 - (thetay1*thetaz1)/12 + (thetax1*thetax12)/12 + (thetax1*thetax12)/40 + (thetax1*(thetax12 + thetax12))/120) + lx0*((thetax12*thetaz1)/120 - (thetax1*thetay1)/6 - (thetaz1*(thetax12 + 3*thetax12 + thetax12))/120 + (thetax12*thetaz1)/60 + (thetaz1*(thetax12 + thetax12))/120),
 ly0*((thetaz1*thetax12)/10 + (thetax1*thetay1)/4 - thetaz1/3 + (thetaz1*(2*thetax12 + 2*thetax12))/120 + (thetaz1*(thetax12 + thetax12))/60) - lz0*((thetay1*(thetax12 + 3*thetax12 + thetax12))/60 - thetay1 - (thetax1*thetaz1)/12 + (thetax12*thetay1)/30 + (thetay1*thetax12)/15 + (thetay1*(4*thetax12 + 12*thetax12 + 2*thetax12))/120 + (thetay1*(thetax12 + thetax12))/60) - lx0*(thetax12/12 + thetax12/2 + thetax12/6 - 1),
 ly0*((thetay1*(thetax12 + thetax12 + 3*thetax12))/60 - thetay1/3 + (thetax1*thetaz1)/12 + (thetax12*thetay1)/60 + (thetay1*thetax12)/20 + thetay13/60) - lz0*((thetaz1*(thetax12 + 3*thetax12 + thetax12))/60 - thetaz1/3 - (thetax1*thetay1)/12 + (thetax12*thetaz1)/60 + (thetax12*thetaz1)/20 + thetaz13/60) - lx0*((thetax1*(thetax12 + 3*thetax12 + thetax12))/120 - (thetax1*(thetax12 + thetax12 + 3*thetax12))/120 + (thetay1*thetaz1)/3 - (thetax1*thetax12)/60 + (thetax1*thetax12)/60) },
{ ly0*((thetax1*(3*thetax12 + thetax12 + thetax12))/120 - thetax1/3 + (thetay1*thetaz1)/12 + (thetax1*thetax12)/40 + (thetax1*thetax12)/12 + (thetax1*(thetax12 + thetax12))/120) + lz0*(thetax12/8 - (thetax1*thetay1*thetaz1)/15 + thetax12/24 + thetax12/8 - 0.5) - lx0*((thetax1*thetaz1)/6 - (thetay1*(thetax12 + thetax12 + 3*thetax12))/120 + (thetax12*thetay1)/120 + (thetay1*thetax12)/60 + (thetay1*(thetax12 + thetax12))/120),
 ly0*((thetay1*(thetax12 + thetax12 + 3*thetax12))/60 - thetay1/3 + (thetax1*thetaz1)/12 + (thetax12*thetay1)/60 + (thetay1*thetax12)/20 + thetay13/60) - lz0*((thetaz1*(thetax12 + 3*thetax12 + thetax12))/60 - thetaz1/3 - (thetax1*thetay1)/12 + (thetax12*thetaz1)/60 + (thetax12*thetaz1)/20 + thetaz13/60) - lx0*((thetax1*(thetax12 + 3*thetax12 + thetax12))/120 - (thetax1*(thetax12 + thetax12 + 3*thetax12))/120 + (thetay1*thetaz1)/3 - (thetax1*thetax12)/60 + (thetax1*thetax12)/60),
 ly0*((thetaz1*(thetax12 + thetax12 + 3*thetax12))/60 - thetaz1 + (thetax1*thetay1)/12 + (thetax12*thetaz1)/30 + (thetax12*thetaz1)/15 + (thetaz1*(4*thetax12 + 2*thetax12 + 12*thetax12))/120 + (thetaz1*(thetax12 + thetax12))/60) - lz0*((thetay1*thetax12)/10 - (thetax1*thetaz1)/4 - thetay1/3 + (thetay1*(2*thetax12 + 2*thetax12))/120 + (thetay1*(thetax12 + thetax12))/60) - lx0*(thetax12/12 + thetax12/6 + thetax12/2 - 1) }
};
for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) J[3+i][3+j] = H[i][j];
} break;

case 1: {
double H[3][3] = {
{ lz0*((thetax1*(3*thetax12 + thetax12 + thetax12))/60 - thetax1 + (thetay1*thetaz1)/12 + (thetax1*thetax12)/30 + (thetax1*thetax12)/15 + (thetax1*(12*thetax12 + 4*thetax12 + 2*thetax12))/120 + (thetax1*(thetax12 + thetax12))/60) - lx0*((thetaz1*thetax12)/10 - (thetay1*thetax1)/4 - thetaz1/3 + (thetaz1*(2*thetax12 + 2*thetax12))/120 + (thetaz1*(thetax12 + thetax12))/60) - ly0*(thetax12/2 + thetax12/12 + thetax12/6 - 1),
 lz0*((thetay1*(thetax12 + 3*thetax12 + thetax12))/120 - thetay1/3 + (thetax1*thetaz1)/12 + (thetax12*thetay1)/12 + (thetay1*thetax12)/40 + (thetay1*(thetax12 + thetax12))/120) + lx0*(thetax12/8 - (thetax1*thetay1*thetaz1)/15 + thetax12/8 + thetax12/24 - 0.5) - ly0*((thetax1*thetay1)/6 - (thetaz1*(3*thetax12 + thetax12 + thetax12))/120 + (thetax12*thetaz1)/60 + (thetax12*thetaz1)/120 + (thetaz1*(thetax12 + thetax12))/120),
 lz0*((thetaz1*(3*thetax12 + thetax12 + thetax12))/60 - thetaz1/3 + (thetax1*thetay1)/12 + (thetax12*thetaz1)/20 + (thetax12*thetaz1)/60 + thetaz13/60) - lx0*((thetax1*(thetax12 + thetax12 + 3*thetax12))/60 - thetax1/3 - (thetay1*thetaz1)/12 + (thetax1*thetax12)/60 + (thetax1*thetax12)/20 + thetax13/60) - ly0*((thetay1*(thetax12 + thetax12 + 3*thetax12))/120 - (thetay1*(3*thetax12 + thetax12 + thetax12))/120 + (thetax1*thetaz1)/3 + (thetax12*thetay1)/60 - (thetay1*thetax12)/60) },
{ lz0*((thetay1*(thetax12 + 3*thetax12 + thetax12))/120 - thetay1/3 + (thetax1*thetaz1)/12 + (thetax12*thetay1)/12 + (thetay1*thetax12)/40 + (thetay1*(thetax12 + thetax12))/120) + lx0*(thetax12/8 - (thetax1*thetay1*thetaz1)/15 + thetax12/8 + thetax12/24 - 0.5) - ly0*((thetax1*thetay1)/6 - (thetaz1*(3*thetax12 + thetax12 + thetax12))/120 + (thetax12*thetaz1)/60 + (thetax12*thetaz1)/120 + (thetaz1*(thetax12 + thetax12))/120),
 lx0*(thetaz1/3 + (thetax1*thetay1)/4 - (thetax12*thetaz1)/60 - (thetaz1*(2*thetax12 + 12*thetax12 + 4*thetax12))/120) - ly0*(thetax12/12 + thetax12/12) + lz0*((thetay1*thetaz1)/4 - thetax1/3 + (thetax1*thetax12)/60 + (thetax1*(4*thetax12 + 12*thetax12 + 2*thetax12))/120),
 lz0*(thetax12/24 + (thetax1*thetay1*thetaz1)/15 + thetax12/8 + thetax12/8 - 0.5) - lx0*((thetay1*(thetax12 + 3*thetax12 + thetax12))/120 - thetay1/3 - (thetax1*thetaz1)/12 + (thetax12*thetay1)/40 + (thetay1*thetax12)/12 + (thetay1*(thetax12 + thetax12))/120) + ly0*((thetax1*thetax12)/120 - (thetay1*thetaz1)/6 - (thetax1*(thetax12 + thetax12 + 3*thetax12))/120 + (thetax1*thetax12)/60 + (thetax1*(thetax12 + thetax12))/120) },
{ lz0*((thetaz1*(3*thetax12 + thetax12 + thetax12))/60 - thetaz1/3 + (thetax1*thetay1)/12 + (thetax12*thetaz1)/20 + (thetax12*thetaz1)/60 + thetaz13/60) - lx0*((thetax1*(thetax12 + thetax12 + 3*thetax12))/60 - thetax1/3 - (thetay1*thetaz1)/12 + (thetax1*thetax12)/60 + (thetax1*thetax12)/20 + thetax13/60) - ly0*((thetay1*(thetax12 + thetax12 + 3*thetax12))/120 - (thetay1*(3*thetax12 + thetax12 + thetax12))/120 + (thetax1*thetaz1)/3 + (thetax12*thetay1)/60 - (thetay1*thetax12)/60),
 lz0*(thetax12/24 + (thetax1*thetay1*thetaz1)/15 + thetax12/8 + thetax12/8 - 0.5) - lx0*((thetay1*(thetax12 + 3*thetax12 + thetax12))/120 - thetay1/3 - (thetax1*thetaz1)/12 + (thetax12*thetay1)/40 + (thetay1*thetax12)/12 + (thetay1*(thetax12 + thetax12))/120) + ly0*((thetax1*thetax12)/120 - (thetay1*thetaz1)/6 - (thetax1*(thetax12 + thetax12 + 3*thetax12))/120 + (thetax1*thetax12)/60 + (thetax1*(thetax12 + thetax12))/120),
 lz0*((thetax1*thetax12)/10 + (thetay1*thetaz1)/4 - thetax1/3 + (thetax1*(2*thetax12 + 2*thetax12))/120 + (thetax1*(thetax12 + thetax12))/60) - lx0*((thetaz1*(thetax12 + thetax12 + 3*thetax12))/60 - thetaz1 - (thetax1*thetay1)/12 + (thetax12*thetaz1)/15 + (thetax12*thetaz1)/30 + (thetaz1*(2*thetax12 + 4*thetax12 + 12*thetax12))/120 + (thetaz1*(thetax12 + thetax12))/60) - ly0*(thetax12/6 + thetax12/12 + thetax12/2 - 1) }
};
for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) J[3+i][3+j] = H[i][j];
} break; 

case 2: {
double H[3][3] = {
{ lx0*((thetay1*thetax12)/10 + (thetaz1*thetax1)/4 - thetay1/3 + (thetay1*(2*thetax12 + 2*thetax12))/120 + (thetay1*(thetax12 + thetax12))/60) - ly0*((thetax1*(3*thetax12 + thetax12 + thetax12))/60 - thetax1 - (thetay1*thetaz1)/12 + (thetax1*thetax12)/15 + (thetax1*thetax12)/30 + (thetax1*(12*thetax12 + 2*thetax12 + 4*thetax12))/120 + (thetax1*(thetax12 + thetax12))/60) - lz0*(thetax12/2 + thetax12/6 + thetax12/12 - 1),
 lx0*((thetax1*(thetax12 + 3*thetax12 + thetax12))/60 - thetax1/3 + (thetay1*thetaz1)/12 + (thetax1*thetax12)/20 + (thetax1*thetax12)/60 + thetax13/60) - ly0*((thetay1*(3*thetax12 + thetax12 + thetax12))/60 - thetay1/3 - (thetax1*thetaz1)/12 + (thetax12*thetay1)/20 + (thetay1*thetax12)/60 + thetay13/60) - lz0*((thetaz1*(3*thetax12 + thetax12 + thetax12))/120 - (thetaz1*(thetax12 + 3*thetax12 + thetax12))/120 + (thetax1*thetay1)/3 - (thetax12*thetaz1)/60 + (thetax12*thetaz1)/60),
 lx0*(thetax12/8 + (thetax1*thetay1*thetaz1)/15 + thetax12/24 + thetax12/8 - 0.5) - ly0*((thetaz1*(thetax12 + thetax12 + 3*thetax12))/120 - thetaz1/3 - (thetax1*thetay1)/12 + (thetax12*thetaz1)/12 + (thetax12*thetaz1)/40 + (thetaz1*(thetax12 + thetax12))/120) + lz0*((thetax12*thetay1)/60 - (thetax1*thetaz1)/6 - (thetay1*(3*thetax12 + thetax12 + thetax12))/120 + (thetay1*thetax12)/120 + (thetay1*(thetax12 + thetax12))/120) },
{ lx0*((thetax1*(thetax12 + 3*thetax12 + thetax12))/60 - thetax1/3 + (thetay1*thetaz1)/12 + (thetax1*thetax12)/20 + (thetax1*thetax12)/60 + thetax13/60) - ly0*((thetay1*(3*thetax12 + thetax12 + thetax12))/60 - thetay1/3 - (thetax1*thetaz1)/12 + (thetax12*thetay1)/20 + (thetay1*thetax12)/60 + thetay13/60) - lz0*((thetaz1*(3*thetax12 + thetax12 + thetax12))/120 - (thetaz1*(thetax12 + 3*thetax12 + thetax12))/120 + (thetax1*thetay1)/3 - (thetax12*thetaz1)/60 + (thetax12*thetaz1)/60),
 lx0*((thetay1*(thetax12 + 3*thetax12 + thetax12))/60 - thetay1 + (thetax1*thetaz1)/12 + (thetax12*thetay1)/15 + (thetay1*thetax12)/30 + (thetay1*(2*thetax12 + 12*thetax12 + 4*thetax12))/120 + (thetay1*(thetax12 + thetax12))/60) - ly0*((thetax1*thetax12)/10 - (thetaz1*thetay1)/4 - thetax1/3 + (thetax1*(2*thetax12 + 2*thetax12))/120 + (thetax1*(thetax12 + thetax12))/60) - lz0*(thetax12/6 + thetax12/2 + thetax12/12 - 1),
 lx0*((thetaz1*(thetax12 + thetax12 + 3*thetax12))/120 - thetaz1/3 + (thetax1*thetay1)/12 + (thetax12*thetaz1)/40 + (thetax12*thetaz1)/12 + (thetaz1*(thetax12 + thetax12))/120) + ly0*(thetax12/24 - (thetax1*thetay1*thetaz1)/15 + thetax12/8 + thetax12/8 - 0.5) - lz0*((thetay1*thetaz1)/6 - (thetax1*(thetax12 + 3*thetax12 + thetax12))/120 + (thetax1*thetax12)/60 + (thetax1*thetax12)/120 + (thetax1*(thetax12 + thetax12))/120) },
{ lx0*(thetax12/8 + (thetax1*thetay1*thetaz1)/15 + thetax12/24 + thetax12/8 - 0.5) - ly0*((thetaz1*(thetax12 + thetax12 + 3*thetax12))/120 - thetaz1/3 - (thetax1*thetay1)/12 + (thetax12*thetaz1)/12 + (thetax12*thetaz1)/40 + (thetaz1*(thetax12 + thetax12))/120) + lz0*((thetax12*thetay1)/60 - (thetax1*thetaz1)/6 - (thetay1*(3*thetax12 + thetax12 + thetax12))/120 + (thetay1*thetax12)/120 + (thetay1*(thetax12 + thetax12))/120),
 lx0*((thetaz1*(thetax12 + thetax12 + 3*thetax12))/120 - thetaz1/3 + (thetax1*thetay1)/12 + (thetax12*thetaz1)/40 + (thetax12*thetaz1)/12 + (thetaz1*(thetax12 + thetax12))/120) + ly0*(thetax12/24 - (thetax1*thetay1*thetaz1)/15 + thetax12/8 + thetax12/8 - 0.5) - lz0*((thetay1*thetaz1)/6 - (thetax1*(thetax12 + 3*thetax12 + thetax12))/120 + (thetax1*thetax12)/60 + (thetax1*thetax12)/120 + (thetax1*(thetax12 + thetax12))/120),
 lx0*((thetax1*thetaz1)/4 - thetay1/3 + (thetax12*thetay1)/60 + (thetay1*(2*thetax12 + 4*thetax12 + 12*thetax12))/120) - lz0*(thetax12/12 + thetax12/12) + ly0*(thetax1/3 + (thetay1*thetaz1)/4 - (thetax1*thetax12)/60 - (thetax1*(4*thetax12 + 2*thetax12 + 12*thetax12))/120) }
};
for(int i=0; i<3; ++i) for(int j=0; j<3; ++j) J[3+i][3+j] = H[i][j];
} break;

}

#else
  J.zero();
#endif
}

