#include <Element.d/Truss.d/TwoNodeTrussRigidMpc.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>
#include <Corotational.d/utilities.h>
#define FOLLOWER_FORCE

TwoNodeTrussRigidMpc::TwoNodeTrussRigidMpc(int* _nn)
 : RigidMpcElement(2, 3, DofSet::XYZdisp, 1, _nn)
{
  /* do nothing */
}

void 
TwoNodeTrussRigidMpc::computeMPCs(CoordSet &cs)
{
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double lx = nd1.x - nd2.x;
  double ly = nd1.y - nd2.y;
  double lz = nd1.z - nd2.z;

  double l = sqrt( lx*lx + ly*ly + lz*lz );

  if(l == 0.0) {
    cerr << " *** ERROR: Rigid truss has zero length, nodes: " << nn[0]+1 << " " << nn[1]+1 << ". Exiting...\n";
    exit(-1);
  }
  else {
    mpcs[0] = new LMPCons(0, 0.0);
 
    double c1[3];
    c1[0] = lx/l;
    c1[1] = ly/l;
    c1[2] = lz/l;

    // translation in x, y, z (1 constraint equation)
    for(int i = 0; i < 3; ++i)
      mpcs[0]->addterm(new LMPCTerm(nn[0], i, c1[i]));
    for(int i = 0; i < 3; ++i)
      mpcs[0]->addterm(new LMPCTerm(nn[1], i, -c1[i]));
  }
}

int
TwoNodeTrussRigidMpc::getTopNumber()
{
  return 101;
}

bool
TwoNodeTrussRigidMpc::isSafe()
{
  return false;
}

int
TwoNodeTrussRigidMpc::buildFrame(CoordSet& cs)
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
TwoNodeTrussRigidMpc::computePressureForce(CoordSet& cs, Vector& elPressureForce,
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

  elPressureForce[0]  = px;
  elPressureForce[1]  = py;
  elPressureForce[2]  = pz;
  elPressureForce[3]  = px;
  elPressureForce[4]  = py;
  elPressureForce[5]  = pz;
  if(mode == 1) elPressureForce[6] = 0.0;
}

void
TwoNodeTrussRigidMpc::updTransMatrix(CoordSet& cs, GeomState* geomState, double t0n[3][3], double& length)
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
TwoNodeTrussRigidMpc::getLength(CoordSet& cs, double& length)
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

/*
void
TwoNodeTrussRigidMpc::getStiffAndForce(GeomState& gState, CoordSet& cs,
                                       FullSquareMatrix& Ktan, double* f)
{
  cerr << "here in TwoNodeTrussRigidMpc::getStiffAndForce\n";
  // nodes' original coordinates
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double l0x = nd1.x-nd2.x;
  double l0y = nd1.y-nd2.y;
  double l0z = nd1.z-nd2.z;
  double l0 = sqrt(l0x*l0x + l0y*l0y + l0z*l0z);

  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  double lx = ns1.x-ns2.x;
  double ly = ns1.y-ns2.y;
  double lz = ns1.z-ns2.z;
  double l = sqrt(lx*lx + ly*ly + lz*lz);

  // current value of constraint's lagrange multiplier
  double lambda = gState[nn[2]].x;

  // let c = the constraint function: c(x) = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 ) - l0
  // let g = gradient of constraint function
  //     H = hessian of constraint function

  // from SQP theory, the contribution of the constraint to the "internal force vector" is [ lambda*g ]
  //                                                                                       [ c        ]
  f[0] = lambda*(lx/l);
  f[1] = lambda*(ly/l);
  f[2] = lambda*(lz/l);
  f[3] = -f[0];
  f[4] = -f[1];
  f[5] = -f[2];
  f[6] = l-l0;

  // and to the "tangent stiffness" is [ lambda*H  g ]
  //                                   [ g^T         ] (note: lambda*h is a scalar multiplication)
  double l2 = l*l;
  double l3 = l*l*l;
  Ktan[0][0] = Ktan[3][3] = lambda*((l2-lx*lx)/l3);
  Ktan[1][1] = Ktan[4][4] = lambda*((l2-ly*ly)/l3);
  Ktan[2][2] = Ktan[5][5] = lambda*((l2-lz*lz)/l3);

  Ktan[0][1] = Ktan[1][0] = Ktan[3][4] = Ktan[4][3] = lambda*(-lx*ly/l3);
  Ktan[0][2] = Ktan[2][0] = Ktan[3][5] = Ktan[5][3] = lambda*(-lx*lz/l3);
  Ktan[1][2] = Ktan[2][1] = Ktan[4][5] = Ktan[5][4] = lambda*(-ly*lz/l3);

  Ktan[0][3] = Ktan[3][0] = -Ktan[0][0];
  Ktan[1][4] = Ktan[4][1] = -Ktan[1][1];
  Ktan[2][5] = Ktan[5][2] = -Ktan[2][2];

  Ktan[0][4] = Ktan[4][0] = Ktan[1][3] = Ktan[3][1] = -Ktan[0][1];
  Ktan[0][5] = Ktan[5][0] = Ktan[2][3] = Ktan[3][2] = -Ktan[0][2];
  Ktan[1][5] = Ktan[5][1] = Ktan[2][4] = Ktan[4][2] = -Ktan[1][2];
  
  Ktan[6][0] = Ktan[0][6] = lx/l;
  Ktan[6][1] = Ktan[1][6] = ly/l;
  Ktan[6][2] = Ktan[2][6] = lz/l;
  Ktan[6][3] = Ktan[3][6] = -Ktan[0][6];
  Ktan[6][4] = Ktan[4][6] = -Ktan[1][6];
  Ktan[6][5] = Ktan[5][6] = -Ktan[2][6];
}
*/

void
TwoNodeTrussRigidMpc::updateLMPCs(GeomState& gState, CoordSet& cs)
{
  // nodes' original coordinates
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  double lx0 = nd1.x - nd2.x;
  double ly0 = nd1.y - nd2.y;
  double lz0 = nd1.z - nd2.z;
  double l0 = sqrt( lx0*lx0 + ly0*ly0 + lz0*lz0 );

  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  double lx = ns1.x - ns2.x;
  double ly = ns1.y - ns2.y;
  double lz = ns1.z - ns2.z;
  double l = sqrt(lx*lx + ly*ly + lz*lz);

  double c1[3];
  c1[0] = lx/l;
  c1[1] = ly/l;
  c1[2] = lz/l;

  for(int i = 0; i < 3; ++i) {
    mpcs[0]->terms[i].coef.r_value = c1[i];
    mpcs[0]->terms[3+i].coef.r_value = -c1[i];
  }
  mpcs[0]->rhs.r_value = l - l0;
}

void 
TwoNodeTrussRigidMpc::getHessian(GeomState& gState, CoordSet&, int, FullSquareMatrix& H)
{
  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  double lx = ns1.x-ns2.x;
  double ly = ns1.y-ns2.y;
  double lz = ns1.z-ns2.z;
  double l = sqrt(lx*lx + ly*ly + lz*lz);
  
  double l2 = l*l;
  double l3 = l*l*l;
  H[0][0] = H[3][3] = (l2-lx*lx)/l3;
  H[1][1] = H[4][4] = (l2-ly*ly)/l3;
  H[2][2] = H[5][5] = (l2-lz*lz)/l3;

  H[0][1] = H[1][0] = H[3][4] = H[4][3] = -lx*ly/l3;
  H[0][2] = H[2][0] = H[3][5] = H[5][3] = -lx*lz/l3;
  H[1][2] = H[2][1] = H[4][5] = H[5][4] = -ly*lz/l3;

  H[0][3] = H[3][0] = -H[0][0];
  H[1][4] = H[4][1] = -H[1][1];
  H[2][5] = H[5][2] = -H[2][2];

  H[0][4] = H[4][0] = H[1][3] = H[3][1] = -H[0][1];
  H[0][5] = H[5][0] = H[2][3] = H[3][2] = -H[0][2];
  H[1][5] = H[5][1] = H[2][4] = H[4][2] = -H[1][2];
}
