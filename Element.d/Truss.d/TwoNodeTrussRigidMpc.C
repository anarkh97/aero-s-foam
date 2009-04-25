#include <Element.d/Truss.d/TwoNodeTrussRigidMpc.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Domain.h>
extern Domain *domain;

TwoNodeTrussRigidMpc::TwoNodeTrussRigidMpc(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
}

void 
TwoNodeTrussRigidMpc::computeMPCs(CoordSet &cs, int &lmpcnum)
{
  double rhs = 0.0;
  lmpcnum++;
  mpc = new LMPCons(lmpcnum, rhs);

  Node &node1 = cs.getNode(nn[0]);
  Node &node2 = cs.getNode(nn[1]);
  double l0[3];
  l0[0] = node1.x-node2.x;
  l0[1] = node1.y-node2.y;
  l0[2] = node1.z-node2.z;
  len0 = sqrt(l0[0]*l0[0] + l0[1]*l0[1] + l0[2]*l0[2]);

  int i;
  for(i=0; i<3; ++i) {
    LMPCTerm term1(nn[0], i, l0[i]);
    mpc->addterm(&term1);
    LMPCTerm term2(nn[1], i, -l0[i]);
    mpc->addterm(&term2);
  }
  glMpcNb = domain->addLMPC(mpc);
}

void
TwoNodeTrussRigidMpc::updateMPCs(GeomState &gState)
{
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  // compute current length
  double l[3];
  l[0] = ns1.x-ns2.x;
  l[1] = ns1.y-ns2.y;
  l[2] = ns1.z-ns2.z;
  double len = sqrt(l[0]*l[0] + l[1]*l[1] + l[2]*l[2]);

  int i;
  for(i=0; i<3; ++i) {
    mpc->rhs.r_value = len0-len;  // check this
    mpc->terms[2*i].coef.r_value = l[i]; 
    mpc->terms[2*i+1].coef.r_value = -l[i];
  }
}

void
TwoNodeTrussRigidMpc::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
}

FullSquareMatrix
TwoNodeTrussRigidMpc::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  FullSquareMatrix elementMassMatrix(6,mel);
  elementMassMatrix.zero();
  return elementMassMatrix;
}

FullSquareMatrix
TwoNodeTrussRigidMpc::stiffness(CoordSet &cs, double *k, int flg)
{
  FullSquareMatrix ret(6,k);
  ret.zero();
  return ret;
}

int
TwoNodeTrussRigidMpc::numNodes()
{
  return 2;
}

int *
TwoNodeTrussRigidMpc::nodes(int *p)
{
  if(p == 0) p = new int[2];
  p[0] = nn[0];
  p[1] = nn[1];
  return p;
}

int
TwoNodeTrussRigidMpc::numDofs()
{
  return 6;
}

int *
TwoNodeTrussRigidMpc::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[6];

  dsa.number(nn[0], DofSet::XYZdisp, p);
  dsa.number(nn[1], DofSet::XYZdisp, p+3);
  return p;
}

void
TwoNodeTrussRigidMpc::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, 2, DofSet::XYZdisp);
}

void
TwoNodeTrussRigidMpc::getStiffAndForce(GeomState &gState, CoordSet &cs,
                                       FullSquareMatrix &Ktan, double *f)
{
  Ktan.zero();
  for(int i=0; i<numDofs(); ++i) f[i] = 0.0;
/* PJSA: this is done by the MpcElements
  // Get Nodes original coordinates
  Node &node1 = cs.getNode(nn[0]);
  Node &node2 = cs.getNode(nn[1]);
  // Get Nodes current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];
  
  // compute current length
  double lx = ns1.x-ns2.x;
  double ly = ns1.y-ns2.y;
  double lz = ns1.z-ns2.z;
  double l = sqrt(lx*lx + ly*ly + lz*lz);
  double coef = 1.0/l;

  // (internal) force vector
  f[0] = coef*lx*lambda;
  f[1] = coef*ly*lambda;
  f[2] = coef*lz*lambda;
  f[3] = -f[0];
  f[4] = -f[1];
  f[5] = -f[2];
 
  // tangent stiffness
  Ktan.zero();

  // geometric stiffness constribution 
  double coef1 = lambda/(l*l*l);
  double l2 = l*l;
  Ktan[0][0] = Ktan[3][3] = coef1*(l2-lx*lx); 
  Ktan[1][1] = Ktan[4][4] = coef1*(l2-ly*ly); 
  Ktan[2][2] = Ktan[5][5] = coef1*(l2-lz*lz); 
 
  Ktan[0][1] = Ktan[1][0] = Ktan[3][4] = Ktan[4][3] = -coef1*lx*ly;
  Ktan[0][2] = Ktan[2][0] = Ktan[3][5] = Ktan[5][3] = -coef1*lx*lz;
  Ktan[1][2] = Ktan[2][1] = Ktan[4][5] = Ktan[5][4] = -coef1*ly*lz;

  Ktan[0][3] = Ktan[3][0] = -Ktan[0][0];
  Ktan[1][4] = Ktan[4][1] = -Ktan[1][1];
  Ktan[2][5] = Ktan[5][2] = -Ktan[2][2];

  Ktan[0][4] = Ktan[4][0] = Ktan[1][3] = Ktan[3][1] = -Ktan[0][1];
  Ktan[0][5] = Ktan[5][0] = Ktan[2][3] = Ktan[3][2] = -Ktan[0][2];
  Ktan[1][5] = Ktan[5][1] = Ktan[2][4] = Ktan[4][2] = -Ktan[1][2];
*/
}

