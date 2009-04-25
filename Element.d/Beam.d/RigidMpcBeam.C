#include <Element.d/Beam.d/RigidMpcBeam.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Domain.h>
extern Domain *domain;

RigidMpcBeam::RigidMpcBeam(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
}

void 
RigidMpcBeam::computeMPCs(CoordSet &cs, int &lmpcnum)
{
  double rhs = 0.0;
  int i;
  for(i=0; i<6; ++i) {
    lmpcnum++;
    mpc[i] = new LMPCons(lmpcnum, rhs);
  }

  Node &node1 = cs.getNode(nn[0]);
  Node &node2 = cs.getNode(nn[1]);
  double l0[3];
  l0[0] = node1.x-node2.x;
  l0[1] = node1.y-node2.y;
  l0[2] = node1.z-node2.z;
  len0 = sqrt(l0[0]*l0[0] + l0[1]*l0[1] + l0[2]*l0[2]);

  // the rotation is equal at both ends (3 constraints)
  for(i=3; i<6; ++i) {
    LMPCTerm term1(nn[0],i,1.0);
    mpc[i]->addterm(&term1);
    LMPCTerm term2(nn[1],i,-1.0);
    mpc[i]->addterm(&term2);
  }
  if(len0 < 1.0e-12) {
    // nodes are coincident, tie them together
    cerr << " *** WARNING: tying coincident nodes in RigidMpcBeam element: length = " << len0 << " < tol (1.0e-12) \n";
    for(i=0; i<3; ++i) {
      LMPCTerm term1(nn[0],i,1.0);
      mpc[i]->addterm(&term1);
      LMPCTerm term2(nn[1],i,-1.0);
      mpc[i]->addterm(&term2);
    }
  }
  else {
    // non-axial rotation and translation are consistent (3 constraints)
    // Ux(2)-Ux(1) + Ly*Rz(1)-Lz*Ry(1) 
    LMPCTerm term01(nn[1],0,1.0);
    mpc[0]->addterm(&term01);
    LMPCTerm term02(nn[0],0,-1.0);
    mpc[0]->addterm(&term02);
    LMPCTerm term03(nn[0],5,l0[1]);
    mpc[0]->addterm(&term03);
    LMPCTerm term04(nn[0],4,-l0[2]);
    mpc[0]->addterm(&term04);
    // Uy(2)-Uy(1) + Lz*Rx(1)-Lx*Rz(1)
    LMPCTerm term11(nn[1],1,1.0);
    mpc[1]->addterm(&term11);
    LMPCTerm term12(nn[0],1,-1.0);
    mpc[1]->addterm(&term12);
    LMPCTerm term13(nn[0],3,l0[2]);
    mpc[1]->addterm(&term13);
    LMPCTerm term14(nn[0],5,-l0[0]);
    mpc[1]->addterm(&term14);
    // Uz(2)-Uz(1) + Lx*Ry(1)-Ly*Rx(1)
    LMPCTerm term21(nn[1],2,1.0);
    mpc[2]->addterm(&term21);
    LMPCTerm term22(nn[0],2,-1.0);
    mpc[2]->addterm(&term22);
    LMPCTerm term23(nn[0],4,l0[0]);
    mpc[2]->addterm(&term23);
    LMPCTerm term24(nn[0],3,-l0[1]);
    mpc[2]->addterm(&term24);
  }

  for(i=0; i<6; ++i) 
    glMpcNb[i] = domain->addLMPC(mpc[i]);
}

void
RigidMpcBeam::updateMPCs(GeomState &gState)
{
  cerr << " *** WARNING: RigidMpcBeam::updateMPCs(...) is not implemented \n";
}

void
RigidMpcBeam::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
}

FullSquareMatrix
RigidMpcBeam::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  FullSquareMatrix elementMassMatrix(12, mel);
  elementMassMatrix.zero();
  return elementMassMatrix;
}

FullSquareMatrix
RigidMpcBeam::stiffness(CoordSet &cs, double *k, int flg)
{
  FullSquareMatrix ret(12, k);
  ret.zero();
  return ret;
}

int
RigidMpcBeam::numNodes()
{
  return 2;
}

int *
RigidMpcBeam::nodes(int *p)
{
  if(p == 0) p = new int[2];
  p[0] = nn[0];
  p[1] = nn[1];
  return p;
}

int
RigidMpcBeam::numDofs()
{
  return 12;
}

int *
RigidMpcBeam::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[12];

  dsa.number(nn[0], DofSet::XYZdisp | DofSet::XYZrot, p  );
  dsa.number(nn[1], DofSet::XYZdisp | DofSet::XYZrot, p+6);
  return p;
}

void
RigidMpcBeam::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot);
}

void
RigidMpcBeam::getStiffAndForce(GeomState &gState, CoordSet &cs,
                                       FullSquareMatrix &Ktan, double *f)
{
  Ktan.zero();
  for(int i=0; i<12; ++i) f[i] = 0.0;
}

