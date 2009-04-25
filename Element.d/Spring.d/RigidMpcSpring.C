#include <Element.d/Spring.d/RigidMpcSpring.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Domain.h>
extern Domain *domain;

// this element is ties two nodes so that they both have the same translations/rotations
RigidMpcSpring::RigidMpcSpring(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
}

void 
RigidMpcSpring::computeMPCs(CoordSet &cs, int &lmpcnum)
{
  double rhs = 0.0;
  int i;
  for(i=0; i<6; ++i) {
    lmpcnum++;
    mpc[i] = new LMPCons(lmpcnum, rhs);
  }

  // the translations & rotation is equal at both ends (6 constraints)
  for(i=0; i<6; ++i) {
    LMPCTerm term1(nn[0],i,1.0);
    mpc[i]->addterm(&term1);
    LMPCTerm term2(nn[0],i,-1.0);
    mpc[i]->addterm(&term2);
  }

  for(i=0; i<6; ++i) 
    glMpcNb[i] = domain->addLMPC(mpc[i]);
}

void
RigidMpcSpring::updateMPCs(GeomState &gState)
{
  /* nothing to update for this element */
}

void
RigidMpcSpring::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
}

FullSquareMatrix
RigidMpcSpring::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  FullSquareMatrix elementMassMatrix(12, mel);
  elementMassMatrix.zero();
  return elementMassMatrix;
}

FullSquareMatrix
RigidMpcSpring::stiffness(CoordSet &cs, double *k, int flg)
{
  FullSquareMatrix ret(12, k);
  ret.zero();
  return ret;
}

int
RigidMpcSpring::numNodes()
{
  return 2;
}

int *
RigidMpcSpring::nodes(int *p)
{
  if(p == 0) p = new int[2];
  p[0] = nn[0];
  p[1] = nn[1];
  return p;
}

int
RigidMpcSpring::numDofs()
{
  return 12;
}

int *
RigidMpcSpring::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[12];

  dsa.number(nn[0], DofSet::XYZdisp | DofSet::XYZrot, p  );
  dsa.number(nn[1], DofSet::XYZdisp | DofSet::XYZrot, p+6);
  return p;
}

void
RigidMpcSpring::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot);
}

void
RigidMpcSpring::getStiffAndForce(GeomState &gState, CoordSet &cs,
                                       FullSquareMatrix &Ktan, double *f)
{
  Ktan.zero();
  for(int i=0; i<12; ++i) f[i] = 0.0;
}

