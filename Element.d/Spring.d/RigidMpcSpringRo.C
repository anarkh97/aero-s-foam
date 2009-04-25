#include <Element.d/Spring.d/RigidMpcSpringRo.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Domain.h>
extern Domain *domain;

// this element is ties two nodes so that they both have the same translations/rotations
RigidMpcSpringRo::RigidMpcSpringRo(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nmpc = 0; 
}

void
RigidMpcSpringRo::setProp(StructProp *p, bool _myProp) 
{ 
  prop = p; 
  myProp = _myProp; 
  if(prop->A != 0.0) nmpc++;
  if(prop->E != 0.0) nmpc++;
  if(prop->nu != 0.0) nmpc++;
  if(nmpc > 0) {
    mpc = new LMPCons * [nmpc];
    lambda = new double[nmpc];
    glMpcNb = new int[nmpc];
  }
}

void 
RigidMpcSpringRo::computeMPCs(CoordSet &cs, int &lmpcnum)
{
  double rhs = 0.0;
  int i;
  for(i=0; i<nmpc; ++i) {
    lmpcnum++;
    mpc[i] = new LMPCons(lmpcnum, rhs);
  }

  // the rotation is equal at both ends (nmpc constraints)
  int count = 0;
  if(prop->A != 0.0) {
    LMPCTerm term1(nn[0],3,1.0);
    mpc[count]->addterm(&term1);
    LMPCTerm term2(nn[0],3,-1.0);
    mpc[count]->addterm(&term2);
    count++;
  }
  if(prop->E != 0.0) {
    LMPCTerm term1(nn[0],4,1.0);
    mpc[count]->addterm(&term1);
    LMPCTerm term2(nn[0],4,-1.0);
    mpc[count]->addterm(&term2);
    count++;
  }
  if(prop->nu != 0.0) {
    LMPCTerm term1(nn[0],5,1.0);
    mpc[count]->addterm(&term1);
    LMPCTerm term2(nn[0],5,-1.0);
    mpc[count]->addterm(&term2);
    count++;
  }

  for(i=0; i<count; ++i) 
    glMpcNb[i] = domain->addLMPC(mpc[i]);
}

void
RigidMpcSpringRo::updateMPCs(GeomState &gState)
{
  /* nothing to update for this element */
}

void
RigidMpcSpringRo::renum(int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
}

FullSquareMatrix
RigidMpcSpringRo::massMatrix(CoordSet &cs, double *mel, int cmflg)
{
  FullSquareMatrix elementMassMatrix(12, mel);
  elementMassMatrix.zero();
  return elementMassMatrix;
}

FullSquareMatrix
RigidMpcSpringRo::stiffness(CoordSet &cs, double *k, int flg)
{
  FullSquareMatrix ret(12, k);
  ret.zero();
  return ret;
}

int
RigidMpcSpringRo::numNodes()
{
  return 2;
}

int *
RigidMpcSpringRo::nodes(int *p)
{
  if(p == 0) p = new int[2];
  p[0] = nn[0];
  p[1] = nn[1];
  return p;
}

int
RigidMpcSpringRo::numDofs()
{
  return 12;
}

int *
RigidMpcSpringRo::dofs(DofSetArray &dsa, int *p)
{
  if(p == 0) p = new int[12];

  dsa.number(nn[0], DofSet::XYZdisp | DofSet::XYZrot, p  );
  dsa.number(nn[1], DofSet::XYZdisp | DofSet::XYZrot, p+6);
  return p;
}

void
RigidMpcSpringRo::markDofs(DofSetArray &dsa)
{
  dsa.mark(nn, 2, DofSet::XYZdisp | DofSet::XYZrot);
}

void
RigidMpcSpringRo::getStiffAndForce(GeomState &gState, CoordSet &cs,
                                       FullSquareMatrix &Ktan, double *f)
{
  Ktan.zero();
  for(int i=0; i<12; ++i) f[i] = 0.0;
}

