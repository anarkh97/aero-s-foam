#include <Element.d/Spring.d/RigidMpcSpring.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>

// this element is ties two nodes so that they both have the same translations/rotations
RigidMpcSpring::RigidMpcSpring(int* _nn)
  : RigidMpcElement(2, 3, DofSet::XYZdisp, 6, _nn)
{
  /* do nothing */
}

void 
RigidMpcSpring::computeMPCs(CoordSet &cs)
{
  for(int i = 0; i < nMpcs; ++i)
    mpcs[i] = new LMPCons(0, 0.0);

  // the translations & rotation is equal at both ends (6 constraints)
  for(int i = 0; i < nMpcs; ++i) {
    mpcs[i]->addterm(new LMPCTerm(nn[0], i, 1.0));
    mpcs[i]->addterm(new LMPCTerm(nn[1], i, -1.0));
  }
}

int
RigidMpcSpring::getTopNumber()
{
  return 2;
}

