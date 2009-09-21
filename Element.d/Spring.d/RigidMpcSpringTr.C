#include <Element.d/Spring.d/RigidMpcSpringTr.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>

// this element is ties two nodes so that they both have the same rotations
RigidMpcSpringTr::RigidMpcSpringTr(int* _nn)
  : RigidMpcElement(2, 3, DofSet::XYZdisp, 3, _nn)
{
  /* do nothing */
}

void 
RigidMpcSpringTr::computeMPCs(CoordSet &cs)
{
  for(int i = 0; i < nMpcs; ++i)
    mpcs[i] = new LMPCons(0, 0.0);

  // the translation is equal at both ends
  int count = 0;
  if(prop->kx != 0.0) {
    mpcs[count]->addterm(new LMPCTerm(nn[0], 0, 1.0));
    mpcs[count]->addterm(new LMPCTerm(nn[1], 0, -1.0));
    count++;
  }
  if(prop->ky != 0.0) {
    mpcs[count]->addterm(new LMPCTerm(nn[0], 1, 1.0));
    mpcs[count]->addterm(new LMPCTerm(nn[1], 1, -1.0));
    count++;
  }
  if(prop->kz != 0.0) {
    mpcs[count]->addterm(new LMPCTerm(nn[0], 2, 1.0));
    mpcs[count]->addterm(new LMPCTerm(nn[0], 2, -1.0));
    count++;
  }
  nMpcs = count;
}

int
RigidMpcSpringTr::getTopNumber() 
{ 
  return 101; 
}

bool 
RigidMpcSpringTr::isSafe() 
{ 
  return false; 
}
