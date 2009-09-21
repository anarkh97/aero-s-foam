#include <Element.d/Spring.d/RigidMpcSpringRo.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>

// this element is ties two nodes so that they both have the same rotations
RigidMpcSpringRo::RigidMpcSpringRo(int* _nn)
  : RigidMpcElement(2, 3, DofSet::XYZrot, 3, _nn)
{
  /* do nothing */
}

void 
RigidMpcSpringRo::computeMPCs(CoordSet &cs)
{
  for(int i = 0; i < nMpcs; ++i)
    mpcs[i] = new LMPCons(0, 0.0);

  // the rotation is equal at both ends
  int count = 0;
  if(prop->kx != 0.0) {
    mpcs[count]->addterm(new LMPCTerm(nn[0], 3, 1.0));
    mpcs[count]->addterm(new LMPCTerm(nn[1], 3, -1.0));
    count++;
  }
  if(prop->ky != 0.0) {
    mpcs[count]->addterm(new LMPCTerm(nn[0], 4, 1.0));
    mpcs[count]->addterm(new LMPCTerm(nn[1], 4, -1.0));
    count++;
  }
  if(prop->kz != 0.0) {
    mpcs[count]->addterm(new LMPCTerm(nn[0], 5, 1.0));
    mpcs[count]->addterm(new LMPCTerm(nn[0], 5, -1.0));
    count++;
  }
  nMpcs = count;
}

int
RigidMpcSpringRo::getTopNumber() 
{ 
  return 101; 
}

bool 
RigidMpcSpringRo::isSafe() 
{ 
  return false; 
}
