#include <Element.d/Joint.d/SphericalJoint.h>
#include <Utils.d/dofset.h>
#include <Driver.d/Mpc.h>
#include <Corotational.d/utilities.h>

SphericalJoint::SphericalJoint(int* _nn)
 : RigidMpcElement(2, 3, DofSet::XYZdisp, 3, _nn)
{
}

void 
SphericalJoint::computeMPCs(CoordSet& cs)
{
  Node &nd1 = cs.getNode(nn[0]);
  Node &nd2 = cs.getNode(nn[1]);

  d0[0] = nd2.x - nd1.x;
  d0[1] = nd2.y - nd1.y;
  d0[2] = nd2.z - nd1.z;

  for(int i = 0; i < nMpcs; ++i)
    mpcs[i] = new LMPCons(0, d0[i]);

  // translation in x, y, z (3 constraint equations)
  for(int i = 0; i < 3; ++i) {
    mpcs[i]->addterm(new LMPCTerm(nn[0], i, -1.0));
    mpcs[i]->addterm(new LMPCTerm(nn[1], i, 1.0));
  }
}

int 
SphericalJoint::getTopNumber() 
{ 
  return 106; 
}

void 
SphericalJoint::updateLMPCs(GeomState& gState, CoordSet& cs)
{
  // nodes' current coordinates
  NodeState ns1 = gState[nn[0]];
  NodeState ns2 = gState[nn[1]];

  double dx = ns2.x - ns1.x;
  double dy = ns2.y - ns1.y;
  double dz = ns2.z - ns1.z;

  double d[3] = { dx, dy, dz };

  // values of constraint functions
  for(int i=0; i<3; ++i) mpcs[i]->rhs.r_value = d[i];
}

