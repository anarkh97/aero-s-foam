#ifndef _SPHERICALJOINT_H_
#define _SPHERICALJOINT_H_

#include <Element.d/MpcElement.d/RigidMpcElement.h>

class SphericalJoint : public RigidMpcElement
{
    double d0[3];
  public:
    SphericalJoint(int*);
    void computeMPCs(CoordSet&);
    int getTopNumber();
    void updateLMPCs(GeomState& gState, CoordSet& cs);
};

#endif
