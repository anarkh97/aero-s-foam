#ifndef _RIGIDSOLID6DOF_H_
#define _RIGIDSOLID6DOF_H_

#include <Element.d/SuperElement.h>

class RigidSolid6Dof : public SuperElement
{
  public:
    RigidSolid6Dof(int numnodes, int *nodenums);

    Element* clone();
    int getTopNumber() { return 101; }
    int numTopNodes() { return 2; }

};

#endif

