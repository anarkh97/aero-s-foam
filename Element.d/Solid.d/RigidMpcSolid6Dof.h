#ifndef _RIGIDMPCSOLID6DOF_H_
#define _RIGIDMPCSOLID6DOF_H_

#include <Element.d/SuperElement.h>

class RigidMpcSolid6Dof : public SuperElement
{
  public:
    RigidMpcSolid6Dof(int, int*);
    int getTopNumber();
    int numTopNodes();
};

#endif

