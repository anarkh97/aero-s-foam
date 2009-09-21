#ifndef _RIGIDMPCSOLID_H_
#define _RIGIDMPCSOLID_H_

#include <Element.d/MpcElement.d/RigidMpcElement.h>

class RigidMpcSolid : public RigidMpcElement
{
 public:
  RigidMpcSolid(int, int*);
  void computeMPCs(CoordSet&);
  int getTopNumber();
  int numTopNodes();
};

#endif

