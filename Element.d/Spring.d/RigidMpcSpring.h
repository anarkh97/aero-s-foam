#ifndef _RIGIDMPCSPRING_H_
#define _RIGIDMPCSPRING_H_

#include <Element.d/MpcElement.d/RigidMpcElement.h>

class RigidMpcSpring : public RigidMpcElement
{
 public:
  RigidMpcSpring(int*);
  void computeMPCs(CoordSet&);
  int getTopNumber();
};

#endif
