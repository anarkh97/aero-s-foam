#ifndef _RIGIDMPCSPRINGRO_H_
#define _RIGIDMPCSPRINGRO_H_

#include <Element.d/MpcElement.d/RigidMpcElement.h>

class RigidMpcSpringRo : public RigidMpcElement
{
 public:
  RigidMpcSpringRo(int*);
  void computeMPCs(CoordSet&);
  int getTopNumber();
  bool isSafe();
};

#endif
