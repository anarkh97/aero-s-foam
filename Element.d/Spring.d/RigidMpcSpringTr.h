#ifndef _RIGIDMPCSPRINGTR_H_
#define _RIGIDMPCSPRINGTR_H_

#include <Element.d/MpcElement.d/RigidMpcElement.h>

class RigidMpcSpringTr : public RigidMpcElement
{
 public:
  RigidMpcSpringTr(int*);
  void computeMPCs(CoordSet&);
  int getTopNumber();
  bool isSafe();
};

#endif
