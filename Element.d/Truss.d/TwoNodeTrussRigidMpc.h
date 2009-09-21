#ifndef _TWONODETRUSSRIGIDMPC_H_
#define _TWONODETRUSSRIGIDMPC_H_

#include <Element.d/MpcElement.d/RigidMpcElement.h>

class TwoNodeTrussRigidMpc : public RigidMpcElement
{
 public:
  TwoNodeTrussRigidMpc(int*);
  void computeMPCs(CoordSet&);
  int getTopNumber();
  bool isSafe();
};

#endif
