#ifndef _RIGIDMPCBEAM_H_
#define _RIGIDMPCBEAM_H_

#include <Element.d/MpcElement.d/RigidMpcElement.h>

class RigidMpcBeam : public RigidMpcElement
{
 public:
  RigidMpcBeam(int*);
  void computeMPCs(CoordSet&);
  int getTopNumber();
};

#endif
