#ifndef _LINEARCONSTRAINTTYPE1_H_
#define _LINEARCONSTRAINTTYPE1_H_

#include <Element.d/MpcElement.d/ConstraintElement.h>

class LinearConstraintType1 : public ConstraintElement
{
  public:
    LinearConstraintType1(int*, DofSet*);
    void computeMPCs(CoordSet&);
    int getTopNumber();
    void update(GeomState& gState, CoordSet& cs);
};

#endif
