#ifndef _TRANSLATIONALSPRING_H_
#define _TRANSLATIONALSPRING_H_

#include <Element.d/Joint.d/ConstantDistanceConstraint.h>

class TranslationalSpring : public ConstantDistanceConstraint
{
  public:
    TranslationalSpring(int*);
    void setProp(StructProp *p, bool _myProp = false);
};

#endif
