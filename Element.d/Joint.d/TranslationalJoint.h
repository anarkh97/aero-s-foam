#ifndef _TRANSLATIONALJOINT_H_
#define _TRANSLATIONALJOINT_H_

#include <Element.d/SuperElement.h>

class TranslationalJoint : public SuperElement
{
  public:
    TranslationalJoint(int*);
    int getTopNumber();
};

#endif
