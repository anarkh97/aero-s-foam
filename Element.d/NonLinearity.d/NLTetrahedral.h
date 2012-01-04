#ifndef _NLTETRAHEDRAL_H_
#define _NLTETRAHEDRAL_H_

#include <Element.d/NonLinearity.d/SolidElementTemplate.h>

typedef SolidElementTemplate<Tetrahedron,4 ,1 > NLTetrahedral4;
typedef SolidElementTemplate<Tetrahedron,10,15> NLTetrahedral10;

#endif
