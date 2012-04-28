#ifndef _NLPENTAHEDRAL_H_
#define _NLPENTAHEDRAL_H_

#include <Element.d/NonLinearity.d/SolidElementTemplate.h>

typedef SolidElementTemplate<Wedge,6, 6 > NLPentahedral6;
typedef SolidElementTemplate<Wedge,15,9 > NLPentahedral15;
typedef SolidElementTemplate<Wedge,26,18> NLPentahedral26;

#endif
