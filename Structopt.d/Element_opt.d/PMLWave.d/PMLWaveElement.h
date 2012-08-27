#ifndef __PML_WAVE_ELEMENT_H__
#ifdef STRUCTOPT

#include <Element.d/Element.h>

inline DComplex eps_x(const StructProp* sp)
{
  return DComplex(sp->P, sp->Ta);
}

inline DComplex eps_y(const StructProp* sp)
{
  return DComplex(sp->Q, sp->W);
}

inline DComplex eps_z(const StructProp* sp)
{
  return DComplex(sp->Ixx, sp->Iyy);
}

#endif
#endif
