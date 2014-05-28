#ifndef _TRUSSELEMENTSEMITEMPLATE_HPP_
#define _TRUSSELEMENTSEMITEMPLATE_HPP_

#include <string>

template<typename doublereal>
class TrussElementTemplate 
{
#ifdef USE_EIGEN3
  public:
    void stiffnessMatrix(doublereal *_stiffness, doublereal *_x,doublereal *_y, doublereal *_z,
                         doublereal E, doublereal A, doublereal preload); 

  protected:
#endif
};

#endif
