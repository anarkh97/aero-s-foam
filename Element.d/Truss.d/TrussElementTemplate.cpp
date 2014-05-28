#ifdef USE_EIGEN3
#ifndef _TRUSSELEMENTSEMITEMPLATE_CPP_
#define _TRUSSELEMENTSEMITEMPLATE_CPP_

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <Element.d/Truss.d/TrussElementTemplate.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

template<typename doublereal>
void
TrussElementTemplate<doublereal>
::stiffnessMatrix(doublereal *_stiffness, doublereal *_x,doublereal *_y, doublereal *_z,
                  doublereal E, doublereal A, doublereal preload)
{
       using std::sqrt;

       Eigen::Map<Eigen::Matrix<doublereal,2,1> > x(_x), y(_y), z(_z);
       Eigen::Map<Eigen::Matrix<doublereal,6,6> > stiffness(_stiffness);

       doublereal dx = x[1] - x[0];
       doublereal dy = y[1] - y[0];
       doublereal dz = z[1] - z[0];

       doublereal length = sqrt( dx*dx + dy*dy + dz*dz );
       doublereal c1[3];

       c1[0] = dx/length;
       c1[1] = dy/length;
       c1[2] = dz/length;

       doublereal elementK = E*A/length;
       int i,j;
       for(i=0; i < 3; ++i) {
         for(j=0; j < 3; ++j) {
            stiffness(i,j)     = elementK*c1[i]*c1[j];
            stiffness(i+3,j+3) = elementK*c1[i]*c1[j];
            stiffness(i+3,j)   = -stiffness(i,j);
            stiffness(i,j+3)   = -stiffness(i,j);
         }
       }

       if(preload != 0.0) {
         for(i=0; i < 3; ++i) {
            stiffness(i,i)     += preload/length;
            stiffness(i+3,i+3) += preload/length;
            stiffness(i+3,i)   = -stiffness(i,i);
            stiffness(i,i+3)   = -stiffness(i,i);
          }
       }
}
#endif
#endif
