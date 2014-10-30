#ifndef _EFFMEMBRANETRIANGLE_HPP_
#define _EFFMEMBRANETRIANGLE_HPP_

#ifdef USE_EIGEN3
#include <Eigen/Core>

template<typename doublereal>
class EffMembraneTriangle
{
  public:
    static Eigen::Matrix<doublereal,9,3>
    L(doublereal x[3], doublereal y[3], doublereal alpha);

    static Eigen::Matrix<doublereal,3,9>
    Bd(doublereal x[3], doublereal y[3], doublereal betam, doublereal zeta[3]);
};
#endif

#endif
