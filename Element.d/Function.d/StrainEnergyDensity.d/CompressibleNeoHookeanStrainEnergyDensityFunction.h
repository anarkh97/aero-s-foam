#ifndef _COMPRESSIBLENEOHOOKEANSTRAINENERGYDENSITYFUNCTION_H_
#define _COMPRESSIBLENEOHOOKEANSTRAINENERGYDENSITYFUNCTION_H_

#include <Element.d/Function.d/StrainEnergyDensity.d/StrainEnergyDensityFunction.h>

// reference: Bonet & Wood, "Nonlinear continuum mechanics for finite element analysis", 2nd edition (2008), p 162 equation 6.27

namespace Simo {

template<typename Scalar>
class CompressibleNeoHookeanStrainEnergyDensityFunction
: public StrainEnergyDensityFunction<Scalar>
{
    double lambda, mu; // Lame constants

  public:
    CompressibleNeoHookeanStrainEnergyDensityFunction(double _lambda, double _mu) : lambda(_lambda), mu(_mu) {}
    CompressibleNeoHookeanStrainEnergyDensityFunction(const Eigen::Array<double,2,1> &sconst, const Eigen::Array<int,0,1>&) {
      lambda = sconst[0];
      mu     = sconst[1];
    }

    Scalar operator() (const Eigen::Matrix<Scalar,3,3> &F) const {
      // inputs: deformation gradient
      // output: strain energy density

      using std::log;
      Scalar lnJ = log(F.determinant());
      return mu/2*((F*F.transpose()).trace() - 3) - mu*lnJ + lambda/2*lnJ*lnJ;
    }

};

} // namespace Simo

#endif
