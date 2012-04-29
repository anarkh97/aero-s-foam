#ifndef _MORTAR_DEFINES_H_
#define _MORTAR_DEFINES_H_

#ifdef USE_EIGEN3
#include <Eigen/Sparse>
#endif

// see also definition of MAX_FFI_DERIVATIVES in acme
#ifndef MAX_MORTAR_DERIVATIVES
#  define MAX_MORTAR_DERIVATIVES 0
#endif

#if (MAX_MORTAR_DERIVATIVES > 0)
#if defined(USE_SACADO) && defined(MORTAR_AUTO_DIFF_SACADO)
#  include "Sacado.hpp"
#  define MORTAR_AUTO_DIFF_SACADO_FAD
   typedef Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES> MadDouble;
   //for the computation of second derivatives using sacado... (not quite working)
   //#  define MORTAR_AUTO_DIFF_SACADO_RAD_FAD
   //typedef Sacado::RadVec::ADvar<Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES> > MadDouble;
#elif defined(USE_EIGEN3) && defined(MORTAR_AUTO_DIFF_EIGEN3)
#  include <Eigen/Core>
#  include <unsupported/Eigen/AutoDiff>
#  define MORTAR_AUTO_DIFF_EIGEN_FAD
   typedef Eigen::Matrix<double, MAX_MORTAR_DERIVATIVES, 1> DerivativeType;
   typedef Eigen::AutoDiffScalar<DerivativeType> MadDouble;
   //for the computation of second derivatives using eigen... (seems to be working but it's slow)
   //TODO: only compute the second derivative if the lagrange multiplier is non-zero.
//#  define MORTAR_AUTO_DIFF_EIGEN_FAD_FAD
//   typedef Eigen::Matrix<double, MAX_MORTAR_DERIVATIVES, 1> DerivativeType;
//   typedef Eigen::Matrix<Eigen::AutoDiffScalar<DerivativeType>, MAX_MORTAR_DERIVATIVES,1> DerivativeType2;
//   typedef Eigen::AutoDiffScalar<DerivativeType2> MadDouble;
#else
#  error "Computation of mortar constraint derivatives requires Sacado or Eigen3"
#endif
#endif

#if (MAX_MORTAR_DERIVATIVES > 0)
   struct MadNode {
     MadDouble x,y,z;
     MadNode(MadDouble _x, MadDouble _y, MadDouble _z) : x(_x), y(_y), z(_z) {}
   }; 
#  include<map>
   typedef std::map<int,MadNode*> MadCoordSet;
#endif

#endif
