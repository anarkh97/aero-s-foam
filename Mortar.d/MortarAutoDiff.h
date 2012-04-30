#ifndef _MORTAR_DEFINES_H_
#define _MORTAR_DEFINES_H_

#ifdef USE_EIGEN3
#  include <Eigen/Core>
#  include <Eigen/Sparse>
#  include <unsupported/Eigen/AutoDiff>
#endif

#ifdef USE_SACADO
#  include <Sacado.hpp>
#endif

// see also definition of MAX_FFI_DERIVATIVES in Acme.d/search/Contact_Defines.h
#if defined(USE_EIGEN3) && defined(EIGEN_AUTODIFF_SCALAR_PJSA_HACK)
#  define MAX_MORTAR_DERIVATIVES 24
#  define MORTAR_AUTO_DIFF_EIGEN3
#elif defined(USE_EIGEN3) && defined(USE_SACADO)
#  define MAX_MORTAR_DERIVATIVES 24
#  define MORTAR_AUTO_DIFF_SACADO
#else
#  define MAX_MORTAR_DERIVATIVES 0
#endif

#if defined(MORTAR_AUTO_DIFF_EIGEN3)
#  ifndef COMPUTE_MORTAR_SECOND_DERIVATIVES
#    define MORTAR_AUTO_DIFF_EIGEN_FAD
     typedef Eigen::Matrix<double, MAX_MORTAR_DERIVATIVES, 1> DerivativeType;
     typedef Eigen::AutoDiffScalar<DerivativeType> ActiveDouble;
#  else
#    define MORTAR_AUTO_DIFF_EIGEN_FAD_FAD
     //for the computation of second derivatives using eigen... (seems to be working but it's slow)
     //TODO: only compute the second derivative if the lagrange multiplier is non-zero.
     typedef Eigen::Matrix<double, MAX_MORTAR_DERIVATIVES, 1> DerivativeType;
     typedef Eigen::Matrix<Eigen::AutoDiffScalar<DerivativeType>, MAX_MORTAR_DERIVATIVES,1> DerivativeType2;
     typedef Eigen::AutoDiffScalar<DerivativeType2> ActiveDouble;
#  endif
#elif defined(MORTAR_AUTO_DIFF_SACADO)
#  ifndef COMPUTE_MORTAR_SECOND_DERIVATIVES
#    define MORTAR_AUTO_DIFF_SACADO_FAD
     typedef Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES> ActiveDouble;
#  else
#    define MORTAR_AUTO_DIFF_SACADO_RAD_FAD
     typedef Sacado::RadVec::ADvar<Sacado::Fad::SFad<double,MAX_MORTAR_DERIVATIVES> > ActiveDouble;
#  endif
#endif

#if (MAX_MORTAR_DERIVATIVES > 0)
   struct MadNode {
     ActiveDouble x,y,z;
     MadNode(ActiveDouble _x, ActiveDouble _y, ActiveDouble _z) : x(_x), y(_y), z(_z) {}
   }; 
#  include<map>
   typedef std::map<int,MadNode*> MadCoordSet;
#endif

#endif
