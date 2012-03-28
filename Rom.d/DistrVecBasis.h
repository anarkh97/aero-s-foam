#ifndef ROM_DISTRVECBASIS_H
#define ROM_DISTRVECBASIS_H

#include "VecBasis.h"

#include <Feti.d/DistrVector.h>

#ifdef _TR1_FUNCTIONAL
  #include <tr1/functional>
#elif defined(USE_BOOST)
  #include <boost/version.hpp>
  #if BOOST_VERSION >= 103400
    #include <boost/tr1/functional.hpp>
  #endif
#endif

namespace Rom {

extern const DistrInfo DEFAULT_DISTR_INFO;

template <typename Scalar>
struct VecTraits<Scalar, GenDistrVector> {
  typedef GenDistrVector<Scalar> Type;
  typedef typename Type::InfoType InfoType;
  typedef typename std::tr1::reference_wrapper<const DistrInfo> InternalInfoType;
  
  static InfoType defaultInfo() { return DEFAULT_DISTR_INFO; }
  static int length(InfoType info) { return info.len; }
  static bool equals(InfoType i, InfoType j) { return &i == &j; }
  static bool not_equals(InfoType i, InfoType j) { return &i != &j; }
};

typedef GenVecBasis<double, GenDistrVector> DistrVecBasis;

} // end namespace Rom

#endif /* ROM_DISTRVECBASIS_H */
