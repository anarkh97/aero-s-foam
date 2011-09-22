#include "DistrVecNodeDof6Conversion.h"

namespace Rom {

DistrVecNodeDof6Conversion::~DistrVecNodeDof6Conversion() {
  typedef ConversionContainer::const_iterator ConversionIt;
  const ConversionIt it_end = subConversions_.end();
  for (ConversionIt it = subConversions_.begin(); it != it_end; ++it) {
    delete *it;
  }

  typedef RestrictedConversionContainer::const_iterator RestrictedConversionIt;
  const RestrictedConversionIt jt_end = subRestrictedConversions_.end();
  for (RestrictedConversionIt jt = subRestrictedConversions_.begin(); jt != jt_end; ++jt) {
    delete *jt;
  } 
}

} // end namespace Rom
