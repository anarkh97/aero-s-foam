#include "BasisId.h"

namespace Rom {

std::string
toString(BasisId::Type t) {
  static const std::string str[]  = { "state", "res", "jac", "for" };
  return str[t];
}

std::string
toString(BasisId::Level l) {
  static const std::string str[] = { "snap", "pod", "gappy" };
  return str[l];
}

} /* end namespace Rom */
