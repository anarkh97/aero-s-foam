#include "BasisId.h"

std::string
toString(BasisId::Type t) {
  static const std::string str[]  = { "state", "res", "jac" };
  return str[t];
}

std::string
toString(BasisId::Level l) {
  static const std::string str[] = { "snap",  "pod", "gappy" };
  return str[l];
}
