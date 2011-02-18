#include "FileNameInfo.h"

#include <Driver.d/GeoSource.h>

#include <sstream>

extern GeoSource *geoSource;

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

inline
std::string
FileNameInfo::getPrefix(const std::string &str) {
  return str.substr(0, str.find('.'));
}

FileNameInfo::FileNameInfo() :
  prefix_(getPrefix(geoSource->getCheckFileInfo()->checkfile))
{}

FileNameInfo::FileNameInfo(const std::string &prefix) :
  prefix_(prefix)
{}

std::string
FileNameInfo::fileName(const BasisId &id) const {
  std::ostringstream builder;

  const char suffix[] = "rob";
  builder << prefix()
          << "." << toString(id.type())
          << "." << toString(id.level())
          << "." << suffix;

  return builder.str();
}
