#include "FileNameInfo.h"

#include <Driver.d/GeoSource.h>

#include <sstream>

extern GeoSource *geoSource;

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
FileNameInfo::basisFileName(const BasisId &id) const {
  std::ostringstream builder;

  const char suffix[] = "rob";
  builder << prefix()
          << "." << toString(id.type())
          << "." << toString(id.level())
          << "." << suffix;

  return builder.str();
}
