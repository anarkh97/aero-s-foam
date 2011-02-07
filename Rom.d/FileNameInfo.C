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
FileNameInfo::fileName(const BasisId &id) const {
  static const std::string typeStr[]  = { "state", "res", "jac"   };
  static const std::string levelStr[] = { "snap",  "pod", "gappy" };

  std::ostringstream builder;
  builder << prefix()             << "."
          << typeStr[id.type()]   << "."
          << levelStr[id.level()] << "."
          << "rob";

  return builder.str();
}
