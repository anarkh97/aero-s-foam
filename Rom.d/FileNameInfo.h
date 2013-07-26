#ifndef ROM_FILENAMEINFO_H
#define ROM_FILENAMEINFO_H

#include "BasisId.h"

#include <string>

namespace Rom {

class FileNameInfo {
public:
  std::string basisFileName(const BasisId &, int index = 0) const;

  const std::string &prefix() const { return prefix_; }
  void prefixIs(const std::string &);

  FileNameInfo(); // Uses check filename in global geoSource to set prefix 
  explicit FileNameInfo(const std::string &prefix);

private:
  static std::string getPrefix(const std::string &);

  std::string prefix_;
};

class BasisFileId {
public:
  BasisFileId(const FileNameInfo &info, const BasisId &id) :
    id_(id), name_(initName(info))
  {}

  BasisFileId(const FileNameInfo &info, BasisId::Type type, BasisId::Level level, int index = 0) :
    id_(type, level), name_(initName(info, index))
  {}

  BasisId id() const { return id_; }
  BasisId::Type type() const { return id_.type(); }
  BasisId::Level level() const { return id_.level(); }

  const std::string &name() const { return name_; }
  operator std::string() const { return name_; }

private:
  std::string initName(const FileNameInfo &info, int index = 0) const { return info.basisFileName(id_, index); }

  const BasisId id_;
  const std::string name_;
};

} /* end namespace Rom */

#endif /* ROM_FILENAMEINFO_H */
