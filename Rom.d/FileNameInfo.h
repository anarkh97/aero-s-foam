#ifndef ROM_FILENAMEINFO_H
#define ROM_FILENAMEINFO_H

#include <string>

class BasisId {
public:
  enum Type  { STATE,     RESIDUAL, JACOBIAN  };
  enum Level { SNAPSHOTS, POD,      GAPPY_POD };

  Type  type()  const { return type_; }
  Level level() const { return level_; }

  BasisId(Type type, Level level) :
    type_(type), level_(level)
  {}

private:
  Type type_;
  Level level_;
};

class FileNameInfo {
public:
  std::string fileName(const BasisId &) const;

  const std::string &prefix() const { return prefix_; }
  void prefixIs(const std::string &);

  FileNameInfo(); // Uses check filename in global geoSource to set prefix 
  explicit FileNameInfo(const std::string &prefix);

private:
  static std::string getPrefix(const std::string &);

  std::string prefix_;
};


#endif /* ROM_FILENAMEINFO_H */
