#ifndef ROM_BASISID_H
#define ROM_BASISID_H

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

std::string toString(BasisId::Type);
std::string toString(BasisId::Level);

#endif /* ROM_BASISID_H */
