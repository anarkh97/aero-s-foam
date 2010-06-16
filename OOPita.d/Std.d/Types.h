#ifndef PITA_STD_TYPES_H
#define PITA_STD_TYPES_H

#include "../Types.h"

namespace Pita { namespace Std {

class TimeSlice;
typedef Fwk::Numeric<TimeSlice, int> SliceCount;

class SliceRank : public Fwk::Interval<TimeSlice, int, SliceCount> {
public:
  explicit SliceRank(int v = 0) : Fwk::Interval<TimeSlice, int, SliceCount>(v) {}

  SliceRank next() const { return SliceRank(value() + 1); }
  SliceRank previous() const { return SliceRank(value() - 1); }
};

inline
SliceRank
operator+(const SliceRank & a, const SliceCount & d) {
  return SliceRank(a.value() + d.value());
}

inline 
SliceRank
operator-(const SliceRank & a, const SliceCount & d) {
  return SliceRank(a.value() - d.value());
}
  
enum SeedType {
  UNDEFINED_SEED = 0,
  MAIN_SEED,
  PROPAGATED_SEED,
  SEED_JUMP,
  SEED_CORRECTION
};

enum ReductionFlag {
  FULL = 0,
  REDUCED
};

inline
Fwk::OStream &
operator<<(Fwk::OStream & out, SeedType t) {
  char c;
  switch (t) {
    case UNDEFINED_SEED:          c = 'U'; break;
    case MAIN_SEED:               c = 'M'; break;
    case PROPAGATED_SEED:         c = 'P'; break;
    case SEED_JUMP:               c = 'J'; break;
    case SEED_CORRECTION:         c = 'C'; break;
    default:                      c = '?'; break;
  }
  out << c;
  return out;
}


template <typename T>
class GenId {
public:
  T type() const { return type_; }
  SliceRank rank() const { return rank_;}

  GenId(T t, SliceRank r) :
    type_(t), rank_(r)
  {}

  bool operator==(const GenId<T> & other) const {
    return (type() == other.type()) && (rank() == other.rank());
  }

  bool operator<(const GenId<T> & other) const {
  return (rank() < other.rank()) ? true : ((rank() == other.rank()) ? (type() < other.type()) : false);
}

private:
  T type_;
  SliceRank rank_;
};

typedef GenId<SeedType> SeedId;

inline
Fwk::OStream &
operator<<(Fwk::OStream & out, const SeedId & id) {
  out << id.type() << id.rank();
  return out;
}

class CommId {
public:
  SeedId seed() const { return seed_; }
  CpuRank cpu() const { return cpu_; }

  CommId(const SeedId & seed, CpuRank cpu) :
    seed_(seed),
    cpu_(cpu)
  {}

  bool operator==(const CommId & other) const {
    return (seed() == other.seed()) && (cpu() == other.cpu()); 
  }

private:
  SeedId seed_;
  CpuRank cpu_;
};

inline
bool
operator<(const CommId & a, const CommId & b) {
  return (a.seed() == b.seed()) ? (a.cpu() < b.cpu()) : (a.seed() < b.seed());
}

inline
Fwk::OStream &
operator<<(Fwk::OStream & out, const CommId & c) {
  out << c.seed().type() << c.seed().rank() << "->" << c.cpu();
  return out;
}

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_TYPES_H */
