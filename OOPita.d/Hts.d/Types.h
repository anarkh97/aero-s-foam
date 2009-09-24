#ifndef PITA_HTS_TYPES_H
#define PITA_HTS_TYPES_H

#include "../Types.h"

namespace Pita { namespace Hts {

// Forward declarations
class PrimalSliceRank;
class DualSliceRank;

class HalfTimeSlice;
typedef Fwk::Numeric<HalfTimeSlice, int> HalfSliceCount;

class HalfSliceRank : public Fwk::Interval<HalfTimeSlice, int, HalfSliceCount> {
public:
  explicit HalfSliceRank(int v = 0) : Fwk::Interval<HalfTimeSlice, int, HalfSliceCount>(v) {}

  PrimalSliceRank primalRank() const;
  DualSliceRank dualRank() const;
};

inline
HalfSliceRank
operator+(const HalfSliceRank & a, const HalfSliceCount & d) {
  return HalfSliceRank(a.value() + d.value());
}

inline 
HalfSliceRank
operator-(const HalfSliceRank & a, const HalfSliceCount & d) {
  return HalfSliceRank(a.value() - d.value());
}
  
enum SeedType {
  UNDEFINED_SEED = 0,
  MAIN_SEED,
  LEFT_SEED,
  RIGHT_SEED,
  SEED_JUMP,
  SEED_CORRECTION
};

enum ReductionFlag {
  FULL = 0,
  REDUCED
};

inline
OStream &
operator<<(OStream & out, SeedType t) {
  char c;
  switch (t) {
    case UNDEFINED_SEED:          c = 'U'; break;
    case MAIN_SEED:               c = 'M'; break;
    case LEFT_SEED:               c = 'L'; break;
    case RIGHT_SEED:              c = 'R'; break;
    case SEED_JUMP:               c = 'J'; break;
    case SEED_CORRECTION:         c = 'C'; break;
    default:                      c = '?'; break;
  }
  out << c;
  return out;
}

enum SliceType {
  UNDEFINED_SLICE = 0,
  FORWARD_HALF_SLICE,
  BACKWARD_HALF_SLICE,
  HEAD_FULL_SLICE,
  TAIL_FULL_SLICE
};

inline
OStream &
operator<<(OStream & out, SliceType t) {
  char c;
  switch (t) {
    case UNDEFINED_SLICE:
      c = 'U';
      break;
    case FORWARD_HALF_SLICE:
      c = 'F';
      break;
    case BACKWARD_HALF_SLICE:
      c = 'B';
      break;
    case HEAD_FULL_SLICE:
      c = 'H';
      break;
    case TAIL_FULL_SLICE:
      c = 'T';
      break;
    default:
      break;
  }
  out << c;
  return out;
}

template <typename T>
class GenId {
public:
  T type() const { return type_; }
  HalfSliceRank rank() const { return rank_;}

  GenId(T t, HalfSliceRank r) :
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
  HalfSliceRank rank_;
};

typedef GenId<SeedType> SeedId;
typedef GenId<SliceType> SliceId;

inline
OStream &
operator<<(OStream & out, const SeedId & id) {
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
OStream &
operator<<(OStream & out, const CommId & c) {
  out << c.seed().type() << c.seed().rank() << "->" << c.cpu();
  return out;
}

class FullTimeSlice;
typedef Fwk::Numeric<FullTimeSlice, int> FullSliceCount;

class FullSliceRank : public Fwk::Interval<FullTimeSlice, int, FullSliceCount> {
public:
  explicit FullSliceRank(int v = 0) : Fwk::Interval<FullTimeSlice, int, FullSliceCount>(v) {}

  operator SliceRank() const { return SliceRank(value()); }
};

inline
FullSliceRank
operator+(const FullSliceRank & a, const FullSliceCount & d) {
  return FullSliceRank(a.value() + d.value());
}

inline 
FullSliceRank
operator-(const FullSliceRank & a, const FullSliceCount & d) {
  return FullSliceRank(a.value() - d.value());
}


class PrimalTimeSlice;
typedef Fwk::Numeric<PrimalTimeSlice, int> PrimalSliceCount;

class PrimalSliceRank : public Fwk::Interval<PrimalTimeSlice, int, PrimalSliceCount> {
public:
  explicit PrimalSliceRank(int v = 0) : Fwk::Interval<PrimalTimeSlice, int, PrimalSliceCount>(v) {}

  HalfSliceRank headRank() const; 
  HalfSliceRank tailRank() const;
};


class DualTimeSlice;
typedef Fwk::Numeric<DualTimeSlice, int> DualSliceCount;

class DualSliceRank : public Fwk::Interval<DualTimeSlice, int, DualSliceCount> {
public:
  explicit DualSliceRank(int v = 0) : Fwk::Interval<DualTimeSlice, int, DualSliceCount>(v) {}

  HalfSliceRank headRank() const; 
  HalfSliceRank tailRank() const;
};


inline
PrimalSliceRank
HalfSliceRank::primalRank() const {
  return PrimalSliceRank(value() >= 0 ? value() / 2 : -((1 - value()) / 2));
}

inline
DualSliceRank
HalfSliceRank::dualRank() const {
  return DualSliceRank(value() >= 1 ? (value() - 1) / 2 : -((-value()) / 2) - 1);
}


inline
HalfSliceRank
PrimalSliceRank::headRank() const {
  return HalfSliceRank(2 * value());
}


inline
HalfSliceRank
PrimalSliceRank::tailRank() const {
  return HalfSliceRank(2 * value() + 1);
}


inline
HalfSliceRank
DualSliceRank::headRank() const {
  return HalfSliceRank(2 * value() + 1);
}


inline
HalfSliceRank
DualSliceRank::tailRank() const {
  return HalfSliceRank(2 * value() + 2);
}

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_TYPES_H */
