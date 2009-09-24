#ifndef PITA_HTS_HALFTIMESLICE_H
#define PITA_HTS_HALFTIMESLICE_H

#include "Fwk.h"
#include "Types.h"
#include "../Seed.h"
#include "../Activity.h"
#include "../NamedTask.h"

namespace Pita { namespace Hts {

// Unique identifier for HalfTimeSlice
class HalfSliceId;

class HalfTimeSlice : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(HalfTimeSlice);

  enum Direction {
    NO_DIRECTION = 0,
    FORWARD,
    BACKWARD
  };
  
  typedef Fwk::GenManagerInterface<HalfTimeSlice *, HalfSliceId> Manager;

  // Read accessors
  HalfSliceRank rank() const { return rank_; }
  Direction direction() const { return direction_; }
  //PhaseRank phase() const { return phase_; }
  const Seed * seed() const { return seed_.ptr(); }
  Seed * propagatedSeed() const { return propagatedSeed_.ptr(); }

  // Mutators
  //virtual void phaseIs(PhaseRank p) { setPhase(p); }
  virtual void seedIs(const Seed * s) { setSeed(s); }
  virtual void propagatedSeedIs(Seed * ps) { setPropagatedSeed(ps); }

protected:
  HalfTimeSlice(HalfSliceRank rank, Direction direction);

  void setSeed(const Seed * s) { seed_ = s; }
  void setPropagatedSeed(Seed * ps) { propagatedSeed_ = ps; }
  //void setPhase(PhaseRank p) { phase_ = p; }

private:
  HalfSliceRank rank_;
  Direction direction_;
  //PhaseRank phase_;

  Seed::PtrConst seed_;
  Seed::Ptr propagatedSeed_;

  DISALLOW_COPY_AND_ASSIGN(HalfTimeSlice);
};

class HalfSliceId {
public:
  HalfSliceId(HalfSliceRank r = HalfSliceRank(-1), HalfTimeSlice::Direction d = HalfTimeSlice::NO_DIRECTION) :
    rank_(r),
    direction_(d)
  {}

  HalfSliceRank rank() const { return rank_; }
  HalfTimeSlice::Direction direction() const { return direction_; }

  bool operator==(const HalfSliceId & other) const {
    return (rank() == other.rank()) && (direction() == other.direction());
  }

  bool operator<(const HalfSliceId & other) const {
    return (rank() < other.rank()) ? true : (rank() == other.rank() && direction() < other.direction());
  }

private:
  HalfSliceRank rank_;
  HalfTimeSlice::Direction direction_;
};

OStream & operator<<(OStream & out, const HalfTimeSlice::Direction & d);
OStream & operator<<(OStream & out, const HalfSliceId & id);

struct HalfSliceComparator :
  public std::binary_function<const HalfTimeSlice::Ptr, const HalfTimeSlice::Ptr, bool> {
public:
  bool operator()(const HalfTimeSlice::Ptr & a, const HalfTimeSlice::Ptr & b) {
    return (HalfSliceId(a->rank(), a->direction()) < HalfSliceId(b->rank(), b->direction()));
  }
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_HALFTIMESLICE_H */