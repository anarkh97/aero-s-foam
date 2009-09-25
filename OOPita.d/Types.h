#ifndef PITA_TYPES_H
#define PITA_TYPES_H

#include "Fwk.d/Fwk.h"

namespace Pita {

class Time;
typedef Fwk::Numeric<Time, double> Seconds;

class TimeStep;
typedef Fwk::Numeric<TimeStep, int> TimeStepCount;

class TimeSlice;
typedef Fwk::Numeric<TimeSlice, int> SliceCount;

class TimeSliceRankFix;
typedef Fwk::Numeric<TimeSliceRankFix, int> SliceRank;

class Cpu;
typedef Fwk::Numeric<Cpu, int> CpuCount;
class CpuTempFix;
typedef Fwk::Numeric<CpuTempFix, int> CpuRank;

class IterationRank : public Fwk::Ordinal<IterationRank, Fwk::U32> {
public:
  explicit IterationRank(Fwk::U32 rank = 0u) : Fwk::Ordinal<IterationRank, Fwk::U32>(rank) {}
  IterationRank next() const { return IterationRank(value() + 1); }
};

class PhaseRank : public Fwk::Ordinal<PhaseRank, U32> {
public:
  explicit PhaseRank(U32 rank = 0u) : Fwk::Ordinal<PhaseRank, U32>(rank) {}

  static PhaseRank none() { return PhaseRank(0); }
  static PhaseRank correction() { return PhaseRank(1); }
  static PhaseRank convergence() { return PhaseRank(2); }
  static PhaseRank basisUpdate() { return PhaseRank(3); }
  static PhaseRank globalCommunication() { return PhaseRank(4); }
  static PhaseRank fineGrid() { return PhaseRank(5); }
};

} /* end namespace Pita */

#endif /* PITA_TYPES_H */
