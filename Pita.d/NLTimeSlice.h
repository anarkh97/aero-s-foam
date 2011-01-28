#ifndef _NLTIMESLICE_H_
#define _NLTIMESLICE_H_

#include <Math.d/Vector.h>
#include <Pita.d/DynamState.h>
#include <Pita.d/DynamStateSet.h>
#include <Pita.d/RankDeficientProjection.h>

class PitaNonLinDynamic;

struct NLTimeSlice
{
  typedef GenVector<double> VecType;
  typedef DynamState<double> State;
  typedef DynamStateSet<double> StateSet;

  NLTimeSlice() : sliceRank(-1), initialTime(0.0), finalTime(0.0), converged(false), active(false), projector(0, defaultTolerance) {}
  NLTimeSlice(const PitaNonLinDynamic & probDesc, int rank);
  void initialize(const PitaNonLinDynamic & probDesc, int rank);
  void clearAllBases();

  int sliceRank;
  double initialTime;
  double finalTime;
  bool converged;
  bool active;
  
  State seedState;
  State propState;
  State jumpState;
  State nextSeedState;
  State oldSeedState;
  RankDeficientProjection projector;
  StateSet propBase;
  StateSet localBase;
  
  // Data necessary to rebuild K when used as metric for orthogonalization
  VecType locDispOG;
  double locTimeOG;

  static const double defaultTolerance = 1.0e-6;
};

#endif
