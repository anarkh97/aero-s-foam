#include <Pita.d/NLTimeSlice.h>
#include <Pita.d/PitaNonLinDynam.h>

void NLTimeSlice::initialize(const PitaNonLinDynamic & probDesc, int rank)
{
  //int Jratio = probDesc.getJratio();
  double sliceSpan = probDesc.getCoarseDt();
  int vectorSize = const_cast<PitaNonLinDynamic &>(probDesc).solVecInfo();
  
  sliceRank = rank;
  initialTime = rank * sliceSpan;
  finalTime = initialTime + sliceSpan;
  converged = false;
  active = false;
  seedBase.reset(vectorSize, 0);
  propBase.reset(vectorSize, 0);
  orthoBase.reset(vectorSize, 0);
  localBase.reset(vectorSize, probDesc.getJratio() + 1);
  seedState.reset(vectorSize);
  propState.reset(vectorSize);
  jumpState.reset(vectorSize);
  nextSeedState.reset(vectorSize, 0.0);
  oldSeedState.reset(vectorSize);
}

void NLTimeSlice::clearAllBases()
{
  seedBase.clear();
  orthoBase.clear();
  propBase.clear();
  localBase.clear();
}

