#include "HomogeneousTaskManager.h"

namespace Pita { namespace Hts {

class InitialSeed : public NamedTask {
public:
  void iterationIs(IterationRank i) {
      targetSeed_->stateIs(initializer_->initialSeed(seedRank_));
      targetSeed_->iterationIs(i);
      targetSeed_->statusIs((seedRank_ == SliceRank(0)) ? Seed::CONVERGED : Seed::ACTIVE);
    }

  InitialSeed(SeedInitializer * initializer, SliceRank seedRank, Seed * targetSeed) :
    NamedTask(buildName(seedRank)),
    initializer_(initializer),
    seedRank_(seedRank),
    targetSeed_(targetSeed)
  {}

protected:
  static String buildName(SliceRank seedRank) {
    return String("Initial Seed ") + toString(seedRank);
  }

private:
  SeedInitializer::Ptr initializer_;
  SliceRank seedRank_;
  Seed::Ptr targetSeed_;
};


HomogeneousTaskManager::HomogeneousTaskManager(LinearLocalNetwork * network,
                                               SeedInitializer * initializer,
                                               LinearProjectionNetworkImpl * correctionMgr,
                                               RemoteState::MpiManager * commMgr) :
  LinearTaskManager(IterationRank(0), network, correctionMgr, commMgr),
  initializer_(initializer)
{
  scheduleIterationZero();
}

void
HomogeneousTaskManager::iterationInc() {
  phases().clear();
  
  IterationRank nextIteration = iteration().next();
  
  correctionMgr()->prepareProjection();
  network()->convergedSlicesInc();
  scheduleNormalIteration();
  
  setIteration(nextIteration);
}

void
HomogeneousTaskManager::scheduleIterationZero() {
  scheduleSeedInitialization();
  scheduleFinePropagation();
}


void
HomogeneousTaskManager::scheduleSeedInitialization() {
  TaskList taskList;
  
  LinearLocalNetwork::MainSeedMap mainSeeds = network()->activeMainSeeds();
  for (LinearLocalNetwork::MainSeedMap::iterator it = mainSeeds.begin();
      it != mainSeeds.end();
      ++it) {
    taskList.push_back(new InitialSeed(initializer_.ptr(), it->first, it->second.ptr()));
  }
  
  Phase::Ptr finePropagation = phaseNew("Seed initialization", taskList);
  phases().push_back(finePropagation);
}

} /* end namespace Hts */ } /* end namespace Pita */
