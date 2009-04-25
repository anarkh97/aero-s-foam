#include "HalfSliceNetwork.h"
#include "Log.h"

namespace Pita {

HalfSliceNetwork::SeedInitializationReactor::SeedInitializationReactor(Activity * notifier, HalfSliceNetwork * parent) :
  Activity::Notifiee(notifier),
  parent_(parent)
{}

void
HalfSliceNetwork::SeedInitializationReactor::onStatus() {
  if (notifier()->status() == Activity::executing) {
    std::map<SliceRank, Seed *> & dualInitialSeed = parent()->dualInitialSeed_;
    if (parent()->dualSeedInitializer_) {
      for (std::map<SliceRank, Seed *>::iterator it = dualInitialSeed.begin(); it != dualInitialSeed.end(); ++it) {
        Seed * seed = it->second;
        seed->stateIs(parent()->dualSeedInitializer_->initialSeed(it->first));
        seed->statusIs(Seed::SPECIAL);
      }
    }
    dualInitialSeed.clear();
    
    std::map<SliceRank, Seed *> & primalInitialSeed = parent()->primalInitialSeed_;
    for (std::map<SliceRank, Seed *>::iterator it = primalInitialSeed.begin(); it != primalInitialSeed.end(); ++it) {
      Seed * seed = it->second;
      seed->stateIs(parent()->primalSeedInitializer_->initialSeed(it->first));
      seed->statusIs(it->first != SliceRank(0) ? Seed::ACTIVE : Seed::CONVERGED);
    }
    primalInitialSeed.clear();
  }
}

HalfSliceNetwork::NewIterationReactor::NewIterationReactor(Activity * notifier, HalfSliceNetwork * parent) :
  Activity::Notifiee(notifier),
  parent_(parent)
{}

void
HalfSliceNetwork::NewIterationReactor::onStatus() {
  switch (notifier()->status()) {
    case Activity::executing:
      //parent()->recvMgr_->iterationIs(notifier()->currentIteration());
      parent()->surrogate_->iterationIs(notifier()->currentIteration()); // HACK TODO
      break;
    case Activity::free:
      parent()->mapping_->convergedSlicesInc(); // HACK
      notifier()->iterationIs(notifier()->currentIteration().next());
      notifier()->statusIs(Activity::scheduled);
      break;
    default:
      break;
  }
}

HalfSliceNetwork::HalfSliceNetwork(size_t vectorSize,
                                   SliceMapping * mapping,
                                   CpuRank localCpu,
                                   Seed::Manager * seedMgr,
                                   HalfTimeSlice::Manager * hsMgr,
                                   FullTimeSliceHead::Manager * fshMgr,
                                   FullTimeSliceTail::Manager * fstMgr,
                                   ScheduledRemoteSeedReader::Manager<CommId> * sendMgr,
                                   ScheduledRemoteSeedWriter::Manager * recvMgr,
                                   RemoteHalfSliceSurrogate * surrogate,// HACK
                                   HalfSliceSchedule * schedule,
                                   SeedInitializer * primalSeedInit,
                                   SeedInitializer * dualSeedInit) : 
  vectorSize_(vectorSize),
  mapping_(mapping),
  localCpu_(localCpu),
  seedMgr_(seedMgr),
  hsMgr_(hsMgr),
  fshMgr_(fshMgr),
  fstMgr_(fstMgr),
  sendMgr_(sendMgr),
  recvMgr_(recvMgr),
  surrogate_(surrogate), // HACK
  schedule_(schedule),
  topology_(HalfSliceNetworkTopology::New()),
  primalSeedInitializer_(primalSeedInit),
  dualSeedInitializer_(dualSeedInit),
  seedInitializationReactor_(new SeedInitializationReactor(activityManagerInstance()->activityNew("Seed Initialization").ptr(), this)),
  newIterationReactor_(new NewIterationReactor(activityManagerInstance()->activityNew("New Iteration").ptr(), this))
{
  init();
}

void
HalfSliceNetwork::init() {
  HalfSliceRank firstActive = mapping_->firstActiveSlice();
  HalfSliceRank firstInactive = mapping_->firstInactiveSlice();
  for (SliceMapping::SliceIdIterator sit = mapping_->hostedSlice(localCpu_, firstActive, firstInactive) ; sit; ++sit) {

    // Create the added slice
    HalfSliceRank rank = (*sit).rank();
    switch ((*sit).type()) {
      case FORWARD_HALF_SLICE:
        {
          HalfSliceId hsId(rank, HalfTimeSlice::FORWARD);
          PhaseRank slicePhase = schedule_->localPropagation(rank);
          //log() << "HalfTimeSliceNew(" << hsId << ") at phase " << slicePhase << "\n";
          HalfTimeSlice::Ptr slice = hsMgr_->instanceNew(hsId);
          slice->phaseIs(slicePhase);
          Seed * seed = getSeed(SeedId(MAIN_SEED, rank));
          if (rank.value() % 2 == 0) { // Add to Seed initialization list if necessary
            //log() << SeedId(MAIN_SEED, rank) << " to be initialized\n";
            primalInitialSeed_[SliceRank((rank.value() / 2))] = seed; 
          } else {
            dualInitialSeed_[SliceRank(((rank.value() + 1) / 2))] = seed;
          } 
          slice->seedIs(seed);
          slice->propagatedSeedIs(getSeed(SeedId(LEFT_SEED, rank + HalfSliceCount(1))));
          //collector_->lastSourceIs(slice.ptr());
        }
        break;
      case BACKWARD_HALF_SLICE:
        {
          HalfSliceId hsId(rank, HalfTimeSlice::BACKWARD);
          PhaseRank slicePhase = schedule_->localPropagation(rank);
          //log() << "HalfTimeSliceNew(" << hsId << ") at phase " << slicePhase << "\n";
          HalfTimeSlice::Ptr slice = hsMgr_->instanceNew(hsId);
          slice->phaseIs(slicePhase);
          HalfSliceRank seedRank(rank + HalfSliceCount(1));
          Seed * seed = getSeed(SeedId(MAIN_SEED, seedRank));
          if (seedRank.value() % 2 == 0) { // Add to Seed initialization list if necessary
            //log() << SeedId(MAIN_SEED, seedRank) << " to be initialized\n";
            primalInitialSeed_[SliceRank((seedRank.value() / 2))] = seed; 
          } else {
            dualInitialSeed_[SliceRank(((rank.value() + 1) / 2))] = seed;
          } 
          slice->seedIs(seed);
          slice->propagatedSeedIs(getSeed(SeedId(RIGHT_SEED, rank)));
          //collector_->lastSourceIs(slice.ptr());
        }
        break;
      case HEAD_FULL_SLICE:
        {
          PhaseRank slicePhase = schedule_->correction(rank);
          //log() << "FullTimeSliceHeadNew(" << rank << ") at phase " << slicePhase << "\n";
          CpuRank tailCpu = mapping_->hostCpu(SliceId(TAIL_FULL_SLICE, rank + HalfSliceCount(1)));

          FullTimeSliceHead::Ptr slice = fshMgr_->instanceNew(rank);

          slice->tailCpuIs(tailCpu); 
          slice->phaseIs(slicePhase);
          slice->updatedSeedIs(getSeed(SeedId(MAIN_SEED, rank)));
          slice->rightPropagatedSeedIs(getSeed(SeedId(RIGHT_SEED, rank)));
        }
        break;
      case TAIL_FULL_SLICE:
        {
          //log() << "FullTimeSliceTailNew(" << rank << ")\n";
          PhaseRank slicePhase = schedule_->correction(rank);
          CpuRank headCpu = mapping_->hostCpu(SliceId(HEAD_FULL_SLICE, rank - HalfSliceCount(1)));

          FullTimeSliceTail::Ptr slice = fstMgr_->instanceNew(rank);
          //log() << "HeadCpu set to " << headCpu << "\n";

          slice->headCpuIs(headCpu);
          slice->phaseIs(slicePhase);
          slice->nextUpdatedSeedIs(getSeed(SeedId(MAIN_SEED, rank + HalfSliceCount(1))));
          slice->nextLeftPropagatedSeedIs(getSeed(SeedId(LEFT_SEED, rank + HalfSliceCount(1))));
        }
        break;
      default:
        throw Fwk::InternalException();
    } 
  }

  seedInitializationReactor_->notifier()->phaseIs(schedule_->mainSeedSynchronization());
  seedInitializationReactor_->notifier()->iterationIs(IterationRank(0));
  seedInitializationReactor_->notifier()->statusIs(Activity::scheduled);
  
  newIterationReactor_->notifier()->phaseIs(PhaseRank(0));
  newIterationReactor_->notifier()->iterationIs(IterationRank(1));
  newIterationReactor_->notifier()->statusIs(Activity::scheduled);

  //collector_->consolidationPhaseIs(schedule_->baseBuilding());
}


Seed *
HalfSliceNetwork::getSeed(const SeedId & id) {
  OStringStream os;
  os << id;
  String seedName = os.str();
  //log() << "Seed(" << seedName << ")\n";
  Seed::Ptr currSeed = seedMgr_->instance(seedName);
  if (!currSeed) {
    currSeed = seedMgr_->instanceNew(seedName);
    //log() << "SeedNew(" << seedName << ")\n";

    // Setup communication framework to account for created seed
    SliceId writerId = *topology_->writer(id); // Assume one and only one writer
    CpuRank writerCpu = mapping_->hostCpu(writerId);
    
    if (writerCpu != localCpu_ && writerCpu != CpuRank(-1)) { // Non local maintainer => Set up RemoteWriter
      //log() << "Maintainer " << writerId.rank() << writerId.type() << " of " << seedName << " is on remote cpu " << writerCpu << "\n";
      String writerName = seedName + "_" + toString(writerCpu);

      // Loop over local consumers to determine scheduling
      for (HalfSliceNetworkTopology::SliceIteratorConst reader_it = topology_->reader(id); reader_it; ++reader_it) {
        SliceId readerId = *reader_it;
        CpuRank readerCpu = mapping_->hostCpu(readerId);
        if (readerCpu == localCpu_) {
          //log() << "Local Reader for " << seedName << " is " << readerId.rank() << readerId.type() << "\n";

          PhaseRank newPhase = getSeedSyncPhase(writerId.type(), readerId.type(), writerId.rank());
          //log() << "RemoteWriter(" << writerName << ") at phase " << newPhase << "\n";

          ScheduledRemoteSeedWriter::Ptr remoteWriter = recvMgr_->instance(id);
          if (!remoteWriter) { // RemoteWriter does not exist, create it
            //log() << "RemoteWriterNew(" << writerName << ")\n";
            ScheduledRemoteSeedWriter::Ptr remoteWriter = recvMgr_->instanceNew(id);
            remoteWriter->originCpuIs(writerCpu);
            remoteWriter->vectorSizeIs(vectorSize_);
            remoteWriter->targetIs(currSeed.ptr());
            remoteWriter->phaseIs(newPhase);
          } else { // RemoteWriter exists, update the scheduling if necessary
            if (newPhase < remoteWriter->phase()) {
              remoteWriter->phaseIs(newPhase);
            }
          }
        }
      }
    } else { // Local maintainer => Set up RemoteReader
      //log() << "Maintainer " << writerId.rank() << writerId.type() << " of " << seedName << " is on local cpu\n";
      for (HalfSliceNetworkTopology::SliceIteratorConst reader_it = topology_->reader(id); reader_it; ++reader_it) {
        SliceId readerId = *reader_it;
        CpuRank readerCpu = mapping_->hostCpu(readerId);
        //log() << "Interested slice is " << readerId.type() << readerId.rank() << " on cpu " << readerCpu << "\n";
        if (readerCpu != localCpu_ && readerCpu != CpuRank(-1)) {
          String readerName = seedName + "_" + toString(readerCpu);
          PhaseRank newPhase = getSeedSyncPhase(writerId.type(), readerId.type(), readerId.rank());
          //log() << "RemoteReader(" << readerName << ") at phase " << newPhase << "\n";
          CommId remoteReaderId(id, readerCpu);
          ScheduledRemoteSeedReader::Ptr remoteReader = sendMgr_->instance(remoteReaderId);
          if (!remoteReader) { // RemoteReader does not exist, create it
            //log() << "RemoteReaderNew(" << readerName << ")\n";
            remoteReader = sendMgr_->instanceNew(remoteReaderId);
            remoteReader->targetCpuIs(readerCpu);
            remoteReader->notifierIs(currSeed.ptr());
            remoteReader->phaseIs(newPhase);
          } else { // RemoteReader exists, update scheduling if necessary
            if (newPhase < remoteReader->phase()) {
              //log() << "Changing phase of " << readerName << " from " << remoteReader->phase() << " to " << newPhase << "\n";
              remoteReader->phaseIs(newPhase);
            }
          }
        }
      }
    }
  }
  return currSeed.ptr();
}

PhaseRank
HalfSliceNetwork::getSeedSyncPhase(SliceType originSliceType, SliceType targetSliceType, HalfSliceRank remoteSliceRank) const {
  if (targetSliceType == FORWARD_HALF_SLICE || targetSliceType == BACKWARD_HALF_SLICE) {
    return schedule_->mainSeedSynchronization(remoteSliceRank);
  }
  if (originSliceType == FORWARD_HALF_SLICE || originSliceType == BACKWARD_HALF_SLICE) {
    return schedule_->propagatedSeedSynchronization(remoteSliceRank);
  }
  return schedule_->correction(remoteSliceRank);
}

} // namespace Pita
