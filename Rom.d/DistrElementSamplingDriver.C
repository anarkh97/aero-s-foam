#include "DistrElementSamplingDriver.h"
#include "SubElementSamplingDriver.h"

#include "DistrDomainUtils.h"

#include "DistrVecBasis.h"

#include "DistrMasterMapping.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>
#include <cassert>
#include <iostream>

namespace Rom {


const DistrInfo&
DistrElementSamplingDriver::vectorSize() const {
  return decDomain_->masterSolVecInfo();
}

DistrElementSamplingDriver::DistrElementSamplingDriver(Domain *domain, Communicator *comm) :
  domain_(domain),
  comm_(comm),
  decDomain_(createDecDomain<double>(domain))
{}

void
DistrElementSamplingDriver::solve() {
  decDomain_->preProcess();

  FileNameInfo fileInfo;

  // Read order reduction data
  DistrVecBasis podBasis;
  DistrBasisInputFile podBasisFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD));

  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, podBasisFile.stateCount()) :
                                     podBasisFile.stateCount();

  filePrint(stderr, " ... Projection subspace of dimension = %d ...\n", projectionSubspaceSize);
  podBasis.dimensionIs(projectionSubspaceSize, vectorSize());

  DistrVecNodeDof6Conversion converter(decDomain_->getAllSubDomains(), decDomain_->getAllSubDomains() + decDomain_->getNumSub());

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain_->getAllSubDomains()),
                                   SubDomIt(decDomain_->getAllSubDomains() + decDomain_->getNumSub()));
  DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  for (DistrVecBasis::iterator it = podBasis.begin(),
                               it_end = podBasis.end();
                               it != it_end; ++it) {
    assert(podBasisFile.validCurrentState());

    podBasisFile.currentStateBuffer(buffer);
    converter.vector(buffer, *it);

    podBasisFile.currentStateIndexInc();
  }

  // Read state snapshots
  // TODO fix skip and offset
  DistrVecBasis snapshots;
  std::vector<double> timeStamps;
  {
    DistrBasisInputFile in(BasisFileId(fileInfo, BasisId::STATE, BasisId::SNAPSHOTS));
    const int skipFactor = 1; //domain->solInfo().skipPodRom;
    const int skipOffSet = 0; //domain->solInfo().skipOffSet;
    const int basisStateCount = 1 + (in.stateCount() - 1) / skipFactor;

    snapshots.dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);

    for (DistrVecBasis::iterator it = snapshots.begin(),
                                 it_end = snapshots.end();
                                 it != it_end; ++it) {
      assert(in.validCurrentState());

      in.currentStateBuffer(buffer);
      converter.vector(buffer, *it);
      timeStamps.push_back(in.currentStateHeaderValue());

      in.currentStateIndexInc();
    }
  }

  // Read velocity snapshots
  // TODO fix skip and offset
  DistrVecBasis *velocSnapshots = 0;
  if(domain->solInfo().velocPodRomFile != "") {
    std::vector<double> timeStamps;
    velocSnapshots = new DistrVecBasis;
    DistrBasisInputFile in(BasisFileId(fileInfo, BasisId::VELOCITY, BasisId::SNAPSHOTS));
    const int skipFactor = 1; //domain->solInfo().skipPodRom;
    const int skipOffSet = 0; //domain->solInfo().skipOffSet;
    const int basisStateCount = 1 + (in.stateCount() - 1) / skipFactor;

    velocSnapshots->dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);

    for (DistrVecBasis::iterator it = velocSnapshots->begin(),
                                 it_end = velocSnapshots->end();
                                 it != it_end; ++it) {
      assert(in.validCurrentState());

      in.currentStateBuffer(buffer);
      converter.vector(buffer, *it);
      timeStamps.push_back(in.currentStateHeaderValue());

      in.currentStateIndexInc();
    }
  }

  // Read acceleration snapshots
  // TODO fix skip and offset
  DistrVecBasis *accelSnapshots = 0;
  if(domain->solInfo().accelPodRomFile != "") {
    std::vector<double> timeStamps;
    accelSnapshots = new DistrVecBasis;
    DistrBasisInputFile in(BasisFileId(fileInfo, BasisId::ACCELERATION, BasisId::SNAPSHOTS));
    const int skipFactor = 1; //domain->solInfo().skipPodRom;
    const int skipOffSet = 0; //domain->solInfo().skipOffSet;
    const int basisStateCount = 1 + (in.stateCount() - 1) / skipFactor;
  
    accelSnapshots->dimensionIs(basisStateCount, vectorSize());
    timeStamps.reserve(basisStateCount);
  
    for (DistrVecBasis::iterator it = accelSnapshots->begin(),
                                 it_end = accelSnapshots->end();
                                 it != it_end; ++it) {
      assert(in.validCurrentState());
  
      in.currentStateBuffer(buffer);
      converter.vector(buffer, *it);
      timeStamps.push_back(in.currentStateHeaderValue());

      in.currentStateIndexInc();
    }
  }

  const int podVectorCount = podBasis.vectorCount();
  const int snapshotCount = snapshots.vectorCount();

  // Temporary buffers shared by all iterations
  Vector podComponents(podVectorCount);

  // Project snapshots on POD basis to get training configurations
  DistrVecBasis displac(snapshotCount, vectorSize());
  for (int iSnap = 0; iSnap != snapshotCount; ++iSnap) {
    expand(podBasis, reduce(podBasis, snapshots[iSnap], podComponents), displac[iSnap]);
  }

  DistrVecBasis *veloc = 0;
  if(velocSnapshots) {
    veloc = new DistrVecBasis(velocSnapshots->vectorCount(), vectorSize());

    // Project velocity snapshots on POD basis to get training configurations
    for (int iSnap = 0; iSnap != velocSnapshots->vectorCount(); ++iSnap) {
      expand(podBasis, reduce(podBasis, (*velocSnapshots)[iSnap], podComponents), (*veloc)[iSnap]);
    }
    delete velocSnapshots;
  }

  DistrVecBasis *accel = 0;
  if(accelSnapshots) {
    accel = new DistrVecBasis(accelSnapshots->vectorCount(), vectorSize());

    // Project acceleration snapshots on POD basis to get training configurations
    for (int iSnap = 0; iSnap != accelSnapshots->vectorCount(); ++iSnap) {
      expand(podBasis, reduce(podBasis, (*accelSnapshots)[iSnap], podComponents), (*accel)[iSnap]);
    }
    delete accelSnapshots;
  }

  SubElementSamplingDriver **subDrivers = new SubElementSamplingDriver * [decDomain_->getNumSub()];
  Vector *solutions = new Vector[decDomain_->getNumSub()];
  int numCPUs = (structCom) ? structCom->numCPUs() : 1;
  int myID = (structCom) ? structCom->myID() : 0;
  bool verboseFlag = (myID == 0); // output to the screen only for subdomains assigned to mpi process with rank 0
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
  for(int i=0; i<decDomain_->getNumSub(); ++i) {
    subDrivers[i] = new SubElementSamplingDriver(decDomain_->getAllSubDomains()[i]);

    VecBasis &subPodBasis = subDrivers[i]->podBasis();
    subPodBasis.dimensionIs(podVectorCount, subDrivers[i]->vectorSize());
    for(int j=0; j<podVectorCount; ++j) {
      subPodBasis[j] = StackVector(podBasis[j].subData(i), podBasis[j].subLen(i));
    }
    VecBasis &subDisplac = subDrivers[i]->displac();
    subDisplac.dimensionIs(snapshotCount, subDrivers[i]->vectorSize());
    for(int j=0; j<snapshotCount; ++j) {
      subDisplac[j] = StackVector(displac[j].subData(i), displac[j].subLen(i));
    }
    subDrivers[i]->timeStampsIs(timeStamps);
    if(veloc) {
      VecBasis *subVeloc = subDrivers[i]->veloc();
      subVeloc->dimensionIs(snapshotCount, subDrivers[i]->vectorSize());
      for(int j=0; j<snapshotCount; ++j) {
        (*subVeloc)[j] = StackVector((*veloc)[j].subData(i), (*veloc)[j].subLen(i));
      }
    }
    if(accel) {
      VecBasis *subAccel = subDrivers[i]->accel();
      subAccel->dimensionIs(snapshotCount, subDrivers[i]->vectorSize());
      for(int j=0; j<snapshotCount; ++j) {
        (*subAccel)[j] = StackVector((*accel)[j].subData(i), (*accel)[j].subLen(i));
      }
    }

    subDrivers[i]->computeSolution(solutions[i], verboseFlag);
  }

  for(int cpu = 0; cpu < numCPUs; ++cpu) {
    if(cpu == myID) {
      for(int i=0; i<decDomain_->getNumSub(); ++i) {
        subDrivers[i]->postProcess(solutions[i], (myID == 0 && i == 0), verboseFlag);
        delete subDrivers[i];
      }
    }
    if(structCom) structCom->sync();
  }

  delete [] subDrivers;
  delete [] solutions;
  if(veloc) delete veloc;
  if(accel) delete accel;
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrElementSamplingDriverNew(Domain *domain) {
  return new Rom::DistrElementSamplingDriver(domain, structCom);
}
