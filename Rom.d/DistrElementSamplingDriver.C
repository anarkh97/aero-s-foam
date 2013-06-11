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

DistrElementSamplingDriver::DistrElementSamplingDriver(Domain *domain, Communicator *comm) :
  domain_(domain),
  comm_(comm)
{}

void
DistrElementSamplingDriver::solve() {
  std::auto_ptr<DecDomain> decDomain(createDecDomain<double>(domain_));
  decDomain->preProcess();

 /* TODO need to read basis into a distributed vector and send each subvector to its SubElementSamplingDriver
         likewise for snapshots
  // read basis
  DistrVecBasis podBasis;
  FileNameInfo fileInfo;
  DistrBasisInputFile podBasisFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD));

  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, podBasisFile.stateCount()) :
                                     podBasisFile.stateCount();

  filePrint(stderr, " ... Projection subspace of dimension = %d ...\n", projectionSubspaceSize);
  podBasis.dimensionIs(projectionSubspaceSize, decDomain->masterSolVecInfo());

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  for (DistrVecBasis::iterator it = podBasis.begin(),
                               it_end = podBasis.end();
                               it != it_end; ++it) {
    assert(podBasisFile.validCurrentState());

    podBasisFile.currentStateBuffer(buffer);
    converter.vector(buffer, *it);

    podBasisFile.currentStateIndexInc();
  }

  // read snapshots
  DistrVecBasis snapshots;
  std::vector<double> timeStamps;
  {
    DistrBasisInputFile in(BasisFileId(fileInfo, BasisId::STATE, BasisId::SNAPSHOTS));
    const int skipFactor = 1; //domain->solInfo().skipPodRom;
    const int skipOffSet = 0; //domain->solInfo().skipOffSet;
    const int basisStateCount = 1 + (in.stateCount() - 1) / skipFactor;
    filePrint(stderr, " ... basisStateCount = %d ...\n", basisStateCount);

    snapshots.dimensionIs(basisStateCount, decDomain->masterSolVecInfo());
    timeStamps.reserve(basisStateCount);

    for (DistrVecBasis::iterator it = snapshots.begin(),
                                 it_end = snapshots.end();
                                 it != it_end; ++it) {
      assert(in.validCurrentState());

      in.currentStateBuffer(buffer);
      converter.vector(buffer, *it);

      in.currentStateIndexInc();
    }
  }
*/
  SubElementSamplingDriver **subDrivers = new SubElementSamplingDriver * [decDomain->getNumSub()];
  Vector *solutions = new Vector[decDomain->getNumSub()];
  for(int i=0; i<decDomain->getNumSub(); ++i) {
    subDrivers[i] = new SubElementSamplingDriver(decDomain->getAllSubDomains()[i]);
    subDrivers[i]->getSolution(solutions[i]);
  }

  int numCPUs = (structCom) ? structCom->numCPUs() : 1;
  int myID = (structCom) ? structCom->myID() : 0;
  for(int cpu = 0; cpu < numCPUs; ++cpu) {
    if(cpu == myID) {
      for(int i=0; i<decDomain->getNumSub(); ++i) {
        subDrivers[i]->postProcess(solutions[i], (myID == 0 && i==0));
        delete subDrivers[i];
      }
    }
    if(structCom) structCom->sync();
  }

  delete [] subDrivers;
  delete [] solutions;
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrElementSamplingDriverNew(Domain *domain) {
  return new Rom::DistrElementSamplingDriver(domain, structCom);
}
