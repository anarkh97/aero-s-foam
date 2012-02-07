#include "DistrExplicitLumpedPodProjectionNonLinDynamic.h"

#include "DistrGalerkinProjectionSolver.h"

#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"
#include "PtrPtrIterAdapter.h"

#include <Driver.d/DecDomain.h>
#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>
#include <Utils.d/DistHelper.h>
#include <Threads.d/PHelper.h>

#include <Driver.d/GeoSource.h>

#include <algorithm>
#include <stdexcept>
#include <cstddef>

extern GeoSource *geoSource;

namespace Rom {

DistrExplicitLumpedPodProjectionNonLinDynamic::DistrExplicitLumpedPodProjectionNonLinDynamic(Domain *domain) :
  MultiDomainDynam(domain)
{}

DistrExplicitLumpedPodProjectionNonLinDynamic::~DistrExplicitLumpedPodProjectionNonLinDynamic() {
  // Nothing to do
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::preProcess() {
  MultiDomainDynam::preProcess();
 
  FileNameInfo fileInfo; 
  DistrBasisInputFile podBasisFile(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD));

  const int projectionSubspaceSize = domain->solInfo().maxSizePodRom ?
                                     std::min(domain->solInfo().maxSizePodRom, podBasisFile.stateCount()) :
                                     podBasisFile.stateCount();

  filePrint(stderr, "Projection subspace of dimension = %d\n", projectionSubspaceSize);
  projectionBasis_.dimensionIs(projectionSubspaceSize, decDomain->masterSolVecInfo());

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());
  
  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  for (DistrVecBasis::iterator it = projectionBasis_.begin(),
                               it_end = projectionBasis_.end();
                               it != it_end; ++it) {
    assert(podBasisFile.validCurrentState());

    podBasisFile.currentStateBuffer(buffer);
    converter.vector(buffer, *it);
    
    podBasisFile.currentStateIndexInc();
  }

  buildPackedElementWeights();
}

MDDynamMat *
DistrExplicitLumpedPodProjectionNonLinDynamic::buildOps(double mCoef, double cCoef, double kCoef) {
  MDDynamMat *result = MultiDomainDynam::buildOps(mCoef, cCoef, kCoef);
  assert(result->M);

  const GenSubDOp<double> &fullMass = *(result->M);
  std::auto_ptr<DistrGalerkinProjectionSolver> solver(new DistrGalerkinProjectionSolver(fullMass));
  solver->projectionBasisIs(projectionBasis_);

  delete result->dynMat;
  result->dynMat = solver.release();

  return result;
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::getInternalForce(DistrVector &d, DistrVector &f, double t) {
  execParal2R(decDomain->getNumSub(),
              this, &DistrExplicitLumpedPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly,
              f, t);
  
  if (domain->solInfo().filterFlags) {
    trProject(f);
  }
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::buildPackedElementWeights() {
  packedElementWeights_.resize(decDomain->getNumSub());
  execParal(decDomain->getNumSub(),
            this, &DistrExplicitLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights);
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double t) {
  SubDomain *sd = decDomain->getSubDomain(iSub);
  Vector residual(f.subLen(iSub), 0.0);
  Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)
  
  sd->getWeightedStiffAndForceOnly(packedElementWeights_[iSub], *(*geomState)[iSub], eIF,
                                   allCorot[iSub], kelArray[iSub], residual,
                                   1.0, t, NULL); // residual -= internal force);
  StackVector subf(f.subData(iSub), f.subLen(iSub));
  subf.linC(residual, -1.0); // f = -residual
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights(int iSub) {
  SubDomain *sd = decDomain->getSubDomain(iSub);
  std::map<int, double> &subElementWeights = packedElementWeights_[iSub];
  
  for (GeoSource::ElementWeightMap::const_iterator it = geoSource->elementLumpingWeightBegin(),
                                                   it_end = geoSource->elementLumpingWeightEnd();
       it != it_end; ++it) {
    const int elemId = it->first;

    const int packedId = sd->glToPackElem(elemId);
    if (packedId < 0) {
      continue;
    }

    const double weight = it->second;
    if (weight != 0.0) {
      subElementWeights.insert(subElementWeights.end(), std::make_pair(packedId, weight));
    }
  }
}

} // end namespace Rom
