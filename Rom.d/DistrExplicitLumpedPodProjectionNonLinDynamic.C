#include "DistrExplicitLumpedPodProjectionNonLinDynamic.h"

#include <Driver.d/DecDomain.h>
#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>
#include <Threads.d/PHelper.h>

#include <Driver.d/GeoSource.h>

#include <utility>
#include <cstddef>

extern GeoSource *geoSource;

namespace Rom {

DistrExplicitLumpedPodProjectionNonLinDynamic::DistrExplicitLumpedPodProjectionNonLinDynamic(Domain *domain) :
  DistrExplicitPodProjectionNonLinDynamicBase(domain)
{}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::preProcess() {
  DistrExplicitPodProjectionNonLinDynamicBase::preProcess();

  buildPackedElementWeights();
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex) {
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
