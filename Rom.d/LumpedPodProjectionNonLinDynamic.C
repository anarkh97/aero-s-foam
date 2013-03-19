#include "LumpedPodProjectionNonLinDynamic.h"

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>

#include <utility>

extern GeoSource *geoSource;

namespace Rom {

LumpedPodProjectionNonLinDynamic::LumpedPodProjectionNonLinDynamic(Domain *domain) :
  PodProjectionNonLinDynamic(domain)
{}

void
LumpedPodProjectionNonLinDynamic::preProcess() {
  PodProjectionNonLinDynamic::preProcess();

  buildPackedElementWeights();
}

void
LumpedPodProjectionNonLinDynamic::getStiffAndForceFromDomain(GeomState &geomState, Vector &elementInternalForce,
                                                             Corotator **allCorot, FullSquareMatrix *kelArray,
                                                             Vector &residual, double lambda, double time, GeomState *refState,
                                                             FullSquareMatrix *melArray) {
  domain->getWeightedStiffAndForceOnly(packedElementWeights_,
                                       geomState, elementInternalForce,
                                       allCorot, kelArray,
                                       residual, lambda, time, refState, melArray);
}

void
LumpedPodProjectionNonLinDynamic::updateStates(GeomState *refState, GeomState& geomState)
{
  domain->updateWeightedElemStatesOnly(packedElementWeights_,
                                       refState, geomState, allCorot);
}

void
LumpedPodProjectionNonLinDynamic::buildPackedElementWeights() {
  for (GeoSource::ElementWeightMap::const_iterator it = geoSource->elementLumpingWeightBegin(),
                                                   it_end = geoSource->elementLumpingWeightEnd();
       it != it_end; ++it) {
    const int elemId = it->first;

    const int packedId = geoSource->glToPackElem(elemId);
    if (packedId < 0) {
      continue;
    }

    const double weight = it->second;
    if (weight != 0.0) {
      packedElementWeights_.insert(packedElementWeights_.end(), std::make_pair(packedId, weight));
    }
  }
}

} /* end namespace Rom */
