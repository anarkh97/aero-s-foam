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
LumpedPodProjectionNonLinDynamic::updateStates(ModalGeomState *refState, ModalGeomState& geomState)
{
  // updateStates is called after midpoint update, so this is a good time to update the velocity and acceleration in geomState_Big
  Vector q_Big(NonLinDynamic::solVecInfo()),
         vel_Big(NonLinDynamic::solVecInfo()),
         acc_Big(NonLinDynamic::solVecInfo());
  const GenVecBasis<double> &projectionBasis = dynamic_cast<GenPodProjectionSolver<double>*>(solver)->projectionBasis();
  projectionBasis.expand(geomState.q, q_Big);
  geomState_Big->explicitUpdate(domain->getNodes(), q_Big);
  projectionBasis.expand(geomState.vel, vel_Big);
  geomState_Big->setVelocity(vel_Big);
  projectionBasis.expand(geomState.acc, acc_Big);
  geomState_Big->setAcceleration(acc_Big);

  domain->updateWeightedElemStatesOnly(packedElementWeights_, refState_Big, *geomState_Big, allCorot);
  *refState_Big = *geomState_Big;
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
