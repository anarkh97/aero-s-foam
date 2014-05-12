#include "LumpedPodProjectionNonLinDynamic.h"
#include "PodProjectionSolver.h"

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
                                                             FullSquareMatrix *melArray, bool forceOnly) {
  if(forceOnly) {
    domain->getWeightedInternalForceOnly(packedElementWeights_,
                                         geomState, elementInternalForce,
                                         allCorot, kelArray,
                                         residual, lambda, time, refState, melArray);
  }
  else {
    domain->getWeightedStiffAndForceOnly(packedElementWeights_,
                                         geomState, elementInternalForce,
                                         allCorot, kelArray,
                                         residual, lambda, time, refState, melArray);
  }
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
      Element *ele = domain->getElementSet()[packedId]; // get weighted element data
      std::vector<int> node_buffer(ele->numNodes());
      packedElementWeights_.insert(packedElementWeights_.end(), std::make_pair(packedId, weight));
      //put nodes for weighted element into dummy vector and insert into packed node vector
      ele->nodes(node_buffer.data());
      packedWeightedNodes_.insert(packedWeightedNodes_.end(), node_buffer.begin(), node_buffer.end());
    }
  }

  // XXX also need to add to packedWeightedNodes the nodes of any elements to which a follower force has been applied
  // if the follower forces are not reduced. See: DistrExplicitLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights

  //sort nodes in ascending order and erase redundant nodes
  std::sort(packedWeightedNodes_.begin(), packedWeightedNodes_.end());
  std::vector<int>::iterator packedNodeIt = std::unique(packedWeightedNodes_.begin(), packedWeightedNodes_.end());
  packedWeightedNodes_.resize(packedNodeIt-packedWeightedNodes_.begin());

  int elemCounter = packedElementWeights_.size();
  filePrint(stderr, " ... # Elems. in Reduced Mesh = %-4d...\n", elemCounter);

  if(elemCounter < domain->numElements()) {
    filePrint(stderr, " ... Compressing Basis              ...\n");
    GenVecBasis<double> &projectionBasis = const_cast<GenVecBasis<double> &>(dynamic_cast<GenPodProjectionSolver<double>*>(solver)->projectionBasis());
    projectionBasis.makeSparseBasis(packedWeightedNodes_, domain->getCDSA());
  }
  else {
    if(!domain->solInfo().useMassNormalizedBasis) {
      // XXX to support this case we would need to precompute offline and the load online the reduced mass matrix,
      //     or we could hyper reduce the entire inertial force vector rather than just the nonlinear part
      filePrint(stderr, " *** ERROR: \"use_mass_normalized_basis off\" is not supported for\n"
                        "     for model III when \"samplmsh.elementmesh.inc\" file is used.\n"
                        "     Aborting ...\n");
      exit(-1);
    }
  }
}

} /* end namespace Rom */
