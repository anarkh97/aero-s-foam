#include "LumpedPodProjectionNonLinDynamic.h"
#include "PodProjectionSolver.h"

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>

#include <utility>

extern GeoSource *geoSource;
extern int verboseFlag;

namespace Rom {

LumpedPodProjectionNonLinDynamic::LumpedPodProjectionNonLinDynamic(Domain *domain) :
  PodProjectionNonLinDynamic(domain),
  localReducedMeshId_(0)
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
    domain->getWeightedInternalForceOnly(packedElementWeights_[localReducedMeshId_],
                                         geomState, elementInternalForce,
                                         allCorot, kelArray,
                                         residual, lambda, time, refState, melArray);
  }
  else {
    domain->getWeightedStiffAndForceOnly(packedElementWeights_[localReducedMeshId_],
                                         geomState, elementInternalForce,
                                         allCorot, kelArray,
                                         residual, lambda, time, refState, melArray);
  }
}

void
LumpedPodProjectionNonLinDynamic::updateStates(ModalGeomState *refState, ModalGeomState& geomState, double time)
{
  if((!domain->solInfo().getNLInfo().linearelastic && (geomState_Big->getHaveRot() || geomState_Big->getTotalNumElemStates() > 0))
     || domain->solInfo().readInROBorModes.size() > 1) {
    // updateStates is called after midpoint update (i.e. once per timestep)
    // so it is a convenient place to update and copy geomState_Big, if necessary
    const GenVecBasis<double> &projectionBasis = dynamic_cast<GenPodProjectionSolver<double>*>(solver)->projectionBasis();
    if(domain->solInfo().readInROBorModes.size() == 1) {
      // note: for local bases method, geomState_Big has already been updated in setLocalBasis
      Vector q_Big(NonLinDynamic::solVecInfo());
      projectionBasis.expand(geomState.q, q_Big);
      geomState_Big->explicitUpdate(domain->getNodes(), q_Big);
    }

    if(geomState_Big->getHaveRot()) {
      Vector vel_Big(NonLinDynamic::solVecInfo()),
             acc_Big(NonLinDynamic::solVecInfo());
      projectionBasis.expand(geomState.vel, vel_Big);
      geomState_Big->setVelocity(vel_Big);
      projectionBasis.expand(geomState.acc, acc_Big);
      geomState_Big->setAcceleration(acc_Big);
    }

    if(geomState_Big->getTotalNumElemStates() > 0)
      domain->updateWeightedElemStatesOnly(packedWeightedElems_, refState_Big, *geomState_Big, allCorot, time);

    *refState_Big = *geomState_Big;
  }
}

void
LumpedPodProjectionNonLinDynamic::buildPackedElementWeights() {
  packedElementWeights_.resize(geoSource->elementLumpingWeightSize()+1);
  for (int j=0; j<geoSource->elementLumpingWeightSize(); ++j) {
    for (GeoSource::ElementWeightMap::const_iterator it = geoSource->elementLumpingWeightBegin(j),
                                                     it_end = geoSource->elementLumpingWeightEnd(j);
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
        packedElementWeights_[j].insert(packedElementWeights_[j].end(), std::make_pair(packedId, weight));
        //put nodes for weighted element into dummy vector and insert into packed node vector
        ele->nodes(node_buffer.data());
        packedWeightedNodes_.insert(packedWeightedNodes_.end(), node_buffer.begin(), node_buffer.end());
        packedWeightedElems_.insert(packedId);
      }
    }
    filePrint(stderr, " ... # Elems. in Reduced Mesh = %-4d...\n", packedElementWeights_[j].size());
  }

  // XXX also need to add to packedWeightedNodes the nodes of any elements to which a follower force has been applied
  // if the follower forces are not reduced. See: DistrExplicitLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights

  //sort nodes in ascending order and erase redundant nodes
  std::sort(packedWeightedNodes_.begin(), packedWeightedNodes_.end());
  std::vector<int>::iterator packedNodeIt = std::unique(packedWeightedNodes_.begin(), packedWeightedNodes_.end());
  packedWeightedNodes_.resize(packedNodeIt-packedWeightedNodes_.begin());

  if(packedWeightedElems_.size() < domain->numElements()) {
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

int
LumpedPodProjectionNonLinDynamic::selectLocalBasis(Vector &q)
{
  // TODO: fast selection alorithm from Section 3.3 of the following reference:
  // Amsallem, Zahr and Farhat. "Nonlinear model order reduction based on local reduced‐order bases."
  // International Journal for Numerical Methods in Engineering 92.10 (2012): 891-916.
  int j = rand()%domain->solInfo().readInROBorModes.size();
  if(verboseFlag) std::cerr << " ... Selecting local basis # " << j << "     ...\n";
  return j;
}

void
LumpedPodProjectionNonLinDynamic::setLocalReducedMesh(int j)
{
  for(std::map<int, double>::const_iterator it = packedElementWeights_[localReducedMeshId_].begin(),
                                        it_end = packedElementWeights_[localReducedMeshId_].end(); it != it_end; ++it) {
    const int iele = it->first;
    kelArray[iele].zero();
  }

  localReducedMeshId_ = std::min(geoSource->elementLumpingWeightSize()-1, j);
}

} /* end namespace Rom */
