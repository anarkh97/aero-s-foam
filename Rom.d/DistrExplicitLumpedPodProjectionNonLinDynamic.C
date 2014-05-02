#include "DistrExplicitLumpedPodProjectionNonLinDynamic.h"

#include <Driver.d/DecDomain.h>
#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>
#include <Threads.d/PHelper.h>

#include <Driver.d/GeoSource.h>

#include <utility>
#include <cstddef>

#include <sys/time.h>

extern GeoSource *geoSource;

namespace Rom {

DistrExplicitLumpedPodProjectionNonLinDynamic::DistrExplicitLumpedPodProjectionNonLinDynamic(Domain *domain) :
  DistrExplicitPodProjectionNonLinDynamicBase(domain), K(NULL)
{}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::preProcess() {
  
  DistrExplicitPodProjectionNonLinDynamicBase::preProcess();

  buildPackedElementWeights();

  if(domain->solInfo().stable) {
    execParal(decDomain->getNumSub(),this,&DistrExplicitLumpedPodProjectionNonLinDynamic::subInitWeightedStiffOnly);
  }
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::updateState(double dt_n_h, DistrVector& v_n_h, DistrVector& d_n1) {

  DistrVector temp1(solVecInfo());
  temp1 = dt_n_h*v_n_h;
  normalizedBasis_.expand( temp1, *d_n);
  execParal1R(decDomain->getNumSub(),this,&DistrExplicitLumpedPodProjectionNonLinDynamic::subUpdateWeightedNodesOnly,*d_n);
  d_n1 += temp1;

  if(haveRot) {
    normalizedBasis_.expand(v_n_h, *v_n);
    execParal1R(decDomain->getNumSub(),this,&DistrExplicitLumpedPodProjectionNonLinDynamic::subSetVelocityWeightedNodesOnly,*v_n);
  }
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex) {

  execParal3R(decDomain->getNumSub(),this,&DistrExplicitLumpedPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly,*fInt,t,tIndex);

  if(domain->solInfo().stable && domain->solInfo().isNonLin() && tIndex%domain->solInfo().stable_freq == 0) {
    GenMDDynamMat<double> ops;
    ops.K = K;
    decDomain->rebuildOps(ops, 0.0, 0.0, 0.0, kelArray);
  }
  
  if (domain->solInfo().filterFlags) {
    trProject(*fInt);
  }

  *tempVec = *fInt - *fExt;
  normalizedBasis_.sparseVecReduce(*tempVec, f);
  //  the residual is computed in this step to avoid projecting into the reduced coordinates twice
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::buildPackedElementWeights() {
  packedElementWeights_.resize(decDomain->getNumSub());
  packedWeightedNodes_.resize(decDomain->getNumSub());  
  execParal(decDomain->getNumSub(),
            this, &DistrExplicitLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights);

  int elemCounter = 0;
  for(int i=0; i<decDomain->getNumSub(); ++i) elemCounter += packedElementWeights_[i].size();
  if(structCom) elemCounter = structCom->globalSum(elemCounter);

  filePrint(stderr, " ... # Elems. in Reduced Mesh = %-4d...\n", elemCounter);

  if(elemCounter < domain->numElements()) {
    filePrint(stderr, " ... Compressing Basis              ...\n");
    DofSetArray **all_cdsa = new DofSetArray * [decDomain->getNumSub()];
    for(int i=0; i<decDomain->getNumSub(); ++i) all_cdsa[i] = decDomain->getSubDomain(i)->getCDSA();
    normalizedBasis_.makeSparseBasis(packedWeightedNodes_, all_cdsa);
    delete [] all_cdsa;
  }
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subUpdateWeightedNodesOnly(int iSub, DistrVector &v) {
  StackVector vec(v.subData(iSub), v.subLen(iSub));
  GeomState *gs = (*geomState).getSubGeomState(iSub);
  gs->update(vec, packedWeightedNodes_[iSub], 2);
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subSetVelocityWeightedNodesOnly(int iSub, DistrVector &v) {
  StackVector vec(v.subData(iSub), v.subLen(iSub));
  GeomState *gs = (*geomState).getSubGeomState(iSub);
  gs->setVelocity(packedWeightedNodes_[iSub].size(), packedWeightedNodes_[iSub].data(), vec, 2);
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subTransformWeightedNodesOnly(int iSub, DistrVector &v, int type) {
  StackVector vec(v.subData(iSub), v.subLen(iSub));
  GeomState *gs = (*geomState).getSubGeomState(iSub);
  gs->transform(vec, packedWeightedNodes_[iSub], type, true);
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double &t, int &tIndex) {
  SubDomain *sd = decDomain->getSubDomain(iSub);
  Vector residual(f.subLen(iSub), 0.0);
  Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)

  if(domain->solInfo().stable && domain->solInfo().isNonLin() && tIndex%domain->solInfo().stable_freq == 0) {
    sd->getWeightedStiffAndForceOnly(packedElementWeights_[iSub], *(*geomState)[iSub], eIF,
                                     allCorot[iSub], kelArray[iSub], residual,
                                     1.0, t, (*geomState)[iSub], melArray[iSub]); // residual -= internal force);
  }
  else {
    sd->getWeightedInternalForceOnly(packedElementWeights_[iSub], *(*geomState)[iSub], eIF,
                                     allCorot[iSub], kelArray[iSub], residual,
                                     1.0, t, (*geomState)[iSub], melArray[iSub]); // residual -= internal force);
  }
  StackVector subf(f.subData(iSub), f.subLen(iSub));
  subf.linC(residual, -1.0); // f = -residual
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subInitWeightedStiffOnly(int iSub) {
  SubDomain *sd = decDomain->getSubDomain(iSub);
  Vector residual(sd->numUncon(), 0.0);
  Vector eIF(sd->maxNumDOF());

  sd->getWeightedStiffAndForceOnly(packedElementWeights_[iSub], *(*geomState)[iSub], eIF,
                                   allCorot[iSub], kelArray[iSub], residual,
                                   1.0, 0.0, (*geomState)[iSub], melArray[iSub]);
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights(int iSub) {
  //need to sort out issue with subdomains for packed node vector
  SubDomain *sd = decDomain->getSubDomain(iSub);
  std::map<int, double> &subElementWeights = packedElementWeights_[iSub];

  std::vector<int> &subWeightedNodes = packedWeightedNodes_[iSub];
  
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
      Element *ele = sd->getElementSet()[packedId]; // get weighted element data
      std::vector<int> node_buffer(ele->numNodes());
      subElementWeights.insert(subElementWeights.end(), std::make_pair(packedId, weight)); //pack element weight
      //put nodes for weighted element into dummy vector and insert into packed node vector
      ele->nodes(node_buffer.data());
      subWeightedNodes.insert(subWeightedNodes.end(), node_buffer.begin(), node_buffer.end());
    }
  }

  if(!domain->solInfo().reduceFollower) {
    // add the nodes for all elements to which follower forces have been applied
    std::vector<int> &followedElemList = sd->getFollowedElemList();
    for(std::vector<int>::iterator it = followedElemList.begin(), it_end = followedElemList.end(); it != it_end; ++it) {
      Element *ele = sd->getElementSet()[*it]; 
      std::vector<int> node_buffer(ele->numNodes());
      ele->nodes(node_buffer.data());
      subWeightedNodes.insert(subWeightedNodes.end(), node_buffer.begin(), node_buffer.end());
    }
  }
/* XXX consider whether to also add the nodes with non-follower external forces
  for(int i = 0; i < sd->nNeumann(); ++i) {
    subWeightedNodes.push_back(sd->getNBC()[i].nnum);
  }
*/
  //sort nodes in ascending order and erase redundant nodes
  std::sort(subWeightedNodes.begin(), subWeightedNodes.end());
  std::vector<int>::iterator packedNodeIt = std::unique(subWeightedNodes.begin(),subWeightedNodes.end());
  subWeightedNodes.resize(packedNodeIt-subWeightedNodes.begin());
}

MDDynamMat *
DistrExplicitLumpedPodProjectionNonLinDynamic::buildOps(double mCoef, double cCoef, double kCoef) {

  MDDynamMat *result = DistrExplicitPodProjectionNonLinDynamicBase::buildOps(mCoef, cCoef, kCoef);
  K = result->K;

  return result;
}

} // end namespace Rom
