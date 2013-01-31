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
  DistrExplicitPodProjectionNonLinDynamicBase(domain)
{}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::preProcess() {
  
  DistrExplicitPodProjectionNonLinDynamicBase::preProcess();

  buildPackedElementWeights();
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::updateDisplacement(DistrVector& temp1, DistrVector& d_n1) {

  normalizedBasis_.projectUp( temp1, *d_n);
  execParal1R(decDomain->getNumSub(),this,&DistrExplicitLumpedPodProjectionNonLinDynamic::subUpdateWeightedNodesOnly,*d_n);
  d_n1 = temp1;
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex) {

  execParal3R(decDomain->getNumSub(),this,&DistrExplicitLumpedPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly,*fInt,t,tIndex);
  
  if (domain->solInfo().filterFlags) {
    trProject(*fInt);
  }

  *a_n = *fInt - *fExt;
  normalizedBasis_.projectDown(*a_n,f);
  //  the residual is computed in this step to avoid projecting into the reduced coordinates twice

}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::buildPackedElementWeights() {
  packedElementWeights_.resize(decDomain->getNumSub());
  packedWeightedNodes_.resize(decDomain->getNumSub());  
  execParal(decDomain->getNumSub(),
            this, &DistrExplicitLumpedPodProjectionNonLinDynamic::subBuildPackedElementWeights);
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subUpdateWeightedNodesOnly(int iSub, DistrVector &v) {
  StackVector vec(v.subData(iSub), v.subLen(iSub));
  GeomState *gs = (*geomState).getSubGeomState(iSub);
  gs->update(vec, packedWeightedNodes_[iSub], 1);
}

void
DistrExplicitLumpedPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly(int iSub, DistrVector &f, double &t, int &tIndex) {
  SubDomain *sd = decDomain->getSubDomain(iSub);
  Vector residual(f.subLen(iSub), 0.0);
  Vector eIF(sd->maxNumDOF()); // eIF = element internal force for one element (a working array)

  if(domain->solInfo().stable && domain->solInfo().isNonLin() && tIndex%domain->solInfo().stable_freq == 0) {
    sd->getWeightedStiffAndForceOnly(packedElementWeights_[iSub], *(*geomState)[iSub], eIF,
                                     allCorot[iSub], kelArray[iSub], residual,
                                     1.0, t, (*geomState)[iSub]); // residual -= internal force);
  }
  else {
    sd->getWeightedInternalForceOnly(packedElementWeights_[iSub], *(*geomState)[iSub], eIF,
                                     allCorot[iSub], kelArray[iSub], residual,
                                     1.0, t, (*geomState)[iSub]); // residual -= internal force);
  }
  StackVector subf(f.subData(iSub), f.subLen(iSub));
  subf.linC(residual, -1.0); // f = -residual
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
      std::vector<int> node_buffer; node_buffer.resize(ele->numNodes()); //resize node buffer
      subElementWeights.insert(subElementWeights.end(), std::make_pair(packedId, weight)); //pack element weight
      //put nodes for weighted element into dummy vector and insert into packed node vector
      ele->nodes(node_buffer.data());
      subWeightedNodes.insert(subWeightedNodes.end(), node_buffer.begin(), node_buffer.end());
    }
  }
  //sort nodes in ascending order and erase redundant nodes
  std::sort(subWeightedNodes.begin(), subWeightedNodes.end());
  std::vector<int>::iterator packedNodeIt = std::unique(subWeightedNodes.begin(),subWeightedNodes.end());
  subWeightedNodes.resize(packedNodeIt-subWeightedNodes.begin());

}

MDDynamMat *
DistrExplicitLumpedPodProjectionNonLinDynamic::buildOps(double mCoef, double cCoef, double kCoef) {

  MDDynamMat *result = DistrExplicitPodProjectionNonLinDynamicBase::buildOps(mCoef, cCoef, kCoef);

  if(decDomain->getNumSub() == 1){
    normalizedBasis_.makeSparseBasis(packedWeightedNodes_[0], decDomain->getSubDomain(0)->getCDSA());}
  else {
    filePrint(stderr,"\n ... Must run 1 subdomain per MPI process! ... \n");
    exit(-1);}

  return result;
}

} // end namespace Rom
