#include "DistrExplicitDEIMPodProjectionNonLinDynamic.h"

#include <Driver.d/DecDomain.h>
#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>
#include <Threads.d/PHelper.h>

#include <Driver.d/GeoSource.h>

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "VecNodeDof6Map.h"
#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"
#include "PtrPtrIterAdapter.h"

#include <utility>
#include <cstddef>
#include <algorithm>

#include <sys/time.h>

extern GeoSource *geoSource;

namespace Rom {

DistrExplicitDEIMPodProjectionNonLinDynamic::DistrExplicitDEIMPodProjectionNonLinDynamic(Domain *domain) :
  DistrExplicitLumpedPodProjectionNonLinDynamic(domain)
{}

void
DistrExplicitDEIMPodProjectionNonLinDynamic::preProcess() {

  DistrExplicitLumpedPodProjectionNonLinDynamic::preProcess();

  lin_fInt     = new DistrVector(MultiDomainDynam::solVecInfo());
  kelArrayCopy = new FullSquareMatrix*[decDomain->getNumSub()];

  buildInterpolationBasis();
}

void
DistrExplicitDEIMPodProjectionNonLinDynamic::getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex) {

  normalizedBasis_.fullExpand( d, *d_n);
  lin_fInt->zero();
  execParal2R(decDomain->getNumSub(),this,&DistrExplicitDEIMPodProjectionNonLinDynamic::subGetKtimesU,*d_n,*lin_fInt);
  execParal3R(decDomain->getNumSub(),this,&DistrExplicitDEIMPodProjectionNonLinDynamic::subGetWeightedInternalForceOnly,*fInt,t,tIndex);
  

  if(domain->solInfo().stable && domain->solInfo().isNonLin() && tIndex%domain->solInfo().stable_freq == 0) {
    GenMDDynamMat<double> ops;
    ops.K = K;
    decDomain->rebuildOps(ops, 0.0, 0.0, 0.0, kelArray);
  }
  
  if (domain->solInfo().filterFlags) {
    trProject(*fInt);
  }

  *a_n = *fInt - *fExt;
//   *a_n = *lin_fInt - *fExt;

  if(haveRot) {
    execParal2R(decDomain->getNumSub(),this,&DistrExplicitDEIMPodProjectionNonLinDynamic::subTransformWeightedNodesOnly,*a_n,3);
    fullMassSolver->reSolve(*a_n);
    execParal2R(decDomain->getNumSub(),this,&DistrExplicitDEIMPodProjectionNonLinDynamic::subTransformWeightedNodesOnly,*a_n,2);
    DistrVector toto(*a_n);
    dynMat->M->mult(toto, *a_n);
  }

   DistrVector dummy(solVecInfo());
   DistrVector nlin_fInt(MultiDomainDynam::solVecInfo());
   DistrVector linMin_fExt(MultiDomainDynam::solVecInfo());

   nlin_fInt   = *fInt - *lin_fInt;
   linMin_fExt = *fExt - *lin_fInt;

//  normalizedBasis_.reduce(*a_n,f);
  normalizedBasis_.reduce(linMin_fExt,dummy);
  deimBasis_.compressedVecReduce(nlin_fInt,f);
  f -= dummy;
  //  the residual is computed in this step to avoid projecting into the reduced coordinates twice

}

void
DistrExplicitDEIMPodProjectionNonLinDynamic::subGetKtimesU(int isub, DistrVector &d, DistrVector &f)
{
  SubDomain *sd = decDomain->getSubDomain(isub);
  StackVector subf(f.subData(isub), f.subLen(isub));
  StackVector subd(d.subData(isub), d.subLen(isub));
  sd->getKtimesU(subd, (double *) 0, subf, 1.0, (kelArrayCopy) ? kelArrayCopy[isub] : (FullSquareMatrix *) 0);
}

void
DistrExplicitDEIMPodProjectionNonLinDynamic::buildInterpolationBasis() {

 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, BasisId::FORCE, BasisId::POD);
 fileName = fileName + ".deim";

 DistrBasisInputFile BasisFile(fileName); //read in mass-normalized basis
 filePrint(stderr, " ... Reading Interpolation basis from file %s ...\n", fileName.c_str());
 const int interpBasisSize = domain->solInfo().maxSizePodRom ?
                             std::min(domain->solInfo().maxSizePodRom, BasisFile.stateCount()) :
                             BasisFile.stateCount();

 filePrint(stderr, " ... Interplation subspace of dimension = %d ...\n", interpBasisSize);
 deimBasis_.dimensionIs(interpBasisSize,decDomain->masterSolVecInfo());

 DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());

 typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
 DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                  SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
 DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

 for (DistrVecBasis::iterator it = deimBasis_.begin(),
                              it_end = deimBasis_.end();
                              it != it_end; ++it) {
   assert(BasisFile.validCurrentState());
   
  
   BasisFile.currentStateBuffer(buffer);

   converter.vector(buffer, *it);
   BasisFile.currentStateIndexInc();
 }
 
 std::vector< std::vector<std::pair<int,int> > > maskedIndicesBuf;
 maskedIndicesBuf.resize(decDomain->getNumSub());
 execParal1R(decDomain->getNumSub(),this,&DistrExplicitDEIMPodProjectionNonLinDynamic::subBuildInterpolationBasis, maskedIndicesBuf); 

 filePrint(stderr," ...Compressing Interpolation Basis...\n");
 DofSetArray **all_cdsa = new DofSetArray * [decDomain->getNumSub()];
 for(int i=0; i<decDomain->getNumSub(); ++i) {all_cdsa[i] = decDomain->getSubDomain(i)->getCDSA();}
 deimBasis_.makeSparseBasis(maskedIndicesBuf, all_cdsa); 

}

void
DistrExplicitDEIMPodProjectionNonLinDynamic::subBuildInterpolationBasis(int iSub, std::vector< std::vector<std::pair<int,int> > > &maskedIndicesBuf) {

  SubDomain *sd = decDomain->getSubDomain(iSub);
  //std::map<int, int> &subMaskedIndicesBuf = maskedIndicesBuf[iSub];
  std::vector<std::pair<int,int> > &subMaskedIndicesBuf = maskedIndicesBuf[iSub];

  for (GeoSource::NodeDofPairVec::const_iterator it = geoSource->nodeDofSlotBegin(),
                                                   it_end = geoSource->nodeDofSlotEnd();
                                                   it != it_end; ++it) {
    
     const int nodeId = it->first;
     const int packedId = sd->globalToLocal(nodeId);
   
     if(packedId < 0) {continue;}

//     subMaskedIndicesBuf.insert(subMaskedIndicesBuf.end(), std::make_pair(packedId,it->second));
     subMaskedIndicesBuf.push_back(std::make_pair(packedId,it->second));

  }

  sd->createKelArray(kelArrayCopy[iSub]);
 
/*  std::cout<<"num row = " << kelArray[iSub][0].numRow()<< std::endl;
  std::cout<<"num col = " << kelArray[iSub][0].numCol()<< std::endl;
  std::cout<<"kelArray = " << std::endl;
  for(int iele = 0; iele != sd->numElements(); iele++){
    for(int col = 0; col != kelArrayCopy[iSub][0].numCol(); col++){
      for(int row = 0; row != kelArrayCopy[iSub][0].numRow(); row++){
         filePrint(stderr,"% 3.6e ",kelArrayCopy[iSub][0][row][col]);}
      filePrint(stderr,"\n");}
    filePrint(stderr,"\n");}*/

}

} // end namespace Rom
