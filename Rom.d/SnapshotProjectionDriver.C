#include "SnapshotProjectionDriver.h"

#include "VecBasis.h"
#include "VecBasisOps.h"
#include "BasisOps.h" 
#include "FileNameInfo.h"
#include "BasisFileStream.h"
#include "VecBasisFile.h"
#include "SimpleBuffer.h"
#include "RenumberingUtils.h"
#include "MeshDesc.h"

#include "VecBasisOps.h"

#include <Driver.d/GeoSource.h>
#include <Driver.d/Domain.h>
#include <Math.d/Vector.h>
#include <Timers.d/StaticTimers.h>
#include <Utils.d/Connectivity.h>
#include <Element.d/Element.h>
#include <Driver.d/SysState.h>

#include <cstddef>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <set>
#include <map>
#include <vector>
#include <memory>
#include <fstream>
#include <string>
#include <utility>

#include <cassert>
#include <iostream>

extern GeoSource *geoSource;

namespace Rom {

// Member functions
// ================
int
SnapshotProjectionDriver::elementCount() const {
  return domain->numElements();
}

int
SnapshotProjectionDriver::vectorSize() const {
  return domain->numUncon();
}

SnapshotProjectionDriver::SnapshotProjectionDriver(Domain *d) :
  SingleDomainDynamic(d),
  velocSnapshots(NULL),
  veloc_(NULL),
  accelSnapshots(NULL),
  accel_(NULL)
{}

SnapshotProjectionDriver::~SnapshotProjectionDriver() {
  if(veloc_) delete veloc_;
  if(accel_) delete accel_;
}

void
SnapshotProjectionDriver::postProcess() {

  const int snapshotCount = displac_.vectorCount();

  SDDynamPostProcessor *postProcessor = getPostProcessor();
  DynamMat *dMat = buildOps(1.0,0.0,0.0);
  Vector zero(solVecInfo(), 0.0);
  Vector *v = &zero, *a = &zero, *vp = &zero, *externalForce = &zero, *aeroForce = NULL;

  std::vector<double>::iterator timeStampIt = timeStamps_.begin();
  for (int iSnap = 0; iSnap != snapshotCount; ++iSnap) {
    filePrint(stderr,"\r %4.2f%% complete", double(iSnap)/double(snapshotCount)*100.);
    geomState->explicitUpdate(domain->getNodes(), displac_[iSnap]);
    if(veloc_) { 
      geomState->setVelocity((*veloc_)[iSnap], 2); 
      v = &(*veloc_)[iSnap];
    }
    if(accel_) {
      geomState->setAcceleration((*accel_)[iSnap], 2);
      a = &(*accel_)[iSnap]; 
    }
    SysState<Vector> systemState(displac_[iSnap], *v, *a, *vp);
    postProcessor->dynamOutput(iSnap, *timeStampIt, *dMat, // XXX note: iSnap is not always the correct timeStepIndex
                               *externalForce, aeroForce, systemState);
    timeStampIt++;
  }
 
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);
}

void
SnapshotProjectionDriver::solve() {
  preProcess();
  postProcess();
  compProjError();
}

void
SnapshotProjectionDriver::preProcess() {

  SingleDomainDynamic::preProcess();

  const FileNameInfo fileInfo;
  
  // Read order reduction data
  const VecNodeDof6Conversion vecDofConversion(*domain->getCDSA());
  assert(vectorSize() == vecDofConversion.vectorSize());

  //Read in non-normalized basis
  {
    BasisInputStream in(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD), vecDofConversion);
    const int podSizeMax = domain->solInfo().maxSizePodRom;
    if (podSizeMax != 0) {
      readVectors(in, podBasis_, podSizeMax);
    } else {
      readVectors(in, podBasis_);
    }
  }

  MGSVectors(podBasis_.data(),podBasis_.numVec(),podBasis_.size());

  // Read state snapshots
  {
    BasisInputStream in(BasisFileId(fileInfo, BasisId::STATE, BasisId::SNAPSHOTS), vecDofConversion);
    const int skipFactor = domain->solInfo().skipPodRom;
    const int skipOffSet = domain->solInfo().skipOffSet;
    const int basisStateCount = (in.size() % 2) + (in.size() - skipOffSet) / skipFactor;

    snapshots.dimensionIs(basisStateCount, in.vectorSize());
    timeStamps_.reserve(basisStateCount);

    int count = 0;
    int skipCounter = skipFactor-skipOffSet;
    while (count < basisStateCount) {
      std::pair<double, double *> data;
      data.second = snapshots[count].data();
      in >> data;
      assert(in);
      if (skipCounter == skipFactor) {
        timeStamps_.push_back(data.first);
        skipCounter = 1;
        ++count;
      } else {
        ++skipCounter;
      }
    }

    assert(timeStamps_.size() == basisStateCount);
  }

  // Read velocity snapshots
  if(domain->solInfo().velocPodRomFile != "") {
    std::vector<double> timeStamps;
    velocSnapshots = new VecBasis;
    BasisInputStream in(BasisFileId(fileInfo, BasisId::VELOCITY, BasisId::SNAPSHOTS), vecDofConversion);
    const int skipFactor = domain->solInfo().skipPodRom;
    const int skipOffSet = domain->solInfo().skipOffSet;
    const int basisStateCount = 1 + (in.size() - 1) / skipFactor;

    velocSnapshots->dimensionIs(basisStateCount, in.vectorSize());
    timeStamps.reserve(basisStateCount);

    int count = 0;
    int skipCounter = skipFactor-skipOffSet;
    while (count < basisStateCount) {
      std::pair<double, double *> data;
      data.second = (*velocSnapshots)[count].data();
      in >> data;
      assert(in);
      if (skipCounter == skipFactor) {
        timeStamps.push_back(data.first);
        skipCounter = 1;
        ++count;
      } else {
        ++skipCounter;
      }
    }

    assert(timeStamps.size() == basisStateCount);
    // TODO: check that timeStamps for velocity snapshots match state snapshots
  }

  // Read acceleration snapshots
  if(domain->solInfo().accelPodRomFile != "") {
    std::vector<double> timeStamps;
    accelSnapshots = new VecBasis;
    BasisInputStream in(BasisFileId(fileInfo, BasisId::ACCELERATION, BasisId::SNAPSHOTS), vecDofConversion);
    const int skipFactor = domain->solInfo().skipPodRom;
    const int skipOffSet = domain->solInfo().skipOffSet;
    const int basisStateCount = 1 + (in.size() - 1) / skipFactor;

    accelSnapshots->dimensionIs(basisStateCount, in.vectorSize());
    timeStamps.reserve(basisStateCount);

    int count = 0;
    int skipCounter = skipFactor-skipOffSet;
    while (count < basisStateCount) {
      std::pair<double, double *> data;
      data.second = (*accelSnapshots)[count].data();
      in >> data;
      assert(in);
      if (skipCounter == skipFactor) {
        timeStamps.push_back(data.first);
        skipCounter = 1;
        ++count;
      } else {
        ++skipCounter;
      }
    }

    assert(timeStamps.size() == basisStateCount);
    // TODO: check that timeStamps for acceleration snapshots match state snapshots
  }
  

  const int podVectorCount = podBasis_.vectorCount();
  const int snapshotCount = snapshots.vectorCount();

  // Temporary buffers shared by all iterations
  Vector podComponents(podVectorCount);

  // Project snapshots on POD basis
  filePrint(stderr," ... Projecting %d displacement snapshots onto basis of size %d ...\n",
            snapshotCount, podVectorCount);
  displac_.dimensionIs(snapshotCount, vectorSize());
  for (int iSnap = 0; iSnap != snapshotCount; ++iSnap) {
    expand(podBasis_, reduce(podBasis_, snapshots[iSnap], podComponents), displac_[iSnap]);
  }

  if(velocSnapshots) {
    veloc_ = new VecBasis(velocSnapshots->vectorCount(), vectorSize());

    // Project velocity snapshots on POD basis
    filePrint(stderr," ... Projecting %d velocity snapshots onto basis of size %d ...\n",
              velocSnapshots->vectorCount(), podVectorCount);
    for (int iSnap = 0; iSnap != velocSnapshots->vectorCount(); ++iSnap) {
      expand(podBasis_, reduce(podBasis_, (*velocSnapshots)[iSnap], podComponents), (*veloc_)[iSnap]);
    }
  //  delete velocSnapshots;
  }

  if(accelSnapshots) {
    accel_ = new VecBasis(accelSnapshots->vectorCount(), vectorSize());

    // Project acceleration snapshots on POD basis
    filePrint(stderr," ... Projecting %d acceleration snapshots onto basis of size %d ...\n",
              accelSnapshots->vectorCount(), podVectorCount);
    for (int iSnap = 0; iSnap != accelSnapshots->vectorCount(); ++iSnap) {
      expand(podBasis_, reduce(podBasis_, (*accelSnapshots)[iSnap], podComponents), (*accel_)[iSnap]);
    }
    //delete accelSnapshots;
  }

}

void
SnapshotProjectionDriver::compProjError() {

  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> dispError(snapshots.vectorSize(),snapshots.numVec());
  Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > dispBuf(snapshots.data(),snapshots.vectorSize(),snapshots.numVec());
  Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > pdispBuf(displac_.data(),displac_.vectorSize(),displac_.numVec());

  dispError = pdispBuf - dispBuf;

  if(domain->solInfo().PODerrornorm.size() == 0) {
     std::cerr << "...No displacement file specified, exiting..." << std::endl;
     exit(-1);
    }

  FILE * dispFile = fopen(domain->solInfo().PODerrornorm[0].c_str(),"w");
  filePrint(dispFile,"Number of Training Configurations: %d\n",snapshots.numVec());

  int numNorms = domain->solInfo().PODerrornorm.size();

  filePrint(dispFile,"Displacement Projection Error:\n");
  filePrint(dispFile,"Frobenius: %1.6e\n", dispError.norm()/dispBuf.norm());
  filePrint(dispFile,"   Snapshot   |      L_inf      |      L1       |        L2   \n");
  for (int i = 0; i != snapshots.numVec(); ++i){
    filePrint(dispFile,"     %d      ",i+1);
    for(int pnorm = 0; pnorm != 3; pnorm ++){
      if((pnorm) == 0){
          double Lperror = dispError.col(i).lpNorm<Eigen::Infinity>();
          Lperror = Lperror/dispBuf.col(i).lpNorm<Eigen::Infinity>()*100.;
          filePrint(dispFile,"     %1.6e ",Lperror);
      } else if(pnorm == 1){
          double Lperror = dispError.col(i).lpNorm<1>();
          Lperror = Lperror/dispBuf.col(i).lpNorm<1>()*100.;
          filePrint(dispFile,"    %1.6e ",Lperror);
      } else if(pnorm == 2) {
          double Lperror = dispError.col(i).lpNorm<2>();
          Lperror = Lperror/dispBuf.col(i).lpNorm<2>()*100.;
          filePrint(dispFile,"    %1.6e   ",(pnorm),Lperror);
      }
    }
    filePrint(dispFile,"\n");
  }

  if(velocSnapshots) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> velError(velocSnapshots->vectorSize(),velocSnapshots->numVec());
    Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > velBuf(velocSnapshots->data(),velocSnapshots->vectorSize(),velocSnapshots->numVec());
    Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > pvelBuf(veloc_->data(),veloc_->vectorSize(),veloc_->numVec());

    velError = pvelBuf - velBuf;

    if(domain->solInfo().PODerrornorm.size() < 2) {
     std::cerr << "...No velocity file specified, exiting..." << std::endl;
     exit(-1);
    }

    FILE * velFile = fopen(domain->solInfo().PODerrornorm[1].c_str(),"w");
    filePrint(velFile,"Number of Training Configurations: %d\n",velocSnapshots->numVec());

    filePrint(velFile,"Velocity Projection Error:\n");
    filePrint(velFile,"Frobenius: %1.6e\n", velError.norm()/velBuf.norm());
    filePrint(velFile,"   Snapshot   |      L_inf      |      L1       |        L2   \n");
    for (int i = 0; i != velocSnapshots->numVec(); ++i){
      filePrint(velFile,"     %d      ",i+1);
      for(int pnorm = 0; pnorm != 3; pnorm ++){
        if((pnorm) == 0){
            double Lperror = velError.col(i).lpNorm<Eigen::Infinity>();
            Lperror = Lperror/velBuf.col(i).lpNorm<Eigen::Infinity>()*100.;
            filePrint(velFile,"     %1.6e ",Lperror);
        } else if(pnorm == 1) {
            double Lperror = velError.col(i).lpNorm<1>();
            Lperror = Lperror/velBuf.col(i).lpNorm<1>()*100.;
            filePrint(velFile,"    %1.6e ",Lperror);
        } else if(pnorm == 2) {
            double Lperror = velError.col(i).lpNorm<2>();
            Lperror = Lperror/velBuf.col(i).lpNorm<2>()*100.;
            filePrint(velFile,"    %1.6e   ",(pnorm),Lperror);
        }
      }
      filePrint(velFile,"\n");
    }
  }

  if(accelSnapshots) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> accelError(accelSnapshots->vectorSize(),accelSnapshots->numVec());
    Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > accelBuf(accelSnapshots->data(),accelSnapshots->vectorSize(),accelSnapshots->numVec());
    Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > paccelBuf(accel_->data(),accel_->vectorSize(),accel_->numVec());

    accelError = paccelBuf - accelBuf;

    if(!domain->solInfo().PODerrornorm.size() < 3) {
     std::cerr << "...No acceleration file specified, exiting..." << std::endl;
     exit(-1);
    }

    FILE * accelFile = fopen(domain->solInfo().PODerrornorm[2].c_str(),"w");
    filePrint(accelFile,"Number of Training Configurations: %d\n",accelSnapshots->numVec());

    filePrint(accelFile,"Acceleration Projection Error\n");
    filePrint(accelFile,"Frobenius: %1.6e\n", accelError.norm()/accelBuf.norm());
    filePrint(accelFile,"   Snapshot   |      L_inf      |      L1       |        L2   \n");
    for (int i = 0; i != accelSnapshots->numVec(); ++i){
      filePrint(accelFile,"     %d      ",i+1);
      for(int pnorm = 0; pnorm != 3; pnorm ++){
        if((pnorm) == 0){
            double Lperror = accelError.col(i).lpNorm<Eigen::Infinity>();
            Lperror = Lperror/accelBuf.col(i).lpNorm<Eigen::Infinity>()*100.;
            filePrint(accelFile,"     %1.6e ",Lperror);
        } else if(pnorm == 1){
            double Lperror = accelError.col(i).lpNorm<1>();
            Lperror = Lperror/accelBuf.col(i).lpNorm<1>()*100.;
            filePrint(accelFile,"    %1.6e ",Lperror);
        } else if(pnorm == 2) {
            double Lperror = accelError.col(i).lpNorm<2>();
            Lperror = Lperror/accelBuf.col(i).lpNorm<2>()*100.;
            filePrint(accelFile,"    %1.6e   ",(pnorm),Lperror);
        }
      }
      filePrint(accelFile,"\n");
    }
  }

}

} // end namespace Rom

Rom::DriverInterface *snapshotProjectionDriverNew(Domain *d) {
  return new Rom::SnapshotProjectionDriver(d);
}
