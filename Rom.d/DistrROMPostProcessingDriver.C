#include "DistrROMPostProcessingDriver.h"

#include "DistrDomainUtils.h"
#include <Driver.d/Domain.h>
#include <Driver.d/DecDomain.h>

#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"
#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include "DistrMasterMapping.h"
#include "DistrVecNodeDof6Conversion.h"
#include "DistrNodeDof6Buffer.h"

#include <Feti.d/DistrVector.h>
#include <Utils.d/DistHelper.h>

#include "PtrPtrIterAdapter.h"

#include <algorithm>
#include <memory>

#include <iostream>
#include <cassert>

class Rbm;
 
namespace Rom {

DistrROMPostProcessingDriver::DistrROMPostProcessingDriver(Domain *domain_) :
MultiDomainDynam(domain_),
normalizedBasis_(),
curState(NULL), fullDispBuffer(NULL), fullVelBuffer(NULL), fullAccBuffer(NULL),
fullVel2Buffer(NULL), fullDummyBuffer(NULL),
dummyDynOps(NULL)
{}

DistrROMPostProcessingDriver::~DistrROMPostProcessingDriver()
{
 if(curState) delete curState;
 if(fullDispBuffer) delete fullDispBuffer;
 if(fullVelBuffer) delete fullVelBuffer;
 if(fullAccBuffer) delete fullAccBuffer;
 if(fullVel2Buffer) delete fullVel2Buffer;
 if(fullDummyBuffer) delete fullDummyBuffer;
 if(dummyDynOps) delete dummyDynOps;
}

void
DistrROMPostProcessingDriver::preProcess() {

  {MultiDomainDynam::preProcess();
  bufferReducedFiles();
  //initialized decDomain class for use in projection basis preprocessing

  // read in distribuited POD basis
  FileNameInfo fileInfo;
  std::string fileName = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
  fileName.append(".normalized");  
  DistrBasisInputFile podBasisFile(fileName);  //read in mass-normalized basis
  filePrint(stderr, " ... Reading basis from file %s ...\n", fileName.c_str());
  filePrint(stderr, " ... Projection subspace of dimension = %d ...\n", projectionSubspaceSize);
  normalizedBasis_.dimensionIs(projectionSubspaceSize, decDomain->masterSolVecInfo());

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(), decDomain->getAllSubDomains() + decDomain->getNumSub());

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));
  DistrNodeDof6Buffer buffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());

  for (DistrVecBasis::iterator it = normalizedBasis_.begin(),
                               it_end = normalizedBasis_.end();
                               it != it_end; ++it) {
    assert(podBasisFile.validCurrentState());

    podBasisFile.currentStateBuffer(buffer);
    converter.vector(buffer, *it);

    podBasisFile.currentStateIndexInc();
  }
  
  VECsize = normalizedBasis_.size();}

  {//initialize multi domain dynamic post processor
  mddPostPro = MultiDomainDynam::getPostProcessor();

  //initialize containers for full coordinates 
  fullDispBuffer  = new GenDistrVector<double>(MultiDomainDynam::solVecInfo());
  fullVelBuffer   = new GenDistrVector<double>(MultiDomainDynam::solVecInfo());
  fullAccBuffer   = new GenDistrVector<double>(MultiDomainDynam::solVecInfo());
  fullVel2Buffer  = new GenDistrVector<double>(MultiDomainDynam::solVecInfo());
  fullDummyBuffer = new GenDistrVector<double>(MultiDomainDynam::solVecInfo());

  //initialize system state vector container for use in mulit domain dynamic post processor
  curState = new SysState<GenDistrVector<double> >( *fullDispBuffer, *fullVelBuffer, *fullAccBuffer, *fullVel2Buffer);}
}  //end preProcessing

void 
DistrROMPostProcessingDriver::bufferReducedFiles(){

  //get output information needed for parsing reduced data
  numConversionFiles = decDomain->getDomain()->solInfo().numRODFile;
  int skipTime = decDomain->getDomain()->solInfo().skipPodRom;
  //loop over reduced coordinate files
  for(int i = 0; i < numConversionFiles; i++) {
    //have all threads parse the reduced coordinate input file
    //there should be plenty of memory per node since projectionSubspaceSize is small
    ifstream reducedCoordFile(decDomain->getDomain()->solInfo().RODConversionFiles[i].c_str());
    if(reducedCoordFile.is_open()) {
      if(skipTime > 1) filePrint(stderr, " ... Skipping every %3d snapshots   ...\n", skipTime);

      double time, dummyVar;
      int datatype, podsize, skipCounter;
      skipCounter = skipTime; // need to include t0 
      reducedCoordFile>>datatype; reducedCoordFile>>podsize;
      DataType.push_back(std::make_pair(datatype, podsize));
          switch(DataType[i].first) {
            case 0 :   // read reduced acceleration data
              {filePrint(stderr, " ... Buffering Reduced Acceleration Data ...\n");
              std::vector<double> timestamps;
              while(reducedCoordFile>>time) {
                if(skipCounter == skipTime) {
                  skipCounter = 1;
                  timestamps.push_back(time);
                  for(int j = 0; j < podsize; j++) {
                    reducedCoordFile>>dummyVar;
                    reducedAccBuffer.push_back(dummyVar);
                  }
                  filePrint(stderr,"\r Timestamp = %f", time);
                } else {
                  for(int j = 0; j < podsize; j++) 
                    reducedCoordFile>>dummyVar;
                  skipCounter += 1;
                }
              }
              TimeStamps.push_back(timestamps);}
              filePrint(stderr,"\n");
              break;
            case 1 :   // read reduced displacement data
              {filePrint(stderr, " ... Buffering Reduced Displacement Data ...\n");
              std::vector<double> timestamps;
              while(reducedCoordFile>>time) {
                if(skipCounter == skipTime) {
                  skipCounter = 1;
                  timestamps.push_back(time);
                  for(int j = 0; j < podsize; j++) {
                    reducedCoordFile>>dummyVar;
                    reducedDispBuffer.push_back(dummyVar);
                  }
                  filePrint(stderr,"\r Timestamp = %f", time);
                } else {
                  for(int j = 0; j < podsize; j++) 
                    reducedCoordFile>>dummyVar;
                  skipCounter += 1;
                }
              }
              TimeStamps.push_back(timestamps);}
              filePrint(stderr,"\n");
              break;
            case 2 :   // read reduced velocity data
              {filePrint(stderr, " ... Buffering Reduced Velocity Data ...\n");
              std::vector<double> timestamps;
              while(reducedCoordFile>>time) {
                if(skipCounter == skipTime){
                  skipCounter = 1;
                  timestamps.push_back(time);
                  for(int j = 0; j < podsize; j++) {
                    reducedCoordFile>>dummyVar;
                    reducedVelBuffer.push_back(dummyVar);
                  }
                  filePrint(stderr,"\r Timestamp = %f", time);
                } else {
                  for(int j = 0; j < podsize; j++)
                    reducedCoordFile>>dummyVar;
                  skipCounter += 1;
                }
              }
              TimeStamps.push_back(timestamps);}
              filePrint(stderr,"\n");
              break;
            default :
              filePrint(stderr, "\n... ROD conversion only supports Acceleration, Displacement, and Velocity ...\n");
          }
    } else {
      filePrint(stderr,"\nFailure to open file \n");
    }

    if(i != 0) { 
      if(DataType[i].second != DataType[i-1].second) {
        filePrint(stderr,"\n *** WARNING: Incompatible Input files %f %f\n", DataType[i-1].second, DataType[i].second);
        //exit(-1);
      }
    }

  } //end loop over input files, finished reading reduced data

  if(DataType.size() > 0) projectionSubspaceSize = DataType[0].second;
}

void
DistrROMPostProcessingDriver::solve() {

   preProcess();

   int counter = 0; //TODO: make this portion more general so it doesn't depend on th assumption
                    //that all files have matching timestamps
   if(TimeStamps.size() > 0)
   for(std::vector<double>::iterator it = TimeStamps[0].begin(); it != TimeStamps[0].end(); it++) {

     // load current state for output 
     for(int i = 0; i < numConversionFiles; i++) {
       std::vector<double> buffer(projectionSubspaceSize);
           switch(DataType[i].first) {
            case 0 :

              for (int j = 0; j < projectionSubspaceSize; j++) 
                buffer[j] = reducedAccBuffer[counter*projectionSubspaceSize+j];
             
              normalizedBasis_.projectUp(buffer, *fullAccBuffer);
              break;
            case 1 :

              for (int j = 0; j < projectionSubspaceSize; j++)
                buffer[j] = reducedDispBuffer[counter*projectionSubspaceSize+j];

              normalizedBasis_.projectUp(buffer, *fullDispBuffer);
              break;
            case 2 :
              if(counter != 0)
                *fullVel2Buffer = *fullVelBuffer;

              for (int j = 0; j < projectionSubspaceSize; j++) 
                buffer[j] = reducedVelBuffer[counter*projectionSubspaceSize+j];

              normalizedBasis_.projectUp(buffer, *fullVelBuffer);
              break;
            default :
              break; 
              //nothing
           }
     }

     geomState->explicitUpdate(decDomain, *fullDispBuffer);
     if(!dummyDynOps) dummyDynOps = new MDDynamMat;
     mddPostPro->dynamOutput(counter, *it, *dummyDynOps, *fullDummyBuffer, fullDummyBuffer, *curState);

     filePrint(stderr,"\r ... ROM Conversion Loop: t = %9.3e, %3d%% complete ...",
                *it, int(*it/(TimeStamps[0].back())*100));

     counter += 1;
   }  //end of loop over time stamps

}


} //end namespace Rom

Rom::DriverInterface *distrROMPostProcessingDriverNew(Domain *D) {
  return new Rom::DistrROMPostProcessingDriver(D);
}
