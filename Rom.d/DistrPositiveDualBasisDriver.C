#include "DistrPositiveDualBasisDriver.h"

#include "DistrDomainUtils.h"
#include "DistrMasterMapping.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"
#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "DistrNonnegativeMatrixFactorization.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>
#include <cassert>
#include <sstream>

namespace Rom {

DistrPositiveDualBasisDriver::DistrPositiveDualBasisDriver(Domain *domain, Communicator *comm) :
  MultiDomainDynam(domain),
  domain_(domain),
  comm_(comm)
{}

//Non-member functions
//===============
void readIntoSolver(DistrNonnegativeMatrixFactorization &solver, DistrNodeDof1Buffer &inputBuffer, int& nodeCount,
                    DistrVecNodeDof1Conversion &converter, const DistrInfo &vectorSize,
                    BasisId::Level type, int numEntries, int& solverCol, int skipFactor=1)
{
  const BasisId::Type workload = BasisId::DUALSTATE;
  FileNameInfo fileInfo;
  for(int i = 0; i < numEntries; i++){
    DistrBasisInputFileTemplate<1> inputFile(BasisFileId(fileInfo,workload,type,i));
    nodeCount = inputFile.nodeCount();
    int basisStateCount = 1+(inputFile.stateCount()-1)/skipFactor;
    {
      int count = 0;
      int skipCounter = skipFactor;
      while (count < basisStateCount) {
        assert(inputFile.validCurrentState());
        inputFile.currentStateBuffer(inputBuffer);

        if (skipCounter >= skipFactor) {
          double *vecBuffer = solver.matrixColBuffer(solverCol);
          GenStackDistVector<double> vec(vectorSize, vecBuffer);

          converter.vector(inputBuffer, vec);
          //converter.paddedMasterVector(inputBuffer, vec);
          //std::fill(vecBuffer + vectorSize.totLen(), vecBuffer + solver.localRows(), 0.0);
          skipCounter = 1;
          ++solverCol;
          ++count;
        } else {
          ++skipCounter;
        }
        inputFile.currentStateIndexInc();
      }
    }
  }
}

void
DistrPositiveDualBasisDriver::solve() {
  
  MultiDomainDynam::preProcess();

  const int blockSize = domain->solInfo().svdBlockSize;

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  for(int i=0; i<decDomain->getNumSub(); ++i) decDomain->getSubDomain(i)->computeMasterFlag(decDomain->getMpcToSub());
  DistrTrivialMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                          SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()),
                                          domain->getNumCTC(), blockSize, comm_, decDomain->getNumSub());

  DistrVecNodeDof1Conversion converter(decDomain->getAllSubDomains(),
                                       decDomain->getAllSubDomains() + decDomain->getNumSub(),
                                       domain->getNumCTC(), blockSize, comm_, decDomain->getNumSub());

  FileNameInfo fileInfo;
  const BasisId::Type workload = BasisId::DUALSTATE;

  const int skipFactor = domain->solInfo().skipPodRom;
  int stateCount = 0;
  int nodeCount = 0;
  int snapBasisStateCount = 0;
  if(!domain->solInfo().snapfiPodRom.empty()) {
    for(int i = 0; i < domain->solInfo().snapfiPodRom.size(); i++) {
      std::string fileName = BasisFileId(fileInfo, workload, BasisId::SNAPSHOTS, i);
      DistrBasisInputFileTemplate<1> inputFile(fileName);
      snapBasisStateCount += 1+(inputFile.stateCount()-1)/skipFactor;
    }
  }

  DistrInfo distrInfo;
  decDomain->makeBlockCyclicDistrInfo(distrInfo, domain->getNumCTC(), blockSize);
  const int localLength = distrInfo.totLen();
  int maxBasisDimension = domain->solInfo().maxSizePodRom + (domain->solInfo().nmfDelROBDim)*(domain->solInfo().nmfNumROBDim-1);

  const int globalProbSize = domain->getNumCTC();
 
  if(domain->solInfo().snapfiPodRom.empty() && domain->solInfo().robfi.empty()) {
    filePrint(stderr, " *** ERROR: no files provided\n");
    exit(-1);
  }

  DistrNonnegativeMatrixFactorization solver(comm_, globalProbSize, snapBasisStateCount, localLength, maxBasisDimension, blockSize,
                                             domain->solInfo().nmfMaxIter, domain->solInfo().nmfTol);
 
  int solverCol = 0;
  DistrNodeDof1Buffer inputBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
  readIntoSolver(solver, inputBuffer, nodeCount, converter, distrInfo, BasisId::SNAPSHOTS,
                 domain->solInfo().snapfiPodRom.size(), solverCol, domain->solInfo().skipPodRom); // read in snapshots

  DistrNodeDof1Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
  for (int iBasis = 0; iBasis < domain->solInfo().nmfNumROBDim; ++iBasis) {
    int orthoBasisDim = domain->solInfo().maxSizePodRom + iBasis*domain->solInfo().nmfDelROBDim;
    filePrint(stderr, " ... Computation of a positive basis of size %d ...\n", orthoBasisDim);
    //solver.basisDimensionIs(orthoBasisDim);
    if (iBasis==0)
      solver.solve(0);
    else
      solver.solve(orthoBasisDim-domain->solInfo().nmfDelROBDim);

    std::string fileName = BasisFileId(fileInfo, workload, BasisId::POD);
    std::ostringstream ss;
    ss << orthoBasisDim;
    fileName.append(ss.str());
    DistrBasisOutputFile outputFile(fileName,
                                    nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(),
                                    comm_, false, 1);
    filePrint(stderr, " ... Writing positive basis to file %s ...\n", fileName.c_str());
    for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
      double * const vecBuffer = const_cast<double *>(solver.basisColBuffer(iVec));
      const GenStackDistVector<double> vec(distrInfo, vecBuffer);
      converter.paddedNodeDof6(vec, outputBuffer);
      outputFile.stateAdd(outputBuffer, 1.0);
    }
  }
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrPositiveDualBasisDriverNew(Domain *domain) {
  return new Rom::DistrPositiveDualBasisDriver(domain, structCom);
}
