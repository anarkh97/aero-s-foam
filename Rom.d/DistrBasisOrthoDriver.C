#include "DistrBasisOrthoDriver.h"

#include "DistrDomainUtils.h"
#include "DistrMasterMapping.h"
#include "DistrNodeDof6Buffer.h"
#include "DistrVecNodeDof6Conversion.h"
#include "PtrPtrIterAdapter.h"
#include "DistrVecBasis.h"
#include "DistrVecBasisOps.h"

#include "DistrSvdOrthogonalization.h"

#include "FileNameInfo.h"
#include "DistrBasisFile.h"

#include <Utils.d/DistHelper.h>

#include <algorithm>
#include <memory>
#include <cassert>

namespace Rom {

DistrBasisOrthoDriver::DistrBasisOrthoDriver(Domain *domain, Communicator *comm) :
  MultiDomainDynam(domain),
  domain_(domain),
  comm_(comm)
{}

void
DistrBasisOrthoDriver::solve() {
  
  MultiDomainDynam::preProcess();

  typedef PtrPtrIterAdapter<SubDomain> SubDomIt;
  DistrMasterMapping masterMapping(SubDomIt(decDomain->getAllSubDomains()),
                                   SubDomIt(decDomain->getAllSubDomains() + decDomain->getNumSub()));

  DistrVecNodeDof6Conversion converter(decDomain->getAllSubDomains(),
                                       decDomain->getAllSubDomains() + decDomain->getNumSub());
  
  FileNameInfo fileInfo;
  const BasisId::Type workload = BasisId::STATE;

  const int skipFactor = domain->solInfo().skipPodRom;
  int stateCount = 0;
  int nodeCount = 0;
  int snapBasisStateCount = 0;
  if(!domain->solInfo().snapfiPodRom.empty()) {
    for(int i = 0; i < domain->solInfo().snapfiPodRom.size(); i++) {
      std::string fileName = BasisFileId(fileInfo, workload, BasisId::SNAPSHOTS, i);
      DistrBasisInputFile inputFile(fileName);
      stateCount += inputFile.stateCount();
      snapBasisStateCount += 1+(stateCount-1)/skipFactor;
    }
  }

  DistrSvdOrthogonalization solver(comm_, comm_->numCPUs(), 1);
  
  const int blockSize = 64; // TODO More flexible
  {
    solver.blockSizeIs(blockSize);
  }

  const int localLength = decDomain->solVecInfo().totLen();
  {
    const int maxLocalLength = comm_->globalMax(localLength);
    const int maxCpuLoad = ((maxLocalLength / blockSize) + (maxLocalLength % blockSize)) * blockSize;
    assert(maxCpuLoad >= localLength);
    const int globalProbSize = maxCpuLoad * solver.rowCpus();
    solver.problemSizeIs(globalProbSize, snapBasisStateCount);
    assert(solver.localRows() == maxCpuLoad);
  }

  double beta = domain->solInfo().newmarkBeta;
  // Assembling mass matrix
  MDDynamMat *dynOps = MultiDomainDynam::buildOps(1.0, 0.0, 0.0);
  assert(dynOps->M);
  assert(dynOps->dynMat);
 
  if(domain->solInfo().snapfiPodRom.empty() && domain->solInfo().robfi.empty()) {
    filePrint(stderr, " *** Error: no files provided\n");
    exit(-1);
  }
   
  int solverCol = 0;
  DistrNodeDof6Buffer inputBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
  for(int i = 0; i < domain->solInfo().snapfiPodRom.size(); i++) {
    DistrBasisInputFile inputFile(BasisFileId(fileInfo, workload, BasisId::SNAPSHOTS, i));
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
          GenStackDistVector<double> vec(decDomain->solVecInfo(), vecBuffer);

          converter.paddedMasterVector(inputBuffer, vec);
          // TODO Multiply by weighting factor if given in input file
          if(beta == 0 && domain->solInfo().normalize == 1) { // new method
            dynOps->dynMat->squareRootMult(vec);
          }

          std::fill(vecBuffer + localLength, vecBuffer + solver.localRows(), 0.0);

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

  for(int i = 0; i < domain->solInfo().robfi.size(); i++) {
    DistrBasisInputFile inputFile(BasisFileId(fileInfo, workload, BasisId::ROB, i));
    int basisStateCount = 1+(inputFile.stateCount()-1)/skipFactor;
    {
      int count = 0;
      int skipCounter = skipFactor;
      while(count < basisStateCount) {
        assert(inputFile.validCurrentState());
        inputFile.currentStateBuffer(inputBuffer);

        if (skipCounter >= skipFactor) {
          double *vecBuffer = solver.matrixColBuffer(solverCol);
          GenStackDistVector<double> vec(decDomain->solVecInfo(), vecBuffer);
            
          converter.paddedMasterVector(inputBuffer, vec);
          // TODO Multiply by weighting factor if given in input file
          if(beta == 0 && domain->solInfo().normalize == 1) { // new method
            dynOps->dynMat->squareRootMult(vec);
          }
          vec *= inputFile.currentStateHeaderValue(); // Multiply by the singular value stored in header
            
          std::fill(vecBuffer + localLength, vecBuffer + solver.localRows(), 0.0);

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

  solver.solve();

  // Output solution
  const int podVectorCount = domain_->solInfo().maxSizePodRom ?
                             std::min(domain_->solInfo().maxSizePodRom, solver.singularValueCount()) :
                             solver.singularValueCount();
  {
    DistrNodeDof6Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
    DistrBasisOutputFile outputFile(BasisFileId(fileInfo, workload, BasisId::POD),
                                    nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(),
                                    comm_, false);

    if(beta != 0 || (beta == 0 && domain->solInfo().normalize == 0))
      filePrint(stderr, " ... Writing orthonormal basis to file %s ...\n", BasisFileId(fileInfo, workload, BasisId::POD).name().c_str());
    for (int iVec = 0; iVec < podVectorCount; ++iVec) {
      double * const vecBuffer = const_cast<double *>(solver.basisColBuffer(iVec));
      const GenStackDistVector<double> vec(decDomain->solVecInfo(), vecBuffer);
      converter.paddedNodeDof6(vec, outputBuffer);
      outputFile.stateAdd(outputBuffer, solver.singularValue(iVec));
    }
  }
  comm_->sync(); // this is important, do not delete.

  // Renormalize basis for explicit
  if(beta == 0) {
    // Read back in output file to perform renormalization
    DistrVecBasis basis;
    {
      DistrBasisInputFile inputFile(BasisFileId(fileInfo, workload, BasisId::POD));
      DistrNodeDof6Buffer inputBuffer(masterMapping.localNodeBegin(), masterMapping.localNodeEnd());
      basis.dimensionIs(podVectorCount, decDomain->masterSolVecInfo()); 
      int i = 0;
      for(DistrVecBasis::iterator it = basis.begin(),
          it_end = basis.end();
          it != it_end; ++it) {
        assert(inputFile.validCurrentState());

        inputFile.currentStateBuffer(inputBuffer);
        converter.vector(inputBuffer, *it);

        inputFile.currentStateIndexInc();
      }
    }

    DistrVecBasis normalizedBasis;
    if(domain->solInfo().normalize == 0) {
      // Old method: renormalize the orthonormal basis
      renormalized_basis(*dynOps->M, basis, normalizedBasis);
    }
    if(domain->solInfo().normalize == 1) {
      // New method: multiply by inverse square root of the mass matrix
      for(int col = 0; col < podVectorCount; col++) {
        dynOps->dynMat->inverseSquareRootMult(basis[col]);
      }
      normalizedBasis = basis;
    }

    // Output the renormalized basis as separate file
    std::string fileName = BasisFileId(fileInfo, workload, BasisId::POD);
    fileName.append(".normalized");
    DistrNodeDof6Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
    DistrBasisOutputFile outputNormalizedFile(fileName, nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(), comm_, false);
    filePrint(stderr, " ... Writing mass-normalized basis to file %s ...\n", fileName.c_str());
    for (int iVec = 0; iVec < podVectorCount; ++iVec) {
      converter.paddedNodeDof6(normalizedBasis[iVec], outputBuffer);
      outputNormalizedFile.stateAdd(outputBuffer, solver.singularValue(iVec));
    }

    // Compute and output identity normalized basis if using new method
    if(domain->solInfo().normalize == 1) {
      MGSVectors(normalizedBasis);
      std::string fileName = BasisFileId(fileInfo, workload, BasisId::POD);
      DistrNodeDof6Buffer outputBuffer(masterMapping.masterNodeBegin(), masterMapping.masterNodeEnd());
      DistrBasisOutputFile outputOrthoNormalFile(fileName, nodeCount, outputBuffer.globalNodeIndexBegin(), outputBuffer.globalNodeIndexEnd(), comm_, false);
      filePrint(stderr, " ... Writing orthonormal basis to file %s ...\n", fileName.c_str());
      for (int iVec = 0; iVec < podVectorCount; ++iVec) {
        converter.paddedNodeDof6(normalizedBasis[iVec], outputBuffer);
        outputOrthoNormalFile.stateAdd(outputBuffer, solver.singularValue(iVec));
      }
    }
  }
}

} /* end namespace Rom */

extern Communicator *structCom;

Rom::DriverInterface *distrBasisOrthoDriverNew(Domain *domain) {
  return new Rom::DistrBasisOrthoDriver(domain, structCom);
}
