#include "BasisOrthoDriver.h"

#include "SvdOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "BasisFileStream.h"
#include "FileNameInfo.h"
#include "SimpleBuffer.h"

#include <Driver.d/Domain.h>
#include <Utils.d/dofset.h>

#include <utility>
#include <algorithm>

BasisOrthoDriver::BasisOrthoDriver(Domain *domain) :
  domain_(domain)
{}

void
BasisOrthoDriver::solve() {
  preProcess();
 
  VecNodeDof6Conversion converter(*domain_->getCDSA());
  FileNameInfo fileInfo;
  SvdOrthogonalization solver;

  std::vector<BasisId::Type> workload;
  if (domain_->solInfo().gaussNewtonPodRom) {
    workload.push_back(BasisId::RESIDUAL);
    workload.push_back(BasisId::JACOBIAN);
  } else {
    workload.push_back(BasisId::STATE);
  }

  for (std::vector<BasisId::Type>::const_iterator it = workload.begin(); it != workload.end(); ++it) {
    BasisId::Type type = *it;

    {
      BasisInputStream input(fileInfo.fileName(BasisId(type, BasisId::SNAPSHOTS)), converter);
      solver.matrixSizeIs(input.vectorSize(), input.size());

      int iCol = 0;
      while (input) {
        input >> solver.matrixCol(iCol++);
      }
    }

    solver.solve();

    BasisOutputStream output(fileInfo.fileName(BasisId(type, BasisId::POD)), converter);
    const int orthoBasisDim = domain->solInfo().maxSizePodRom ?
                              std::min(domain->solInfo().maxSizePodRom, solver.singularValueCount()) :
                              solver.singularValueCount();

    for (int iVec = 0; iVec < orthoBasisDim; ++iVec) {
      output << std::make_pair(solver.singularValue(iVec), solver.matrixCol(iVec));
    }
  }
}

void
BasisOrthoDriver::preProcess() {
  domain_->preProcessing();
 
  // Build the constrained DofSetArray incorporating the boundary conditions 
  const int numdof = domain_->numdof();
  SimpleBuffer<int> bc(numdof);
  SimpleBuffer<double> bcx(numdof);

  domain_->make_bc(bc.array(), bcx.array());
  domain_->make_constrainedDSA(bc.array());
}

