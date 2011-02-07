#include "BasisOrthoDriver.h"

#include "BasisOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "BasisFileStream.h"
#include "FileNameInfo.h"
#include "SimpleBuffer.h"

#include <Driver.d/Domain.h>
#include <Utils.d/dofset.h>

BasisOrthoDriver::BasisOrthoDriver(Domain *domain) :
  domain_(domain)
{}

void
BasisOrthoDriver::solve() {
  preProcess();
 
  VecNodeDof6Conversion converter(*domain_->getCDSA());
  FileNameInfo fileInfo;
  BasisOrthogonalization solver;

  std::vector<BasisId::Type> workload;
  workload.push_back(BasisId::RESIDUAL);
  workload.push_back(domain_->solInfo().gaussNewtonPodRom ? BasisId::JACOBIAN : BasisId::STATE);

  for (std::vector<BasisId::Type>::const_iterator it = workload.begin(); it != workload.end(); ++it) {
    BasisId::Type type = *it;
    BasisInputStream input(fileInfo.fileName(BasisId(type, BasisId::SNAPSHOTS)), converter);
    BasisOutputStream output(fileInfo.fileName(BasisId(type, BasisId::POD)), converter);

    solver.basisNew(input, output);
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

