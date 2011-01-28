#include "BasisOrthoDriver.h"

#include "BasisOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "BasisFileStream.h"
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
  BasisInputStream input("RawBasis", converter); //TODO filename
  BasisOutputStream output("OrthoBasis", converter); //TODO filename

  BasisOrthogonalization solver;
  solver.basisNew(input, output);
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

