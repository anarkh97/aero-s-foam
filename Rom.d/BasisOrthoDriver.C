#include "BasisOrthoDriver.h"

#include "BasisOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "BasisFileIterator.h"

#include <Driver.d/Domain.h>
#include <Utils.d/dofset.h>

#include <OOPita.d/SimpleBuffer.h>

BasisOrthoDriver::BasisOrthoDriver(Domain *domain) :
  domain_(domain)
{}

void
BasisOrthoDriver::solve() {
  preProcess();
 
  VecNodeDof6Conversion converter(*domain_->getCDSA());
  BasisInputRange input("RawBasis", converter); //TODO filename
  BasisOutputRange output("OrthoBasis", converter); //TODO filename

  BasisOrthogonalization solver;
  solver.basisNew(input, output.begin());
}

void
BasisOrthoDriver::preProcess() {
  domain_->preProcessing();
 
  // Build the constrained DofSetArray incorporating the boundary conditions 
  const int numdof = domain_->numdof();
  ::Pita::SimpleBuffer<int> bc(numdof);
  ::Pita::SimpleBuffer<double> bcx(numdof);

  domain_->make_bc(bc.array(), bcx.array());
  domain_->make_constrainedDSA(bc.array());
}

