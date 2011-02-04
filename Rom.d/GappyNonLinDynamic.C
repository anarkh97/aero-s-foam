#include "GappyNonLinDynamic.h"

#include "BasisFileStream.h"
#include "VecBasisFile.h"

#include <Driver.d/GeoSource.h> 

#include <stdexcept>

extern GeoSource *geoSource;

GappyNonLinDynamic::GappyNonLinDynamic(Domain * d) :
  NonLinDynamic(d)
{}

void
GappyNonLinDynamic::fillBasisFromInput(const std::string &fileName, VecBasis &target) {
  BasisInputStream input(fileName, *vecNodeDof6Conversion_);
  input >> target;
}

void
GappyNonLinDynamic::fillRestrictedBasisFromInput(const std::string &fileName, VecBasis &target) {
  BasisInputStream input(fileName, *vecNodeDof6Conversion_);

  VecBasis result(input.size(), restrictionMapping_->restrictedInfo());
  Vector temp(restrictionMapping_->originInfo());

  for (VecBasis::iterator vecIt = result.begin(); vecIt != result.end(); ++vecIt) {
    input >> temp;
    restrictionMapping_->restriction(temp, *vecIt);
  }

  target.swap(result);
}

void
GappyNonLinDynamic::preProcess() {
  NonLinDynamic::preProcess();

  if (!dynamic_cast<GappyProjectionSolver *>(NonLinDynamic::getSolver())) {
    throw std::runtime_error("Solver must be a GappyProjectionSolver");
  }
 
  // Read system approximation parameters
  const GeoSource &gs = *geoSource;
  restrictionMapping_.reset(new NodalRestrictionMapping(*this->domain->getCDSA(),
                                                        gs.sampleNodeBegin(), gs.sampleNodeEnd()));
  
  vecNodeDof6Conversion_.reset(new VecNodeDof6Conversion(*this->domain->getCDSA()));
  
  fillBasisFromInput("GappyReducedBasis",  reducedBasis_);       // TODO filename
  fillRestrictedBasisFromInput("GappyJacobianBasis", jacobianProjection_); // TODO filename
  fillRestrictedBasisFromInput("GappyResidualBasis", residualProjection_); // TODO filename

  // Setup solver
  getSolver()->systemApproximationIs(*restrictionMapping_,
                                     reducedBasis_,
                                     jacobianProjection_,
                                     residualProjection_);
  getSolver()->factor(); // Delayed factorization
}

GappyProjectionSolver *
GappyNonLinDynamic::getSolver() {
  return static_cast<GappyProjectionSolver *>(NonLinDynamic::getSolver());
}

bool
GappyNonLinDynamic::factorWhenBuilding() const {
  return false; // Delayed factorization
}
 
