#include "GappyNonLinDynamic.h"

#include "BasisFileStream.h"
#include "VecBasisFile.h"
#include "FileNameInfo.h"

#include <Driver.d/GeoSource.h> 

#include <stdexcept>

extern GeoSource *geoSource;

GappyNonLinDynamic::GappyNonLinDynamic(Domain * d) :
  NonLinDynamic(d)
{}

void
GappyNonLinDynamic::fillBasisFromInput(const std::string &fileName, VecBasis &target) const {
  BasisInputStream input(fileName, *vecNodeDof6Conversion_);
  input >> target;
}

void
GappyNonLinDynamic::fillRestrictedBasisFromInput(const std::string &fileName, VecBasis &target) const {
  BasisInputStream input(fileName, *vecNodeDof6Conversion_);
  readRestrictedVectors(input, target, *restrictionMapping_);
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

  {
    FileNameInfo fileInfo;

    fillBasisFromInput(fileInfo.fileName(BasisId(BasisId::STATE, BasisId::GAPPY_POD)), reducedBasis_);

    fillRestrictedBasisFromInput(fileInfo.fileName(BasisId(BasisId::RESIDUAL, BasisId::GAPPY_POD)), residualProjection_);
    fillRestrictedBasisFromInput(fileInfo.fileName(BasisId(BasisId::JACOBIAN, BasisId::GAPPY_POD)), jacobianProjection_);
  }

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

double
GappyNonLinDynamic::getResidualNorm(const Vector &residual) {
  return getSolver()->projectAndComputeNorm(residual);
}

bool
GappyNonLinDynamic::factorWhenBuilding() const {
  return false; // Delayed factorization
}
 
