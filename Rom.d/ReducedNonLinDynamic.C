#include "ReducedNonLinDynamic.h"

#include "FileNameInfo.h"
#include "BasisFileStream.h"
#include "NodeDof6Buffer.h"
#include "VecNodeDof6Conversion.h"

#include "VecBasis.h"
#include "VecBasisFile.h"

#include "BasisOps.h"
#include "CholeskyUtils.h"

#include <Driver.d/Domain.h>
#include <Driver.d/GeoSource.h>
#include <Utils.d/SolverInfo.h>

#include <Corotational.d/GeomState.h>
#include <Math.d/Vector.h>

#include <Corotational.d/Corotator.h>
#include <Element.d/Element.h>

#include <Math.d/FullSquareMatrix.h>
#include <Math.d/SparseMatrix.h>

#include <Utils.d/DistHelper.h>
#include <Utils.d/dofset.h>
#include <Utils.d/vector_to_array.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <cstddef>

extern GeoSource *geoSource;
extern int verboseFlag;

namespace Rom {

namespace { // anonymous

void
make_reduced_matrix(const VecBasis &basis, const SparseMatrix &fullMatrix, FullSquareMatrix &reducedMatrix) {
  const int vectorCount = basis.vectorCount();
  reducedMatrix.setSize(vectorCount);

  Vector temp(basis.vectorInfo());
  for (int iRow = 0; iRow < vectorCount; ++iRow) {
    const_cast<SparseMatrix &>(fullMatrix).transposeMult(basis[iRow], temp);
    Vector result(reducedMatrix[iRow], vectorCount, false); // View matrix row as vector
    reduce(basis, temp, result);
  }

  reducedMatrix.symmetrize();
}

} // end anonymous namespace

// Public member functions

ReducedNonLinDynamic::ReducedNonLinDynamic(Domain *domain) :
  domain_(domain),
  timeStepCount_(0),
  fullMatrix_(NULL),
  reducedMass_(new FullSquareMatrix),
  reducedDamping_(new FullSquareMatrix),
  reducedJacobian_(new FullSquareMatrix),
  reducedSolver_(*reducedJacobian_),
  firstResidualNorm_(0.0),
  firstIncrementNorm_(0.0),
  firstEnergyNorm_(0.0),
  secondResidualNorm_(0.0),
  zeroRotMetric_(new FullSquareMatrix)
{}

ReducedNonLinDynamic::~ReducedNonLinDynamic() {
  for (std::vector<Corotator *>::const_iterator it = corotators_.begin(),
                                                it_end = corotators_.end();
                                                it != it_end; ++it) {
    const Corotator *c = *it;
    if (!dynamic_cast<const Element *>(c)) {
      delete c;
    }
  }

  delete[] melArray_;
  delete[] kelArray_;
}

void
ReducedNonLinDynamic::preProcess() {
  // Makes renumbering, connectivities and dofsets
  domain_->preProcessing();

  const int unconstrainedDofCount = domain_->numdof();
  {
    std::vector<int> bc(unconstrainedDofCount);
    std::vector<double> bcx(unconstrainedDofCount);

    // Make the boundary conditions info
    domain_->make_bc(toArray(bc), toArray(bcx));
    domain_->make_constrainedDSA(toArray(bc));
  }
  domain_->makeAllDOFs();
  
  // vcx_ stores the prescribed velocities, which are associated with prescribed displacements.
  // If a user defined displacement is used, then vcx_ will contain the user defined velocities also.
  vcx_.resize(unconstrainedDofCount);
  std::fill_n(vcx_.begin(), vcx_.size(), 0.0);
/* PJSA 3/5/2012: why should the initial velocities be in vcx??
  {
    const BCond* velBc = domain_->getInitVelocity();
    for (int iBc = 0; iBc < domain_->numInitVelocity(); ++iBc) {
      const int dof = domain_->getCDSA()->locate(velBc[iBc].nnum, 1 << velBc[iBc].dofnum);
      if (dof >= 0) {
        vcx_[dof] = velBc[iBc].val;
      }
    }
  }
*/

  fullMatrix_.reset(domain_->constructDBSparseMatrix<double>());
  std::auto_ptr<SparseMatrix> fullDampingMatrix;
  if (hasDamping()) {
    fullDampingMatrix.reset(domain_->constructDBSparseMatrix<double>());
  }

  AllOps<double> operators;
  operators.M = fullMatrix_.get();
  operators.C = fullDampingMatrix.get();

  domain_->createKelArray(kelArray_, melArray_); // Also zero out the rotational dofs of the mass matrix
  domain_->makeSparseOps<double>(operators, 0.0, 0.0, 0.0, NULL, kelArray_, melArray_);

  // Read reduced basis
  {
    FileNameInfo fileInfo;
    VecNodeDof6Conversion converter(*domain_->getCDSA());
    BasisInputStream reducedBasisInput(BasisFileId(fileInfo, BasisId::STATE, BasisId::POD), converter);

    if (reducedBasisInput.vectorSize() != fullSolVecInfo()) {
      throw std::domain_error("Projection basis has incorrect #rows");
    }

    const int userSubspaceSize = domain_->solInfo().maxSizePodRom;
    const int reducedSubspaceSize = userSubspaceSize ?
      std::min(userSubspaceSize, reducedBasisInput.size()) :
      reducedBasisInput.size();

    readVectors(reducedBasisInput, reducedBasis_, reducedSubspaceSize);

    filePrint(stderr, "Projection subspace of dimension = %d\n", reducedBasis_.vectorCount());
  }

  {
    zeroRotMetric_->setSize(solVecInfo());
    Vector temp(solVecInfo());
    for (int iRow = 0; iRow < solVecInfo(); ++iRow) {
      temp = reducedBasis_[iRow];
      zeroRotDofs(temp);
      Vector result((*zeroRotMetric_)[iRow], solVecInfo(), false); // View matrix row as vector
      reduce(reducedBasis_, temp, result);
    }
    zeroRotMetric_->symmetrize();
    cholesky_factor_upper(*zeroRotMetric_);
  }

  // Assemble reduced matrices
  make_reduced_matrix(reducedBasis_, *fullMatrix_, *reducedMass_);
  if (hasDamping()) {
    make_reduced_matrix(reducedBasis_, *fullDampingMatrix, *reducedDamping_);
  }

  corotators_.resize(domain_->numElements());
  domain_->createCorotators(toArray(corotators_));

  if (domain->solInfo().iacc_switch) {
    // Build solver for M^{-1}
    reducedJacobian_->copy(*reducedMass_);
    cholesky_factor_upper(*reducedJacobian_);
  }
}

void
ReducedNonLinDynamic::computeTimeInfo() {
  // Time integration information
  const double finalTime = domain_->solInfo().tmax;
  const double dt = getDt();

  timeStepCount_ = static_cast<int>((finalTime + 0.49 * dt) / dt);

  // Compute time remainder
  const double remainder = finalTime - timeStepCount_ * dt;
  if (std::abs(remainder) > 0.01 * dt) {
    const double newFinalTime = domain_->solInfo().tmax = timeStepCount_ * dt;
    filePrint(stderr, " Warning: Total time is being changed to : %e\n", newFinalTime);
  }
}

int
ReducedNonLinDynamic::getInitState(Vector &d_n, Vector &v_n, Vector &a_n, Vector &v_p) const {
  VecBasis fullVectors(4, fullSolVecInfo());
  
  domain_->initDispVeloc(fullVectors[0], fullVectors[1], fullVectors[2], fullVectors[3]);
 
  if (domain_->solInfo().zeroRot) {
    zeroRotDofs(fullVectors[1]);
    zeroRotDofs(fullVectors[2]);
    zeroRotDofs(fullVectors[3]);
  }

  reduce(reducedBasis_, fullVectors[0], d_n);
  reduce(reducedBasis_, fullVectors[1], v_n);
  reduce(reducedBasis_, fullVectors[2], a_n);
  reduce(reducedBasis_, fullVectors[3], v_p);

  return getAeroAlg();
}

void
ReducedNonLinDynamic::getInitialTime(int &initTimeIndex, double &initTime) const {
  initTimeIndex = domain_->solInfo().initialTimeIndex;
  initTime = domain_->solInfo().initialTime;
}

void
ReducedNonLinDynamic::readRestartFile(Vector &, Vector &, Vector &, Vector &, GeomState &) {
  // Not supported
}

int
ReducedNonLinDynamic::solVecInfo() const {
  return reducedBasis_.vectorCount();
}

int
ReducedNonLinDynamic::sysVecInfo() const {
  return solVecInfo();
}

int
ReducedNonLinDynamic::elemVecInfo() const {
  return domain_->maxNumDOF();
}

double
ReducedNonLinDynamic::getDt() const {
  return domain_->solInfo().getTimeStep();
}

double
ReducedNonLinDynamic::getDelta() const {
  return getDt() / 2.0;
}

int
ReducedNonLinDynamic::getMaxStep() const {
  return timeStepCount_;
}

double
ReducedNonLinDynamic::getTolerance() const {
  return domain_->solInfo().getNLInfo().tolRes * firstResidualNorm_;
}

int
ReducedNonLinDynamic::getMaxit() const {
  return domain_->solInfo().getNLInfo().maxiter;
}

void
ReducedNonLinDynamic::getNewmarkParameters(double &beta, double &gamma, double &alpha_f, double &alpha_m) const {
  beta  = domain_->solInfo().newmarkBeta;
  gamma = domain_->solInfo().newmarkGamma;
  alpha_f = domain_->solInfo().newmarkAlphaF;
  alpha_m = domain_->solInfo().newmarkAlphaM;
}

int
ReducedNonLinDynamic::getAeroAlg() const {
  return -1; // Disabled
}

int
ReducedNonLinDynamic::getThermoeFlag() const {
  return -1; // Disabled
}

int
ReducedNonLinDynamic::getThermohFlag() const {
  return -1; // Disabled
}

int
ReducedNonLinDynamic::getAeroheatFlag() const {
  return -1; // Disabled
}

ReducedNonLinDynamic::Solver *
ReducedNonLinDynamic::getSolver() {
  return const_cast<Solver *>(const_cast<const ReducedNonLinDynamic *>(this)->getSolver());
}

const ReducedNonLinDynamic::Solver *
ReducedNonLinDynamic::getSolver() const {
  return &reducedSolver_;
}

GeomState *
ReducedNonLinDynamic::createGeomState() const {
  return new GeomState(*domain_->getDSA(), *domain_->getCDSA(), domain_->getNodes(), &domain_->getElementSet());
}

GeomState *
ReducedNonLinDynamic::copyGeomState(const GeomState* geomState) const {
  return new GeomState(*geomState);
}

void
ReducedNonLinDynamic::updateStates(GeomState *refState, GeomState &geomState) {
  domain_->updateStates(refState, geomState, toArray(corotators_));
}

void
ReducedNonLinDynamic::updatePrescribedDisplacement(GeomState *geomState) const {
  // Initial conditions
  if (domain_->numInitDisp6() > 0) {
    geomState->updatePrescribedDisplacement(domain_->getInitDisp6(), domain_->numInitDisp6(), domain_->getNodes());
  } else if (domain_->numInitDisp() > 0) {
    geomState->updatePrescribedDisplacement(domain_->getInitDisp(), domain_->numInitDisp(), domain_->getNodes());
  }
  
  // Boundary conditions
  if (domain_->nDirichlet() > 0) {
    geomState->updatePrescribedDisplacement(domain_->getDBC(), domain_->nDirichlet(), domain_->getNodes());
  }
}

void
ReducedNonLinDynamic::getConstForce(Vector &constantForce) const {
  Vector fullConstantForce(fullSolVecInfo());
  domain_->computeConstantForce(fullConstantForce);
  reduce(reducedBasis_, fullConstantForce, constantForce);
}

void
ReducedNonLinDynamic::getExternalForce(Vector &externalForce, const Vector &constantForce, int /*timeIndex*/, double time,
                                       const GeomState * /*geomState*/, Vector & /*elemNonConForce*/, const Vector & /*aeroForce*/, double /*delta*/) const {
  externalForce = constantForce;

  // Add time-dependent component of external forces to time-independent (constant) component
  Vector fullExternalForce(fullSolVecInfo());
  domain_->computeExtForce(fullExternalForce, time);
  reduceAdd(reducedBasis_, fullExternalForce, externalForce);
}

double
ReducedNonLinDynamic::getStiffAndForce(GeomState &geomState, Vector &residual, Vector &elementInternalForce, double time, GeomState *refState) {
  getStiffAndReducedForceFromDomain(residual, geomState, time, refState, elementInternalForce);
  return residual.norm();
}

void
ReducedNonLinDynamic::reBuild(GeomState &geomState, int iteration, double delta) {
  // Rebuild every updateK iterations
  if (iteration % domain_->solInfo().getNLInfo().updateK == 0) {
    // Recompute tangent stiffness matrix
    {
      fullMatrix_->zeroAll();

      AllOps<double> operators;
      operators.K = fullMatrix_.get();
      domain_->makeSparseOps<double>(operators, 0.0, 0.0, 0.0, NULL, kelArray_, NULL);
    }

    // Coefficients for the time-integrator
    double beta, gamma, alpha_f, alpha_m;
    const double dt = getDt();
    getNewmarkParameters(beta, gamma, alpha_f, alpha_m);
    const double Kcoef = dt * dt * beta;
    const double Ccoef = dt * gamma;
    const double Mcoef = (1.0 - alpha_m) / (1.0 - alpha_f);

    // Jacobian_r <- coef_M * M_r + coef_C * C_r + coef_K * K_r
    make_reduced_matrix(reducedBasis_, *fullMatrix_, *reducedJacobian_);
    *reducedJacobian_ *= Kcoef;

    if (hasDamping()) {
      reducedJacobian_->linAdd(Ccoef, *reducedDamping_);
    }

    reducedJacobian_->linAdd(Mcoef, *reducedMass_);

    // Factor reduced system
    cholesky_factor_upper(*reducedJacobian_);
  }
}

void
ReducedNonLinDynamic::formRHSinitializer(Vector &, Vector &, Vector &, GeomState &, Vector &) {
  throw std::logic_error("Not implemented");
}

void
ReducedNonLinDynamic::formRHSinitializer(Vector &externalForce, Vector &velocity, Vector &elementInternalForce, GeomState &geomState, Vector &rhs, GeomState *refState) {
  // rhs <- externalForce - C * velocity - internalForce
  rhs = externalForce;
  
  if (hasDamping()) {
    Vector temp(solVecInfo());
    reducedDamping_->mult(velocity, temp);
    rhs -= temp;
  }
 
  getStiffAndReducedForceFromDomain(rhs, geomState, domain_->solInfo().initialTime, refState, elementInternalForce);
}

double
ReducedNonLinDynamic::formRHScorrector(Vector &inc_displacement, Vector &velocity, Vector &acceleration, Vector &residual, Vector &rhs) {
  double beta, gamma, alpha_f, alpha_m;
  getNewmarkParameters(beta, gamma, alpha_f, alpha_m);
  const double dt = getDt();
  
  // rhs = dt*dt*beta*residual - ((1-alpha_m)/(1-alpha_f)*M+dt*gamma*C)*inc_displacement
  //       + (dt*(1-alpha_m)*M - dt*dt*(beta-(1-alpha_f)*gamma)*C)*velocity
  //       + (dt*dt*((1-alpha_m)/2-beta)*M - dt*dt*dt*(1-alpha_f)*(2*beta-gamma)/2*C)*acceleration
  Vector temp(solVecInfo());
  temp.linC(-(1-alpha_m)/(1-alpha_f), inc_displacement, dt*(1-alpha_m), velocity, dt*dt*((1-alpha_m)/2-beta), acceleration);
  reducedMass_->mult(temp, rhs);
  if (hasDamping()) {
    temp.linC(-dt*gamma, inc_displacement, -dt*dt*(beta-(1-alpha_f)*gamma), velocity, -dt*dt*dt*(1-alpha_f)*(2*beta-gamma)/2, acceleration);
    reducedDamping_->multAdd(temp, rhs);
  }
  rhs.linAdd(dt*dt*beta, residual);
  
  return rhs.norm();
}

int
ReducedNonLinDynamic::checkConvergence(int iteration, double residualNorm, Vector &residual, Vector &dv, double time) {
  const double incrementNorm = dv.norm();
  const double energyNorm = residual * dv;

  if (iteration == 0) {
    firstResidualNorm_ = residualNorm;
    firstIncrementNorm_  = incrementNorm;
    firstEnergyNorm_ = energyNorm;
  } else if (iteration == 1) {
    secondResidualNorm_ = residualNorm;
  }

  const double relativeResidualNorm = residualNorm / firstResidualNorm_;
  const double relativeIncrementNorm  = incrementNorm / firstIncrementNorm_;
  const double relativeEnergyNorm = energyNorm / firstEnergyNorm_;

  int converged = 0;

  if (residualNorm <= getTolerance() && incrementNorm <= getIncrementTolerance()) { 
    converged = 1;
  }

  // Check for divergence
  const double divergenceGrowthFactor = 1.0e10;
  if (residualNorm >= divergenceGrowthFactor * firstResidualNorm_ && residualNorm > secondResidualNorm_) {
    converged = -1;
  }

  if (verboseFlag) {
    filePrint(stderr, " Iteration # %d\n", iteration);
    filePrint(stderr, " r      = %e dv      = %e energy      = %e\n"
                      " rel. r = %e rel. dv = %e rel. energy = %e\n",
                      residualNorm, incrementNorm, energyNorm,
                      relativeResidualNorm, relativeIncrementNorm, relativeEnergyNorm);
  }

  return converged;
}


void
ReducedNonLinDynamic::dynamOutput(GeomState *geomState, Vector &velocity, Vector & /*velocity_p*/,
                                  double time, int step,
                                  Vector &externalForce, Vector &aeroForce, Vector &acceleration,
                                  GeomState *refState) {
  Vector fullVelocity(fullSolVecInfo());
  Vector fullAcceleration(fullSolVecInfo());
  Vector fullExternalForce(fullSolVecInfo());
  Vector fullAeroForce(fullSolVecInfo(), 0.0); // Disabled
 
  expand(reducedBasis_, velocity, fullVelocity);
  expand(reducedBasis_, acceleration, fullAcceleration);
  expand(reducedBasis_, externalForce, fullExternalForce);

  domain_->postProcessing(geomState, fullExternalForce, fullAeroForce, time, step + 1,
                          fullVelocity.data(), toArray(vcx_),
                          toArray(corotators_), melArray_, fullAcceleration.data(),
                          NULL, refState, NULL);
}

void
ReducedNonLinDynamic::processLastOutput() {
  domain_->solInfo().lastIt = true;
  
  OutputInfo *oinfo = geoSource->getOutputInfo();
  const int outputInfoCount = geoSource->getNumOutInfo();

  for (int iOut = 0; iOut < outputInfoCount; ++iOut) {
    oinfo[iOut].interval = 1;
  }
}

void 
ReducedNonLinDynamic::dynamCommToFluid(GeomState *, GeomState *, Vector &, Vector &,
                                       Vector &, Vector &, int, int, int) {
  throw std::logic_error("Not implemented");
}

void
ReducedNonLinDynamic::printTimers(double) {
  // Not implemented
}


// Private member functions

int
ReducedNonLinDynamic::fullSolVecInfo() const {
  return domain_->numUncon();
}

bool
ReducedNonLinDynamic::hasDamping() const {
  return (domain_->solInfo().alphaDamp != 0.0) || (domain_->solInfo().betaDamp != 0.0);
}

void
ReducedNonLinDynamic::getStiffAndReducedForceFromDomain(Vector &force,
                                                        GeomState &geomState, double time, GeomState *refState,
                                                        Vector &elementInternalForce) {
  Vector fullForce(fullSolVecInfo(), 0.0);
  domain_->getStiffAndForce(geomState, elementInternalForce, toArray(corotators_), kelArray_, fullForce, 1.0, time, refState);
  reduceAdd(reducedBasis_, fullForce, force);
}

double
ReducedNonLinDynamic::getIncrementTolerance() const {
  return domain_->solInfo().getNLInfo().tolInc * firstIncrementNorm_;
}

void
ReducedNonLinDynamic::zeroRotDofs(Vector &vec) const {
  ::zeroRotDofs(*domain_->getCDSA(), vec);
}

// Solver class

ReducedNonLinDynamic::Solver::Solver(const GenFullSquareMatrix<double> &matrix) :
  matrix_(matrix)
{}

void
ReducedNonLinDynamic::Solver::reSolve(Vector &v) {
  cholesky_solve_upper(matrix_, v.data());
}

// Updater class

ReducedNonLinDynamic::Updater::RefState *
ReducedNonLinDynamic::Updater::initRef(const GeomType *u) {
  return new RefState(*u);
}

ReducedNonLinDynamic::Updater::StateIncr *
ReducedNonLinDynamic::Updater::initInc(const GeomType *, const VecType *v) {
  return new VecType(*v);
}

void
ReducedNonLinDynamic::Updater::copyState(const GeomType *gn, RefState *gp) {
  *gp = *gn;
}

void
ReducedNonLinDynamic::Updater::zeroInc(StateIncr *du) {
  du->zero();
}

void
ReducedNonLinDynamic::Updater::get_inc_displacement(ProbDescr *pbd,
                                                    GeomType *geomState, StateIncr &du, GeomType *refState,
                                                    bool zeroRot) {
  VecType fullDu(pbd->fullSolVecInfo());
  geomState->get_inc_displacement(fullDu, *refState, zeroRot);
  reduce(pbd->reducedBasis_, fullDu, du);
  if (zeroRot) {
    cholesky_solve_upper(*pbd->zeroRotMetric_, du.data());
  }
}

double
ReducedNonLinDynamic::Updater::integrate(ProbDescr *pbd, RefState *refState, GeomType *geomState,
	                     	                 StateIncr *du, VecType &residual, 
                                         VecType &elementInternalForce, VecType &, VecType &,
                                         VecType &, double midTime) {
  {
    VecType fullDu(pbd->fullSolVecInfo());
    expand(pbd->reducedBasis_, *du, fullDu);
    geomState->update(fullDu);
  }

  return pbd->getStiffAndForce(*geomState, residual, elementInternalForce, midTime, refState);
}

void
ReducedNonLinDynamic::Updater::midpointIntegrate(ProbDescr *pbd,
                                                 VecType &velN, double delta,
                                                 GeomType *refState, GeomType *geomState,
                                                 StateIncr *, VecType &, 
                                                 VecType &, VecType &, VecType &acceleration,
                                                 bool zeroRot) {
  VecType fullVelocity(pbd->fullSolVecInfo()), fullAcceleration(pbd->fullSolVecInfo());
  expand(pbd->reducedBasis_, velN, fullVelocity);
  expand(pbd->reducedBasis_, acceleration, fullAcceleration);

  {
    double beta, gamma, alpha_f, alpha_m;
    pbd->getNewmarkParameters(beta, gamma, alpha_f, alpha_m);

    geomState->midpoint_step_update(fullVelocity, fullAcceleration, delta, *refState, beta, gamma, alpha_f, alpha_m, zeroRot);
  }
  
  if (zeroRot) {
    geomState->zeroRotDofs(fullVelocity);
    geomState->zeroRotDofs(fullAcceleration);
  }

  reduce(pbd->reducedBasis_, fullVelocity, velN);
  reduce(pbd->reducedBasis_, fullAcceleration, acceleration);

  if (zeroRot) {
    cholesky_solve_upper(*pbd->zeroRotMetric_, velN.data());
    cholesky_solve_upper(*pbd->zeroRotMetric_, acceleration.data());
  }
}

void
ReducedNonLinDynamic::Updater::updateIncr(StateIncr *du, VecType &ddu) {
  *du = ddu;
}

double
ReducedNonLinDynamic::Updater::formRHScorrector(ProbDescr *pbd, VecType &inc_displac, VecType &vel_n,
                                                VecType &accel, VecType &residual, VecType &rhs, GeomType *, double delta) {
  return pbd->formRHScorrector(inc_displac, vel_n, accel, residual, rhs); 
}

void
ReducedNonLinDynamic::Updater::copyTo(RefState *, GeomType *, GeomType *, StateIncr *, VecType &, VecType &, VecType &, VecType &,
                                      RefState *, GeomType *, GeomType *, StateIncr *, VecType &, VecType &, VecType &, VecType &) {
  throw std::logic_error("Not implemented");
}

} // end namespace Rom
