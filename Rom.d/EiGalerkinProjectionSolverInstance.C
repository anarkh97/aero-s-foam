#ifdef USE_EIGEN3
#include "EiGalerkinProjectionSolver.C"

namespace Rom {

template
GenEiSparseGalerkinProjectionSolver<double>
::GenEiSparseGalerkinProjectionSolver(Connectivity*, DofSetArray*, ConstrainedDSA*, bool);

template
GenEiSparseGalerkinProjectionSolver<std::complex<double> >
::GenEiSparseGalerkinProjectionSolver(Connectivity*, DofSetArray*, ConstrainedDSA*, bool);

template
void
GenEiSparseGalerkinProjectionSolver<double>
::zeroAll();

template
void
GenEiSparseGalerkinProjectionSolver<std::complex<double> >
::zeroAll();

template
void
GenEiSparseGalerkinProjectionSolver<double>
::addReducedMass(double);

template
void
GenEiSparseGalerkinProjectionSolver<std::complex<double> >
::addReducedMass(double);

template
void
GenEiSparseGalerkinProjectionSolver<double>
::projectionBasisIs(const GenVecBasis<double>&);

template
void
GenEiSparseGalerkinProjectionSolver<std::complex<double> >
::projectionBasisIs(const GenVecBasis<std::complex<double> >&);

template
void
GenEiSparseGalerkinProjectionSolver<double>
::factor();

template
void
GenEiSparseGalerkinProjectionSolver<std::complex<double> >
::factor();
   
template
void
GenEiSparseGalerkinProjectionSolver<double>
::reSolve(GenVector<double>&);

template
void
GenEiSparseGalerkinProjectionSolver<std::complex<double> >
::reSolve(GenVector<std::complex<double> >&);

template
void
GenEiSparseGalerkinProjectionSolver<double>
::solve(GenVector<double>&, GenVector<double>&);

template
void
GenEiSparseGalerkinProjectionSolver<std::complex<double> >
::solve(GenVector<std::complex<double> >&, GenVector<std::complex<double> >&);

} /* end namespace Rom */
#endif
