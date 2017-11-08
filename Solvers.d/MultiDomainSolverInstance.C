#include <Solvers.d/MultiDomainSolver.C>
#include <complex>

#define MULTIDOMAINSOLVER_INSTANTIATION_HELPER(Scalar) \
template \
void \
MultiDomainSolver<Scalar>::reSolve(GenDistrVector<Scalar>&);\
\
template \
void \
MultiDomainSolver<Scalar>::solve(GenDistrVector<Scalar>&, GenDistrVector<Scalar>&);\
\
template \
double \
MultiDomainSolver<Scalar>::getFNormSq(GenDistrVector<Scalar>&);\
\
template \
void \
MultiDomainSolver<Scalar>::multLTinv(Scalar*, GenDistrVector<Scalar>&);\
\
template \
void \
MultiDomainSolver<Scalar>::multL(Scalar*, GenDistrVector<Scalar>&);\
\
template \
void \
MultiDomainSolver<Scalar>::multLinv(GenDistrVector<Scalar>&, Scalar*);\
\
template \
void \
MultiDomainSolver<Scalar>::multLT(GenDistrVector<Scalar>&, Scalar*);

MULTIDOMAINSOLVER_INSTANTIATION_HELPER(double);
MULTIDOMAINSOLVER_INSTANTIATION_HELPER(std::complex<double>);

