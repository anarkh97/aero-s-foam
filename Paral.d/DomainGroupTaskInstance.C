#include <Paral.d/DomainGroupTask.C>

#define DGTASK_INSTANTIATION_HELPER(Scalar) \
template \
GenDomainGroupTask<Scalar>::GenDomainGroupTask(int _nsub, GenSubDomain<Scalar> **_sd, double _cm,  \
                                               double _cc, double _ck, Rbm **_rbms, FullSquareMatrix **_kelArray, \
                                               double _alpha, double _beta, int _numSommer, int _solvertype, \
                                               FSCommunicator *_com, FullSquareMatrix **_melArray, \
                                               FullSquareMatrix **_celArray, bool elemsetHasDamping); \
template \
GenDomainGroupTask<Scalar>::~GenDomainGroupTask(); \
template \
void \
GenDomainGroupTask<Scalar>::runFor(int isub, bool make_feti); \

DGTASK_INSTANTIATION_HELPER(double);
DGTASK_INSTANTIATION_HELPER(complex<double>);
