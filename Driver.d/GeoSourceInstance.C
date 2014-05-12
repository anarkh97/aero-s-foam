#include <Driver.d/GeoSource.C>

template
void GeoSource::distributeBCs<double>(GenSubDomain<double>*&, int*, int*);
template
void GeoSource::distributeBCs<complex<double> >(GenSubDomain<complex<double> >*&, int*, int*);

#ifdef DISTRIBUTED
template
void GeoSource::distributeOutputNodes<double>(GenSubDomain<double>*&, int*, int*);
template
void GeoSource::distributeOutputNodes<complex<double> >(GenSubDomain<complex<double> >*&, int*, int*);

template
void GeoSource::distributeOutputNodesX<double>(GenSubDomain<double>*, Connectivity*);
template
void GeoSource::distributeOutputNodesX<complex<double> >(GenSubDomain<complex<double> >*, Connectivity*);
#endif

template
GenSubDomain<double> * GeoSource::readDistributedInputFiles<double>(int, int);
template
GenSubDomain<complex<double> > * GeoSource::readDistributedInputFiles<complex<double> >(int, int);

#ifdef USE_EIGEN3
template
void GeoSource::outputEigenScalars<double>(int, Eigen::Matrix<double, Eigen::Dynamic, 1>*, double);
template
void GeoSource::outputEigenScalars<complex<double> >(int, Eigen::Matrix<complex<double>, Eigen::Dynamic, 1>*, double);

template
void GeoSource::outputEigenVectors<double>(int, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>*, double);
template
void GeoSource::outputEigenVectors<complex<double> >(int, Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic>*, double);
#endif

#define GEOSOURCE_INSTANTIATION_HELPER(dim) \
template void GeoSource::outputNodeVectors<dim>(int, double (*)[dim], int, double);\
template void GeoSource::outputNodeVectors<dim>(int, DComplex (*)[dim], int, double);\
template void GeoSource::outputNodeVectors6<dim>(int, double (*)[dim], int, double);\
template void GeoSource::outputNodeVectors6<dim>(int, DComplex (*)[dim], int, double);\
template void GeoSource::outputNodeVectors9<dim>(int, double (*)[dim], int, double);\
template void GeoSource::outputNodeVectors4<dim>(int, double (*)[dim], int, double);

GEOSOURCE_INSTANTIATION_HELPER(3);
GEOSOURCE_INSTANTIATION_HELPER(4);
GEOSOURCE_INSTANTIATION_HELPER(6);
GEOSOURCE_INSTANTIATION_HELPER(9);
GEOSOURCE_INSTANTIATION_HELPER(11);
