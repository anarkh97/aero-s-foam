#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.C>

#define ELASPLASKINHARDMAT_INSTANTIATION_HELPER(e) \
template \
ElasPlasKinHardMat<e>::ElasPlasKinHardMat(StructProp *p); \
 \
template \
void \
ElasPlasKinHardMat<e>::getStress(Tensor *stress, Tensor &strain, double *state, double temp); \
 \
template \
void  \
ElasPlasKinHardMat<e>::getElasticity(Tensor *_tm); \
 \
template \
void  \
ElasPlasKinHardMat<e>::getTangentMaterial(Tensor *tm, Tensor &strain, double *state, double temp); \
 \
template \
void  \
ElasPlasKinHardMat<e>::getStressAndTangentMaterial(Tensor *stess, Tensor *tm, Tensor &strain, double *state, double temp); \
 \
template \
void  \
ElasPlasKinHardMat<e>::updateStates(Tensor &en, Tensor &enp, double *state, double temp); \
 \
template \
void \
ElasPlasKinHardMat<e>::integrate(Tensor *_stress, Tensor *_tm, Tensor &_en, Tensor  &_enp, \
                                 double *staten, double *statenp, double temp); \
 \
template \
void \
ElasPlasKinHardMat<e>::integrate(Tensor *_stress, Tensor &_en, Tensor  &_enp, \
                                 double *staten, double *statenp, double temp); \
 \
template \
void  \
ElasPlasKinHardMat<e>::initStates(double *st); \
 \
template \
StrainEvaluator * \
ElasPlasKinHardMat<e>::getStrainEvaluator(); \
 \
template \
bool \
ElasPlasKinHardMat<e>::getBackStress(double *statenp, Tensor *_backstress); \
 \
template \
bool \
ElasPlasKinHardMat<e>::getPlasticStrain(double *statenp, Tensor *_plasticstrain); \
 \
template \
double \
ElasPlasKinHardMat<e>::getStrainEnergyDensity(Tensor &_enp, double *statenp, double temp); \
 \
template \
double \
ElasPlasKinHardMat<e>::getDissipatedEnergy(double *statenp);

ELASPLASKINHARDMAT_INSTANTIATION_HELPER(0);
ELASPLASKINHARDMAT_INSTANTIATION_HELPER(1);
ELASPLASKINHARDMAT_INSTANTIATION_HELPER(2);

