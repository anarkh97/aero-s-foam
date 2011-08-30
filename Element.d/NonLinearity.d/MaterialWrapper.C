#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>

//#define PJSA_CONVERT_TO_PK2
//#define PJSA_CONVERT_TO_M

template<typename Material>
MaterialWrapper<Material>::MaterialWrapper(Material *_mat)
{
  mat = _mat->Clone();
}

template<typename Material>
int
MaterialWrapper<Material>::getNumStates()
{
  return 0;
}

template<>
inline int
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getNumStates()
{
  return 19;
}

template<typename Material>
void
MaterialWrapper<Material>::getStress(Tensor *_stress, Tensor &_strain, double*)
{

}

template<typename Material>
void 
MaterialWrapper<Material>::getTangentMaterial(Tensor *_tm, Tensor &, double*)
{

}

template<typename Material>
void 
MaterialWrapper<Material>::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*)
{

}

template<typename Material>
void 
MaterialWrapper<Material>::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                                     double *staten, double *statenp, double)
{
  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor
  std::vector<double> ltangents; // first elasticity tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = enp[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, &ltangents);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];
}

template<>
inline void 
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                                                                    double *staten, double *statenp, double)
{
  // clone material for thread-safety reasons
  IsotropicLinearElasticJ2PlasticMaterial *clone = mat->Clone();

  std::vector<double> PlasticStrain;
  std::vector<double> BackStress;
  if(staten) {
    // set the internal variables
    // these should be the converged values at the beginning of the time step 
    for (int i = 0; i < 9; ++i) PlasticStrain.push_back(staten[i]);
    clone->SetMaterialPlasticStrain(PlasticStrain);
    for (int i = 0; i < 9; ++i) BackStress.push_back(staten[9+i]);
    clone->SetMaterialBackStress(BackStress);
    clone->SetMaterialEquivalentPlasticStrain(staten[18]);
  }

  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // Cauchy stress tensor
  std::vector<double> ltangents; // consistent tangent modulus

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = enp[3*i+j];

  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, &ltangents);

#ifndef CONVERT_TO_PK1
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];
#else
  Tensor_d0s2 sigma;
  // copy lstress into sigma (this is the Cauchy stress tensor)
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      sigma[3*i+j] = lstress[3*i+j];

  // compute the first Piola-Kirchhoff stress tensor P = J*sigma*F^{-T}
  double J;
  Tensor_d0s2 Ft, Ftinv;
  enp.getDeterminant(J);
  enp.getTranspose(Ft);
  Ft.getInverse(Ftinv);
  (*stress) = (sigma|Ftinv);
  (*stress) = J*(*stress);
#endif

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];

  if(statenp) {
    // get the updated internal variables
    PlasticStrain = clone->GetMaterialPlasticStrain();
    for (int i = 0; i < 9; ++i) statenp[i] = PlasticStrain[i];
    BackStress = clone->GetMaterialBackStress();
    for (int i = 0; i < 9; ++i) statenp[9+i] = BackStress[i];
    statenp[18] = clone->GetMaterialEquivalentPlasticStrain();
  }

  delete clone;
}

template<typename Material>
void
MaterialWrapper<Material>::initStates(double *state)
{
  for(int i = 0; i < getNumStates(); ++i) state[i] = 0;
}

template<typename Material>
double 
MaterialWrapper<Material>::getDensity()
{ 
  return mat->GetDensityInReference();
} 

template<>
inline double
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getDensity()
{
  std::cerr << "Warning IsotropicLinearElasticJ2PlasticMaterial doesn't implement GetDensityInReference function\n";
  return 0.0;
}

extern DeformationGradient deformationGradient;

template<typename Material>
StrainEvaluator *
MaterialWrapper<Material>::getStrainEvaluator() 
{
  return &deformationGradient;
}

template<typename Material>
double
MaterialWrapper<Material>::getEquivPlasticStrain(double *statenp)
{ 
  return 0;
}

template<>
inline
double MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getEquivPlasticStrain(double *statenp)
{
  return statenp[18];
}

