#include <Utils.d/NodeSpaceArray.h>
#include <Element.d/NonLinearity.d/StrainEvaluator.h>
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif

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
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, NULL);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];
}

template<>
inline void
MaterialWrapper<NeoHookean>::getStress(Tensor *_stress, Tensor &_strain, double* state)
{
  // Note: this function is called for post-processing.
  // In this case we prefer to output the PK2 stress to be consistent with other
  // materials, and also because the second invariant of the deviatoric PK1 stress
  // can be negative.
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, NULL);

#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix3d,Eigen::RowMajor> F(&lstrain[0]), P(&lstress[0]), s(&(*stress)[0]);
  s = F.inverse()*P; // symmetric 2nd Piola-Kirchhoff stress tensor, S = F^{-1}*P
#else
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];
#endif
}

template<>
inline void
MaterialWrapper<MooneyRivlin>::getStress(Tensor *_stress, Tensor &_strain, double* state)
{
  // Note: this function is called for post-processing.
  // In this case we prefer to output the PK2 stress to be consistent with other
  // materials, and also because the second invariant of the deviatoric PK1 stress
  // can be negative.
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, NULL);

#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix3d,Eigen::RowMajor> F(&lstrain[0]), P(&lstress[0]), s(&(*stress)[0]);
  s = F.inverse()*P; // symmetric 2nd Piola-Kirchhoff stress tensor, S = F^{-1}*P
#else
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];
#endif
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getStress(Tensor *_stress, Tensor &_strain, double* state)
{
  // clone material for thread-safety reasons
  IsotropicLinearElasticJ2PlasticMaterial *clone = mat->Clone();

  std::vector<double> PlasticStrain;
  std::vector<double> BackStress;
  if(state) {
    // set the internal variables    
    for (int i = 0; i < 9; ++i) PlasticStrain.push_back(state[i]);
    clone->SetMaterialPlasticStrain(PlasticStrain);
    for (int i = 0; i < 9; ++i) BackStress.push_back(state[9+i]);
    clone->SetMaterialBackStress(BackStress);
    clone->SetMaterialEquivalentPlasticStrain(state[18]);
  }
  
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);
    
  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // Cauchy stress tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];
  
  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, NULL);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];

  delete clone;
}

template<typename Material>
void 
MaterialWrapper<Material>::getTangentMaterial(Tensor *_tm, Tensor &_strain, double*)
{
  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor
  std::vector<double> ltangents; // first elasticity tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, &ltangents);

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getTangentMaterial(Tensor *_tm, Tensor &_strain, double* state)
{
  // clone material for thread-safety reasons
  IsotropicLinearElasticJ2PlasticMaterial *clone = mat->Clone();

  std::vector<double> PlasticStrain;
  std::vector<double> BackStress;
  if(state) {
    // set the internal variables    
    for (int i = 0; i < 9; ++i) PlasticStrain.push_back(state[i]);
    clone->SetMaterialPlasticStrain(PlasticStrain);
    for (int i = 0; i < 9; ++i) BackStress.push_back(state[9+i]);
    clone->SetMaterialBackStress(BackStress);
    clone->SetMaterialEquivalentPlasticStrain(state[18]);
  }
  
  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);
    
  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // Cauchy stress tensor
  std::vector<double> ltangents; // consistent tangent modulus

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];
  
  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, &ltangents);

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];

  delete clone;
}

template<typename Material>
void 
MaterialWrapper<Material>::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*)
{
  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor
  std::vector<double> ltangents; // first elasticity tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];

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
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double* state)
{
  // clone material for thread-safety reasons
  IsotropicLinearElasticJ2PlasticMaterial *clone = mat->Clone();

  std::vector<double> PlasticStrain;
  std::vector<double> BackStress;
  if(state) {
    // set the internal variables
    for (int i = 0; i < 9; ++i) PlasticStrain.push_back(state[i]);
    clone->SetMaterialPlasticStrain(PlasticStrain);
    for (int i = 0; i < 9; ++i) BackStress.push_back(state[9+i]);
    clone->SetMaterialBackStress(BackStress);
    clone->SetMaterialEquivalentPlasticStrain(state[18]);
  }

  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &strain = static_cast<Tensor_d0s2 &>(_strain);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // Cauchy stress tensor
  std::vector<double> ltangents; // consistent tangent modulus

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = strain[3*i+j];

  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, &ltangents);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];

  delete clone;
}

template<typename Material>
void 
MaterialWrapper<Material>::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                                     double *, double *statenp, double)
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

template<typename Material>
void
MaterialWrapper<Material>::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                                     double *, double *statenp, double)
{
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // first P-K stress tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = enp[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, NULL);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];
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

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];

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

template<>
inline void 
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
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

  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);

  std::vector<double> lstrain; // deformation gradient
  std::vector<double> lstress; // Cauchy stress tensor

  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = enp[3*i+j];

  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, NULL);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];

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

