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
MaterialWrapper<Material>::getStress(Tensor *_stress, Tensor &_strain, double*, double temp)
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
MaterialWrapper<NeoHookean>::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
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
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > F(&lstrain[0]), P(&lstress[0]), s(&(*stress)[0]);
  s = F.inverse()*P; // symmetric 2nd Piola-Kirchhoff stress tensor, S = F^{-1}*P
#else
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];
#endif
}

template<>
inline void
MaterialWrapper<MooneyRivlin>::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
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
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > F(&lstrain[0]), P(&lstress[0]), s(&(*stress)[0]);
  s = F.inverse()*P; // symmetric 2nd Piola-Kirchhoff stress tensor, S = F^{-1}*P
#else
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];
#endif
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
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

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getStress(Tensor *_stress, Tensor &_strain, double* state, double temp)
{
  std::cerr << "ERROR: MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getStress is not implemented\n";
}

template<typename Material>
void 
MaterialWrapper<Material>::getTangentMaterial(Tensor *_tm, Tensor &_strain, double*, double temp)
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
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getTangentMaterial(Tensor *_tm, Tensor &_strain, double* state, double temp)
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

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getTangentMaterial(Tensor *_tm, Tensor &_strain, double* state, double temp)
{
  std::cerr << "ERROR: MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getTangentMaterial is not implemented\n";
}

template<typename Material>
void 
MaterialWrapper<Material>::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*, double temp)
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
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double* state, double temp)
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

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double* state, double temp)
{
  std::cerr << "ERROR: MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getStressAndTangentMaterial is not implemented\n";
}

template<typename Material>
void 
MaterialWrapper<Material>::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                                     double *, double *statenp, double, double)
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
                                     double *, double *statenp, double, double)
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
                                                                    double *staten, double *statenp, double, double dt)
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

  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, &ltangents, true, dt);

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
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                                                                               double *staten, double *statenp, double, double)
{
  std::cerr << "ERROR: MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::integrate is not implemented\n";
}

template<>
inline void 
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                                                                    double *staten, double *statenp, double, double dt)
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

  clone->ComputeElastoPlasticConstitutiveResponse(lstrain, &lstress, NULL, true, dt);

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

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::integrate(Tensor *_stress, Tensor &, Tensor &_enp,
                                                                               double *staten, double *statenp, double, double)
{
  std::cerr << "ERROR: MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::integrate is not implemented\n";
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
  return 0.0;
}

template<>
inline double
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::getDensity()
{
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
inline double
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::getEquivPlasticStrain(double *statenp)
{
  return statenp[18];
}

template<typename Material>
double
MaterialWrapper<Material>::getStrainEnergyDensity(Tensor &_enp, double *, double)
{
  std::cerr << "WARNING: MaterialWrapper<Material>::getStrainEnergyDensity is not implemented\n";
  return 0.0;
}

template<typename Material>
void
MaterialWrapper<Material>::setSDProps(MFTTData *ysst)
{
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::setSDProps(MFTTData *ysst)
{
  double SigmaY = mat->GetYieldStressFromTensionTest();
  if(SigmaY < 0 && ysst && ysst->getID() == -int(SigmaY)) {
    for(int i=0; i<ysst->getNumPoints(); ++i) {
      mat->SetExperimentalCurveData(ysst->getT(i), ysst->getV(i));
    }
  }
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::setSDProps(MFTTData *ysst)
{
  double SigmaY = mat->GetYieldStressFromTensionTest();
  if(SigmaY < 0 && ysst && ysst->getID() == -int(SigmaY)) {
    for(int i=0; i<ysst->getNumPoints(); ++i) {
      mat->SetExperimentalCurveData(ysst->getT(i), ysst->getV(i));
    }
  }
}

template<typename Material>
void
MaterialWrapper<Material>::setSRDProps(MFTTData *yssrt)
{
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::setSRDProps(MFTTData *yssrt)
{
  if(yssrtid > 0 && yssrt && yssrt->getID() == yssrtid) {
    for(int i=0; i<yssrt->getNumPoints(); ++i) {
      mat->SetExperimentalCurveData2(yssrt->getT(i), yssrt->getV(i));
    }
  }
}

template<>
inline void
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::setSRDProps(MFTTData *yssrt)
{
  if(yssrtid > 0 && yssrt && yssrt->getID() == yssrtid) {
    for(int i=0; i<yssrt->getNumPoints(); ++i) {
      mat->SetExperimentalCurveData2(yssrt->getT(i), yssrt->getV(i));
    }
  }
}

#ifdef USE_EIGEN3
template<>
inline double
MaterialWrapper<NeoHookean>::getStrainEnergyDensity(Tensor &_enp, double *, double temp)
{
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > F(&enp[0]);

  using std::log;
  double lnJ = log(F.determinant());
  return mu/2*((F*F.transpose()).trace() - 3) - mu*lnJ + lambda/2*lnJ*lnJ;
}

template<>
inline double
MaterialWrapper<MooneyRivlin>::getStrainEnergyDensity(Tensor &_enp, double *, double temp)
{
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);
  Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> > F(&enp[0]);

  using std::log;
  Eigen::Matrix<double,3,3> C = F.transpose()*F;
  double I1 = C.trace();
  double I2 = 0.5*(I1*I1-(C*C).trace());
  double J = F.determinant();
  double d = 2*(mu1+2*mu2);
  return mu1*(I1-3) + mu2*(I2-3) + kappa*(J-1)*(J-1) - d*log(J);
}
#endif

template<>
inline void
MaterialWrapper<IsotropicLinearElastic>::print(std::ostream &out) const
{
  double rho = mat->GetDensityInReference();
  double E = mu*(3*lambda+2*mu)/(lambda+mu);
  double nu = lambda/(2*(lambda+mu));
  out << "IsotropicLinearElastic " << rho << " " << E << " " << nu;
}

template<>
inline void 
MaterialWrapper<NeoHookean>::print(std::ostream &out) const
{
  double rho = mat->GetDensityInReference();
  double E = mu*(3*lambda+2*mu)/(lambda+mu);
  double nu = lambda/(2*(lambda+mu));
  out << "NeoHookean " << rho << " " << E << " " << nu;
}

template<>
inline void 
MaterialWrapper<IsotropicLinearElasticJ2PlasticMaterial>::print(std::ostream &out) const
{
  double rho = 0.0;
  double E = mu*(3*lambda+2*mu)/(lambda+mu);
  double nu = lambda/(2*(lambda+mu));
  double sigmaY = mat->GetYieldStressFromTensionTest();
  double K = mat->GetIsotropicHardeningModulus();
  double H = mat->GetKinematicHardeningModulus();
  double Tol = mat->GetTolerance();
  double epsF = mat->GetEquivalentPlasticStrainAtFailure();
  out << "IsotropicLinearElasticJ2PlasticMaterial " << rho << " " << " " << E << " " << nu << " " << sigmaY << " " << K << " " << H
      << " " << Tol << " " << epsF << " " << yssrtid;
}

template<>
inline void 
MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial>::print(std::ostream &out) const
{
  double rho = 0.0;
  double E = mu*(3*lambda+2*mu)/(lambda+mu);
  double nu = lambda/(2*(lambda+mu));
  double sigmaY = mat->GetYieldStressFromTensionTest();
  double K = mat->GetIsotropicHardeningModulus();
  double H = mat->GetKinematicHardeningModulus();
  double Tol = mat->GetTolerance();
  double epsF = mat->GetEquivalentPlasticStrainAtFailure();
  out << "IsotropicLinearElasticJ2PlasticPlaneStressMaterial " << rho << " " << E << " " << nu << " " << sigmaY << " " << K << " " << H
      << " " << Tol << " " << epsF << " " << yssrtid;
}

template<>
inline void 
MaterialWrapper<MooneyRivlin>::print(std::ostream &out) const
{
  double rho = mat->GetDensityInReference();
  out << "MooneyRivlin " << rho << " " << mu1 << " " << mu2 << " " << kappa;
}

template<typename Material>
void
MaterialWrapper<Material>::getMaterialConstants(std::vector<double> &c)
{
  std::cerr << "material law does not implement getMaterialConstants function\n";
}

template<>
inline void
MaterialWrapper<NeoHookean>::getMaterialConstants(std::vector<double> &c)
{
  c.resize(2);
  c[0] = lambda;
  c[1] = mu;
}

template<>
inline void
MaterialWrapper<MooneyRivlin>::getMaterialConstants(std::vector<double> &c)
{
  c.resize(3);
  c[0] = mu1;
  c[1] = mu2;
  c[2] = kappa;
}
