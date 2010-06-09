#include <Element.d/NonLinearity.d/NLMaterial.h>
#include <Element.d/NonLinearity.d/SimpleMat.h>
#include <Material.d/Material.h>

//#define PJSA_CONVERT_TO_PK2
//#define PJSA_CONVERT_TO_M

SimpleMat::SimpleMat(int type, double rho, double E, double nu)
{
  double lambda = E*nu/((1.+nu)*(1.-2.*nu));
  double mu = E/(2.*(1.+nu));
  switch(type) {
   case 1:
     mat = new IsotropicLinearElastic(lambda, mu, rho);
     break;
   case 2:
     mat = new NeoHookean(lambda, mu, rho);
     break;
  }
}

void
SimpleMat::getStress(Tensor *_stress, Tensor &_strain, double*)
{

}

void 
SimpleMat::getTangentMaterial(Tensor *_tm, Tensor &, double*)
{

}

void 
SimpleMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*)
{

}

void 
SimpleMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                     double *staten, double *statenp, double)
{
  Tensor_d0s4 *tm = static_cast<Tensor_d0s4 *>(_tm);
  Tensor_d0s2 *stress = static_cast<Tensor_d0s2 *>(_stress);
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);

  std::vector<double> lstrain; // gradU + I
  std::vector<double> lstress; // first P-K stress tensor
  std::vector<double> ltangents; // first elasticity tensor

  // copy emp into lstrain
  lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      lstrain[3*i+j] = enp[3*i+j];

  mat->GetConstitutiveResponse(&lstrain, &lstress, &ltangents);

#ifndef PJSA_CONVERT_TO_PK2
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      (*stress)[3*i+j] = lstress[3*i+j];
#else
  Tensor_d0s2 P;
  // copy lstress into P (this is the first P-K stress tensor)
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      P[3*i+j] = lstress[3*i+j];

  // compute the second P-K stress tensor
  Tensor_d0s2 Finv;
  enp.getInverse(Finv);
  (*stress) = Finv|P;
#endif

#ifndef PJSA_CONVERT_TO_M
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          (*tm)[i*27+j*9+k*3+l] = ltangents[i*27+j*9+k*3+l];

#else
  // XXXX also need to convert first elasticity tensor (ltangents) to second elasticity tensor
  std::cerr << "in SimpleMat::integrate, conversion from first elasticity tensor to material elasticity tensor not implemented\n";
#endif
}

double 
SimpleMat::getDensity()
{ 
  return mat->GetDensityInReference();
} 
