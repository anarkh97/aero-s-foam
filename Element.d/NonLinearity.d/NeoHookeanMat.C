#include <Element.d/NonLinearity.d/NeoHookeanMat.h>

NeoHookeanMat::NeoHookeanMat(StructProp *p)
{
  double E = p->E;
  double nu = p->nu;
  rho = p->rho; // PJSA
  lambda = E*nu/((1.+nu)*(1.-2.*nu));
  mu = E/(2.*(1.+nu));
}

NeoHookeanMat::NeoHookeanMat(double _rho, double E, double nu)
{
  rho = _rho;
  lambda = E*nu/((1.+nu)*(1.-2.*nu));
  mu = E/(2.*(1.+nu));
}

void
NeoHookeanMat::getStress(Tensor *_stress, Tensor &_strain, double*)
{

}

void 
NeoHookeanMat::getTangentMaterial(Tensor *_tm, Tensor &, double*)
{

}

void 
NeoHookeanMat::getStressAndTangentMaterial(Tensor *_stress, Tensor *_tm, Tensor &_strain, double*)
{

}

void 
NeoHookeanMat::integrate(Tensor *_stress, Tensor *_tm, Tensor &, Tensor &_enp,
                         double *staten, double *statenp, double)
{
  Tensor_d0s4_Ss12s34 *tm = static_cast<Tensor_d0s4_Ss12s34 *>(_tm);
  Tensor_d0s2_Ss12 *stress = static_cast<Tensor_d0s2_Ss12 *>(_stress);
  Tensor_d0s2 &enp = static_cast<Tensor_d0s2 &>(_enp);

  std::vector<double> lstrain; // Green-Lagrange strain
  std::vector<double> lstress; // first P-K stress tensor
  std::vector<double> ltangents; // second (material) elasticity tensor 

  // copy emp into lstrain
  if(lstrain.size()!=9) lstrain.resize(9);
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      lstrain[3*i+j] = lstrain[3*j+i] = enp[i*(5-i)/2+j];

  GetConstitutiveResponse(&lstrain, &lstress, &ltangents);

  // copy lstress into stress
  for (int i = 0; i < 3; i++)
    for (int j = i; j < 3; j++)
      (*stress)[i*(5-i)/2+j] = lstress[3*i+j];

  // copy ltangents into tm
  for(int i=0; i<3; i++)
    for(int j=i; j<3; j++)
      for(int k=0; k<3; k++)
        for(int l=k; l<3; l++)
          (*tm)[i*(5-i)/2+j][k*(5-k)/2+l] = ltangents[i*27+j*9+k*3+l];
}

static double matlib_determinant(double *A)
{
  double det;

  det = A[0]*(A[4]*A[8]-A[5]*A[7])
    -A[1]*(A[3]*A[8]-A[5]*A[6])
    +A[2]*(A[3]*A[7]-A[4]*A[6]);

  return det;
}

static double matlib_inverse(double *A,double *Ainv)
{
  double det,detinv;

  det = matlib_determinant(A);
  if (fabs(det) < 1e-10) return 0.e0;

  detinv = 1./det;
  Ainv[0] = detinv*( A[4]*A[8]-A[5]*A[7]);
  Ainv[1] = detinv*(-A[1]*A[8]+A[2]*A[7]);
  Ainv[2] = detinv*( A[1]*A[5]-A[2]*A[4]);
  Ainv[3] = detinv*(-A[3]*A[8]+A[5]*A[6]);
  Ainv[4] = detinv*( A[0]*A[8]-A[2]*A[6]);
  Ainv[5] = detinv*(-A[0]*A[5]+A[2]*A[3]);
  Ainv[6] = detinv*( A[3]*A[7]-A[4]*A[6]);
  Ainv[7] = detinv*(-A[0]*A[7]+A[1]*A[6]);
  Ainv[8] = detinv*( A[0]*A[4]-A[1]*A[3]);

  return det;
}

bool
NeoHookeanMat::GetConstitutiveResponse(const std::vector<double> * strain,
                                       std::vector<double> * stress,
                                       std::vector<double> * tangents) const
{
  int i,j,k,l,m,n,ij,jj,kl,jk,il,ik,im,jl,kj,kn,mj,nl,ijkl;
  double coef,defVol,detC,p,trace;
  double delta[]  = {1., 0., 0.,
                     0., 1., 0.,
                     0., 0., 1.};
  double C[9],Cinv[9];
  int J;
  int ndf=3;
  int ndm=3;

  /* compute right Cauchy-Green tensor C */
  for(i=0;i<3;i++) for(J=0;J<3;J++) C[i*3+J] = 2*(*strain)[i*3+J] + delta[i*3+J];

  /* compute PK2 stresses and derivatives wrt C*/
  detC = matlib_inverse(C,Cinv);

  if (detC < 1.e-10) {
    std::cerr << "NeoHookean::GetConstitutiveResponse:  close to negative jacobian\n";
    return false;
  }

  defVol = 0.5*log(detC);
  p = lambda*defVol;

  trace = C[0]+C[4]+C[8];

  if (stress) {
    if(stress->size()!=9) stress->resize(9);

    coef = p-mu;
    for (j=0,ij=0,jj=0; j < 3; j++,jj+=4) {
      for (i=0; i < 3; i++,ij++)
        (*stress)[ij] = coef*Cinv[ij];
      (*stress)[jj] += mu;
    }
  }

  if (tangents) {
    if(tangents->size()!=81) tangents->resize(81);

    coef = mu-p;
    for (l=0,kl=0,ijkl=0; l < 3; l++)
      for (k=0,jk=0; k < 3; k++,kl++)
        for (j=0,ij=0,jl=l*3; j < 3; j++,jk++,jl++)
          for (i=0,ik=k*3,il=l*3; i < 3; i++,ij++,ik++,il++,ijkl++)
            (*tangents)[ijkl] = lambda*Cinv[ij]*Cinv[kl]
              +coef*(Cinv[ik]*Cinv[jl]+Cinv[il]*Cinv[jk]);
  }

  return true;
}
