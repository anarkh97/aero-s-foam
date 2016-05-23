/*
 * NeoHookean.cpp
 * DG++
 *
 * Created by Adrian Lew on 10/24/06.
 *  
 * Copyright (c) 2006 Adrian Lew
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */ 


#include <Material.d/Material.h>
#include <cmath>
#include <iostream>


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




static void matlib_mults(double *A, double *B,double *C)
{
  C[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2] ; 
  C[1] = A[3] * B[0] + A[4] * B[1] + A[5] * B[2] ; 
  C[2] = A[6] * B[0] + A[7] * B[1] + A[8] * B[2] ; 
  C[3] = A[0] * B[3] + A[1] * B[4] + A[2] * B[5] ; 
  C[4] = A[3] * B[3] + A[4] * B[4] + A[5] * B[5] ; 
  C[5] = A[6] * B[3] + A[7] * B[4] + A[8] * B[5] ; 
  C[6] = A[0] * B[6] + A[1] * B[7] + A[2] * B[8] ; 
  C[7] = A[3] * B[6] + A[4] * B[7] + A[5] * B[8] ; 
  C[8] = A[6] * B[6] + A[7] * B[7] + A[8] * B[8] ; 
}





bool NeoHookean::GetConstitutiveResponse(const std::vector<double> * strain,
					 std::vector<double> * stress,
					 std::vector<double> * tangents) const
{
  int i,j,k,l,m,n,ij,jj,kl,jk,il,ik,im,jl,kj,kn,mj,nl,ijkl,indx;
  double coef,defVol,detC,p,trace;
  double F[]  = {1., 0., 0.,
		 0., 1., 0.,
		 0., 0., 1.};
  double C[9],Cinv[9],S[9],M[81];
  //double Finv[9],detF;
  int J;
  int ndf=3;
  int ndm=3;

  /*Fill in the deformation gradient*/
  for(i=0;i<ndf;i++) for(J=0;J<ndm;J++) F[i*3+J]=(*strain)[i*ndm+J];

  /* compute right Cauchy-Green tensor C */
  matlib_mults(F,F,C);

  /* compute PK2 stresses and derivatives wrt C*/
  detC = matlib_inverse(C,Cinv);
  //detF = matlib_inverse(F,Finv);

  if (detC < 1.e-10) {
    std::cerr << "NeoHookean::GetConstitutiveResponse:  close to negative jacobian\n";
    return false;
  }
  
  defVol = 0.5*log(detC);
  p = Lambda*defVol;

  trace = C[0]+C[4]+C[8];

  coef = p-Mu;

  for (j=0,ij=0,jj=0; j < 3; j++,jj+=4) {
    for (i=0; i < 3; i++,ij++)
      S[ij] = coef*Cinv[ij];
    S[jj] += Mu;
  }


  if (tangents) {
    coef = Mu-p;
    for (l=0,kl=0,ijkl=0; l < 3; l++)
      for (k=0,jk=0; k < 3; k++,kl++)
	for (j=0,ij=0,jl=l*3; j < 3; j++,jk++,jl++)
	  for (i=0,ik=k*3,il=l*3; i < 3; i++,ij++,ik++,il++,ijkl++)
	    M[ijkl] = Lambda*Cinv[ij]*Cinv[kl]
	      +coef*(Cinv[ik]*Cinv[jl]+Cinv[il]*Cinv[jk]);
  }
  
  if(stress->size()!=9) stress->resize(9);

  /* PK2 -> PK1 */
  for (j=0,ij=0; j < ndm; j++)
    for(i=0; i < ndm; i++,ij++) {
      (*stress)[ij] = 0.e0;
      for (k=0,ik=i,kj=j*3; k < 3; k++,ik+=3,kj++)
	(*stress)[ij] += F[ik]*S[kj];
    }


  if (tangents) {
    if(tangents->size()!=81) tangents->resize(81);

    /* apply partial push-forward and add geometrical term */
    for (l=0,ijkl=0; l < ndm; l++)
      for (k=0; k < ndf; k++)
	for (j=0,jl=l*3; j < ndm; j++,jl++)
	  for (i=0; i < ndf; i++,ijkl++) {

	    (*tangents)[ijkl] = 0.e0;
	    /* push-forward */
	    for (n=0,kn=k,nl=l*3; n < 3; n++,kn+=3,nl++) {
	      indx = nl*9;
	      for (m=0,im=i,mj=j*3; m < 3; m++,im+=3,mj++)
		(*tangents)[ijkl] += F[im]*M[mj+indx]*F[kn];
	    }
	    
	    /* geometrical term */
	    if (i == k)
	      (*tangents)[ijkl] += S[jl];
	  }
  }

  return true;
}
