// The following specialization is for the Pade-Lanczos
// UH 05/20/08 
// revised by PJSA 10/15/08
//

#include <cassert>
#include <cstdlib>
#include <sys/time.h>

template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::PadeLanczos_BuildSubspace(int nRHS, VecType **sol_prev, VecType **OrthoVec, int numOrthoVec)
{
//
// This routine builds a M-orthonormal basis for the Krylov subspace
// { (K - sigma M)^{-1} rhs, (K - sigma M)^{-1} M (K - sigma M)^{-1} rhs, ... , 
//  [(K - sigma M)^{-1} M]^n (K - sigma M)^{-1} rhs}
//
//--- UH (05/22/08) The random value should be modified for complex problems.
//
// UH (05/21/2008)
//

  if (allOps->C_deriv) {
    std::cerr << "\n !!! The Pade-Lanczos algorithm is not implemented yet";
    std::cerr << " for the damped case !!!\n\n";
    assert(0 > 1);
  }

  srand((unsigned int) time(0));

//  startTimerMemory(times->formRhs, times->memoryRhs);

  //------
  VecType *u = new VecType(probDesc->solVecInfo());
  VecType *w = new VecType(probDesc->solVecInfo());
  //------
  for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
    filePrint(stderr," ... Compute Basis Vector #%3d        ... \n",iRHS);
    //---
    if (iRHS == 0) {
      probDesc->getRHS(*u);
      allOps->sysSolver->solve(*u, *w);
    }
    else {
      //--- Get the previous vector
      *w = *(sol_prev[iRHS-1]);
      //--- Compute M * V(:, iRHS-1) 
      allOps->M->mult(*w, *u);
      if(domain->solInfo().modeFilterFlag) probDesc->project(*u); // PJSA 10-15-08
      //--- Solve the linear system
      allOps->sysSolver->solve(*u, *w);
    }
    //---
    bool hasExpansionForV = false;
    while (hasExpansionForV == false) {
      double mNorm = 0.0;
      std::vector<Scalar> beta(iRHS + numOrthoVec);
      //--- Perform two steps of Gram-Schmidt
      for (int iijj = 0; iijj < 2; ++iijj) {
        allOps->M->mult(*w, *u);
        mNorm = sqrt(ScalarTypes::Real( (*w) * (*u) ));
        (*w) *= 1.0 / mNorm;
        if ((iRHS == 0) && (numOrthoVec == 0)) {
          hasExpansionForV = true;
          break;
        }
        double invMNorm = 1.0 / mNorm;
        //--------
        for (int jj = 0; jj < numOrthoVec; ++jj) {
          beta[jj] = ( (*u) * (*(OrthoVec[jj])) ) * invMNorm;
        }
        for (int jj = 0; jj < iRHS; ++jj) {
          beta[jj + numOrthoVec] = ( (*u) * (*(sol_prev[jj])) ) * invMNorm;
        }
        //--------
        for (int jj = 0; jj < numOrthoVec; ++jj) {
          *u = *(OrthoVec[jj]);
          *u *= -beta[jj];
          (*w) += *u;
        }
        for (int jj = 0; jj < iRHS; ++jj) {
          *u = *(sol_prev[jj]);
          *u *= -beta[jj + numOrthoVec];
          (*w) += *u;
        }
      } // for (int iijj = 0; iijj < 2; ++iijj)
      //---
      allOps->M->mult(*w, *u);
      mNorm = sqrt(ScalarTypes::Real( (*w) * (*u) ));
      if (mNorm < 1e-15) {
        //--- Reset w to random values
        for (int iirr = 0; iirr < w->size(); ++iirr)
          (*w)[iirr] = ((double) rand()) / RAND_MAX;
        continue;
      }
      else {
        hasExpansionForV = true;
        *w *= 1.0 / mNorm;
      }
    } // while (hasExpansionForV == false)

    *(sol_prev[iRHS]) = *w;

  } // for (int iRHS = 0; iRHS < nRHS; ++iRHS)
  //------
  
  delete u;
  delete w;

//  stopTimerMemory(times->formRhs, times->memoryRhs);

//--- Check the M-orthonormality ---  
/*
for (int iVec = 0; iVec < numOrthoVec; ++iVec) {
  VecType *Mv = new VecType(probDesc->solVecInfo());
  allOps->M->mult(*(OrthoVec[iVec]), *Mv);
  for (int jVec = 0; jVec < numOrthoVec; ++jVec) {
    double dot = ScalarTypes::Real( (*OrthoVec[jVec]) * (*Mv) );
    std::cout.precision(3);
    std::cout << dot << " ";
  }
  for (int jRHS = 0; jRHS < nRHS; ++jRHS) {
    double dot = ScalarTypes::Real( (*sol_prev[jRHS]) * (*Mv) );
    std::cout.precision(3);
    std::cout << dot << " ";
  }
  std::cout << std::endl;
  delete Mv;
}
//---
for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
  VecType *Mv = new VecType(probDesc->solVecInfo());
  allOps->M->mult(*(sol_prev[iRHS]), *Mv);
  for (int jVec = 0; jVec < numOrthoVec; ++jVec) {
    double dot = ScalarTypes::Real( (*OrthoVec[jVec]) * (*Mv) );
    std::cout.precision(3);
    std::cout << dot << " ";
  }
  for (int jRHS = 0; jRHS < nRHS; ++jRHS) {
    double dot = ScalarTypes::Real( (*sol_prev[jRHS]) * (*Mv) );
    std::cout.precision(3);
    std::cout << dot << " ";
  }
  std::cout << std::endl;
  delete Mv;
}
*/
//----------------------------------

}


extern "C" {
extern void _FORTRAN(dsysv)(const char &, const int &, const int &,
                 double *, const int &, int *, double *, const int &,
                 double *, int &, int &);
extern void _FORTRAN(zsysv)(const char &, const int &, const int &,
                 complex<double> *, const int &, int *, complex<double> *,
                 const int &, complex<double> *, int &, int &);
void _FORTRAN(dspev)(const char &JOBZ, const char &UPLO,
                     const int &N, double *AP, double *W,
                     double *Z, const int &LDZ, double *WORK, int &INFO);
}

#ifndef _TSYSV__
#define _TSYSV__
inline void Tsysv(const char &a, const int &b, const int &c,
                 double *d, const int &e, int *f, double *g, const int &h,
                 double *i, int &j, int &k) {
 _FORTRAN(dsysv)(a,b,c,d,e,f,g,h,i,j,k);
}
inline void Tsysv(const char &a, const int &b, const int &c,
                 complex<double> *d, const int &e, int *f, complex<double> *g,
                 const int &h,
                 complex<double> *i, int &j, int &k) {
 _FORTRAN(zsysv)(a,b,c,d,e,f,g,h,i,j,k);
}
#endif
/*
void _FORTRAN(dgesv)(const int &, const int &, double *, const int &, int *, double *,
                       const int &, int &);

  void _FORTRAN(zgesv)(const int &, const int &, complex<double> *, const int &, int *, complex<double> *,
                       const int &, int &);


inline void Tgesv(const int &a, const int &b, double *c, const int &d, int *e, double *f,
                  const int &g, int &h)
{
 _FORTRAN(dgesv)(a,b,c,d,e,f,g,h);
}

inline void Tgesv(const int &a, const int &b, complex<double> *c, const int &d, int *e, complex<double> *f,
                  const int &g, int &h)
{
 _FORTRAN(zgesv)(a,b,c,d,e,f,g,h);
}*/




template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::PadeLanczos_Evaluate(int nRHS, VecType **BasisVector, double* &VtKV, double* &Vtb, double w, VecType *sol, double shiftForK)
{
// 
// This routine builds a Pade approximation to H(w) = (K - w^2 M)^{-1} b
// We assume that the basis stored in V is M-orthonormal.
// (The vectors of V are stored in BasisVector).
// The Pade approximation is H_{Pade}(w) = V * (V^T * K * V - w^2 * V^T * M * V)^{-1} * V^T * b
// Note that V^T M V is not computed and set to the identity matrix.
//
// On the first call of this routine, the matrix V^T K V and the vector V^T b are computed.
//
// UH (05/21/2008)
//

  if (VtKV == 0) {
    VtKV = new double[nRHS * nRHS];
    //----
    VecType *u = new VecType(probDesc->solVecInfo());
    for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
      //--
      if (allOps->K == 0) {
        std::cerr << " !!! Stiffness matrix is not created !!! " << std::endl;
        assert(0 > 1);
      }
      allOps->K->mult(*(BasisVector[iRHS]), *u);
      //--
      for (int jRHS = iRHS; jRHS < nRHS; ++jRHS) {
        VtKV[jRHS + iRHS * nRHS] = ScalarTypes::Real( (*u) * (*(BasisVector[jRHS])) );
        VtKV[iRHS + jRHS * nRHS] = VtKV[jRHS + iRHS * nRHS];
      } // for (int jRHS = 0; jRHS < nRHS; ++jRHS)
      //--
    } // for (int iRHS = 0; iRHS < nRHS; ++iRHS)
    //--- 
    //--- allOps->K is assembled with a shift proportional to M.
    //--- The next lines cancel the shift.
    //--- 
    for (int jj = 0; jj < nRHS; ++jj)
      VtKV[jj * (1 + nRHS)] += shiftForK;
    //----
    Vtb = new double[nRHS];
    probDesc->getRHS(*u);
    for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
      Vtb[iRHS] = ScalarTypes::Real( (*u) * (*(BasisVector[iRHS])) );
    }
    //----

////--- Check VtKV ---  
//for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
//  for (int jRHS = 0; jRHS < nRHS; ++jRHS) {
//    std::cout.precision(2);
//    std::cout << VtKV[iRHS + jRHS * nRHS] << " ";
//  }
//  std::cout << std::endl;
//}
////----------------------------------

    // PJSA 10-14-08 solve the reduced eigenvalue problem to get the poles of the Pade approximant
    bool xxxx = false; // change to true if eigenvectors are required
                       // this idea needs a bit of work and testing. esp output freqsweep modes and eigenmodes go to same file
    if(domain->solInfo().pade_poles && nRHS == domain->solInfo().nFreqSweepRHS*domain->solInfo().padeN) {
      char jobz = (xxxx) ? 'V' : 'N';
      double *ap = new double[nRHS*(nRHS+1)/2];
      double *w = new double[nRHS];
      double *z = (xxxx) ? new double[nRHS*nRHS] : 0;
      double *work = new double[3*nRHS];
      int info;
      // fill the ap matrix
      int index = 0;
      for(int j = 0; j < nRHS; ++j)
        for(int i = 0; i <= j; ++i)
          ap[index++] = VtKV[i + j*nRHS];

      _FORTRAN(dspev)(jobz, 'U', nRHS, ap, w, z, nRHS, work, info);

      if(info < 0) filePrint(stderr, "Error in dspev: the %d-th argument had an illegal value.\n", -info);
      else if(info > 0) filePrint(stderr, "Warning in dspev: the algorithm failed to converge; %d\n"
                                          "                  off-diagonal elements of an intermediate tridiagonal"
                                          "                  form did not converge to zero.\n", info);
      filePrint(stderr," --------------------------------------\n");
      filePrint(stderr," Mode\tPoles of Pade Approximant\n");
      filePrint(stderr," --------------------------------------\n");
      int imode = 0;
      for(int i = 0; i < nRHS; ++i) {
        if(w[i] >= domain->solInfo().pade_poles_sigmaL && w[i] <= domain->solInfo().pade_poles_sigmaU) {
          filePrint(stderr, " %d\t%e\n", ++imode, w[i]);
          if(xxxx) {
            u->zero();
            for(int j = 0; j < nRHS; ++j) u->linAdd(z[nRHS*i+j], (*(BasisVector[j])));
            //XXXX domain->postProcessing<Scalar>(*u, (Scalar *)0, *u, 0, w[i]); 
          }
        }
      }
      filePrint(stderr," --------------------------------------\n");
      delete [] ap; delete [] w; delete [] work;
    }

    delete u;
  }

  //--- Shift the projected stiffness matrix
  std::vector<double> copyVtKV(nRHS * nRHS);
  for (int ii = 0; ii < nRHS * nRHS; ++ii)
    copyVtKV[ii] = VtKV[ii];
  for (int jj = 0; jj < nRHS; ++jj)
    copyVtKV[jj * (1 + nRHS)] += - w * w;

  std::vector<double> alpha(nRHS);
  for (int ii = 0; ii < nRHS; ++ii)
    alpha[ii] = Vtb[ii];

/*
      FILE *ff = fopen("pm","a");
      for(int ii=0;ii<nRHS;ii++) for(int jj=0;jj<nRHS;jj++) {
        fprintf(ff,"%d %d %.20e %.20e\n",ii+1,jj+1,
           ScalarTypes::Real(copyVtKV[ii*(nRHS)+jj]),
           ScalarTypes::Imag(copyVtKV[ii*(nRHS)+jj]));
      }
      fclose(ff);
*/


  //--- Solve the reduced linear system
  std::vector<int> ipiv(nRHS);
  int lwork = 3 * nRHS;
  std::vector<double> work(lwork);
  int info = 0;
  _FORTRAN(dsysv)('U', nRHS, 1, &copyVtKV[0], nRHS, &ipiv[0], &alpha[0], nRHS, &work[0], lwork, info);

//for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
//  std::cout << alpha[iRHS] << std::endl;
//}

  //--- Compute the Pade approximant vector
  for (int jj = 0; jj < sol->size(); ++jj) {
    Scalar value = 0.0;
    for (int iRHS = 0; iRHS < nRHS; ++iRHS) {
      value += (*(BasisVector[iRHS]))[jj] * alpha[iRHS];
    } // for (int iRHS = 0; iRHS < nRHS; ++iRHS)
    (*sol)[jj] = value;
  } // for (int jj = 0; jj < sol->size(); ++jj)
    
}




template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::galProjection(bool gpReorthoFlag,
                  int nRHS, VecType *sol, VecType **u, VecType **v,
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  double w, double deltaw)
{
// fprintf(stderr,"w deltaw size: %f %f %d\n",w,deltaw,sol->size());
/*     for(int i=0;i<nRHS;i++) 
       for(int j=0;j<sol->size();j++) 
         fprintf(stderr,"sss %d %d %e %e\n",i,j,
ScalarTypes::Real((*u[i])[j]),
ScalarTypes::Imag((*u[i])[j]));
*/

 VecType *f = new VecType(probDesc->solVecInfo());

 if (gpReorthoFlag) {
   VecType *a = new VecType(probDesc->solVecInfo());

   for(int i=0;i<nRHS;i++) *u[i] = *v[i];

   if(domain->solInfo().isCoupled) {
     for(int i=0;i<nRHS;i++) {
        scaleInvDisp(*u[i]);
     }
   }
 // Orthogonalize
   int ngs = 5;
   for(int m=0;m<ngs;m++) { 
     for(int i=0;i<nRHS;i++) {
       for(int j=0;j<i;j++) {
         Scalar dotp = *u[i] *  *u[j];
         (*u[i]).linAdd(-dotp,*u[j]);
//         if (m==ngs-1) fprintf(stderr,"dot %d %d %e %e\n",
//                 i,j,ScalarTypes::Real(dotp),ScalarTypes::Imag(dotp));
       }
       Scalar nrm = *u[i] * *u[i];
//       if (m==0) fprintf(stderr,"nrm %d %.20e\n",
//                 i,sqrt(ScalarTypes::Real(nrm)));
       *u[i] *= 1.0/sqrt(ScalarTypes::Real(nrm));
     }
   }


   VecType *b = new VecType(probDesc->solVecInfo());
   VecType *c = new VecType(probDesc->solVecInfo());
   if (VhMV==0) VhMV = new Scalar[(nRHS)*(nRHS)];
   if (VhKV==0) VhKV = new Scalar[(nRHS)*(nRHS)];
   if (VhCV==0) VhCV = new Scalar[(nRHS)*(nRHS)];
   for(int i=0;i<nRHS;i++) {
     
// for(int k=0;k<sol->size();k++)
//   fprintf(stderr,"u[%d][%d]= %.18e %.18e\n",i,k+1,ScalarTypes::Real((*(u[i]))[k]),ScalarTypes::Imag((*(u[i]))[k]));

     allOps->K->mult(*(u[i]), *a);
// for(int k=0;k<sol->size();k++)
//   fprintf(stderr,"a[%d]= %.18e %.18e\n",k+1,ScalarTypes::Real((*a)[k]),ScalarTypes::Imag((*a)[k]));

     allOps->M->mult(*(u[i]), *b);
// for(int k=0;k<sol->size();k++)
//   fprintf(stderr,"b[%d]= %.18e %.18e\n",k+1,ScalarTypes::Real((*b)[k]),ScalarTypes::Imag((*b)[k]));


     for(int k=0;k<sol->size();k++) (*c)[k] = 0;
     if (allOps->C_deriv) {
       if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*(u[i]), *c);
     } else (*c).zero(); 
     //else if (allOps->C) allOps->C->mult(*(u[i]), *c);  

     for(int j=0;j<nRHS;j++) {
       VhMV[i*(nRHS)+j] = *b * *u[j];
       VhCV[i*(nRHS)+j] = *c * *u[j];
       VhKV[i*(nRHS)+j] = *a * *u[j] + (w-deltaw)*(w-deltaw) * VhMV[i*(nRHS)+j];
       ScalarTypes::addComplex(VhKV[i*(nRHS)+j], -(w-deltaw)*VhCV[i*(nRHS)+j] );
     }
   }
   delete b;
   delete c;
/*
   fprintf(stderr,"K\n");
   for(int i=0;i<sol->size();i++) {
     for(int j=0;j<sol->size();j++) (*f)[j] = 0;
     (*f)[i] = 1;
     allOps->K->mult(*f, *a);
     for(int j=0;j<sol->size();j++) {
      if (sqrt(ScalarTypes::Real((*a)[j])*ScalarTypes::Real((*a)[j])+
               ScalarTypes::Imag((*a)[j])*ScalarTypes::Imag((*a)[j])) > 0.0 )
      fprintf(stderr,"%d %d: %.18e %.18e \n",i,j,ScalarTypes::Real((*a)[j]),ScalarTypes::Imag((*a)[j]));
     }
   }
   fprintf(stderr,"M\n");
   for(int i=0;i<sol->size();i++) {
     for(int j=0;j<sol->size();j++) (*f)[j] = 0;
     (*f)[i] = 1;
     allOps->M->mult(*f, *a);
     for(int j=0;j<sol->size();j++) {
      if (sqrt(ScalarTypes::Real((*a)[j])*ScalarTypes::Real((*a)[j])+
               ScalarTypes::Imag((*a)[j])*ScalarTypes::Imag((*a)[j])) > 0.0 )
      fprintf(stderr,"%d %d: %.18e %.18e \n",i,j,ScalarTypes::Real((*a)[j]),ScalarTypes::Imag((*a)[j]));
     }
   }
   if (allOps->C_deriv) {
    if (allOps->C_deriv[0])  {
     fprintf(stderr,"allOps->C_deriv[0]\n");
     for(int i=0;i<sol->size();i++) {
       for(int j=0;j<sol->size();j++) (*f)[j] = 0;
       (*f)[i] = 1;
       allOps->C_deriv[0]->mult(*f, *a);
       for(int j=0;j<sol->size();j++) {
        if (sqrt(ScalarTypes::Real((*a)[j])*ScalarTypes::Real((*a)[j])+
                 ScalarTypes::Imag((*a)[j])*ScalarTypes::Imag((*a)[j])) > 0.0 )
        fprintf(stderr,"%d %d: %e %e \n",i,j,ScalarTypes::Real((*a)[j]),ScalarTypes::Imag((*a)[j]));
       }
     }
    }
   }
   if (allOps->C) {
     fprintf(stderr,"allOps->C\n");
     for(int i=0;i<sol->size();i++) {
       for(int j=0;j<sol->size();j++) (*f)[j] = 0;
       (*f)[i] = 1;
       allOps->C->mult(*f, *a);
       for(int j=0;j<sol->size();j++) {
        if (sqrt(ScalarTypes::Real((*a)[j])*ScalarTypes::Real((*a)[j])+
                 ScalarTypes::Imag((*a)[j])*ScalarTypes::Imag((*a)[j])) > 0.0 )
        fprintf(stderr,"%d %d: %e %e \n",i,j,ScalarTypes::Real((*a)[j]),ScalarTypes::Imag((*a)[j]));
       }
     }
   }
*/
   delete a;
 }

 // Project 
 probDesc->getRHS(*f, w,deltaw);
// for(int k=0;k<sol->size();k++)
//   fprintf(stderr,"f[%d]= %d %e %e\n",k+1,sizeof((*f)[0]),ScalarTypes::Real((*f)[k]),ScalarTypes::Imag((*f)[k]));


 Scalar *z = new Scalar[nRHS];
 for(int i=0;i<nRHS;i++) {
   z[i] = *f * *u[i];
 }

 delete f;
 
// fprintf(stderr,"allOps->C_deriv %p\n",allOps->C_deriv);

 Scalar *zz = new Scalar[(nRHS)*(nRHS)];
 for(int i=0;i<nRHS;i++)
   for(int j=0;j<nRHS;j++) {
     zz[i*(nRHS)+j] = VhKV[i*(nRHS)+j] - w*w * VhMV[i*(nRHS)+j];
     ScalarTypes::addComplex(zz[i*(nRHS)+j], 
                      w * VhCV[i*(nRHS)+j]);
   }

/*
      FILE *ff = fopen("pm","a");
      for(int ii=0;ii<nRHS;ii++) for(int jj=0;jj<nRHS;jj++) {
        fprintf(ff,"%d %d %.20e %.20e\n",ii+1,jj+1,
           ScalarTypes::Real(zz[ii*(nRHS)+jj]),
           ScalarTypes::Imag(zz[ii*(nRHS)+jj]));
      }
      fclose(ff);
*/

//--- Solve the reduced linear system
 std::vector<int> ipiv(nRHS);
 int lwork = 3 * nRHS;
 std::vector<Scalar> work(lwork);
 int info = 0;
 Tgesv(nRHS, 1, &zz[0], nRHS, &ipiv[0], &z[0], nRHS, info);

 fprintf(stderr,"Coefs:");
 for(int i=0;i<nRHS;i++)
   fprintf(stderr," %e %e",ScalarTypes::Real(z[i]),
                          ScalarTypes::Imag(z[i]));
 fprintf(stderr,"\n");

//--- Compute the approximant vector
 (*sol).zero();
 for(int i=0;i<nRHS;i++) 
   (*sol).linAdd(z[i],*u[i]);

 if(domain->solInfo().isCoupled) scaleDisp(*sol);

 delete[] z;
 delete[] zz;
}


template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::krylovGalProjection(int nRHS, int nOrtho,
                  VecType *sol, VecType **u, VecType **v,
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  double w, double deltaw)
{
// fprintf(stderr,"w deltaw size: %f %f %d %d %d\n",w,deltaw,nRHS,nOrtho,sol->size());

 VecType *f = new VecType(probDesc->solVecInfo());

 if (nRHS>0) {
   VecType *a = new VecType(probDesc->solVecInfo());

   for(int i=0;i<nOrtho;i++) *u[i] = *v[i];
   if(domain->solInfo().isCoupled) {
//     for(int i=0;i<nOrtho;i++)
//fprintf(stderr,"gogo %d %e\n",i,ScalarTypes::Real((*u[i])[0]));
     for(int i=0;i<nOrtho;i++)
        scaleInvDisp(*u[i]);
     for(int i=0;i<nOrtho;i++) {
       *a = *u[i];
       int ngs = 2;
       for(int m=0;m<ngs;m++) { 
          for(int j=0;j<i;j++) {
            Scalar dotp = *a *  *u[j];
            (*a).linAdd(-dotp,*u[j]);
          }
          Scalar nrm = *a * *a;
          *a *= 1.0/sqrt(ScalarTypes::Real(nrm));
       }
       *u[i] = *a;
     }
   }

   for(int i=0;i<nRHS;i++) {
     if (i == 0) {
       probDesc->getRHS(*f);
     }
     else {
       allOps->M->mult(*u[nOrtho+i-1], *f);
     }
     filePrint(stderr,"\n ... Solving RHS #%3d               ...\n",i);
     allOps->sysSolver->solve(*f, *a);
     int ngs = 2;
     for(int m=0;m<ngs;m++) { 
        for(int j=0;j<nOrtho+i;j++) {
          Scalar dotp = *a *  *u[j];
          (*a).linAdd(-dotp,*u[j]);
        }
        Scalar nrm = *a * *a;
        *a *= 1.0/sqrt(ScalarTypes::Real(nrm));
     }
     *u[nOrtho+i] = *a;
   }

   VecType *b = new VecType(probDesc->solVecInfo());
   VecType *c = new VecType(probDesc->solVecInfo());

   if (VhMV!=0) delete[] VhMV; 
   if (VhKV!=0) delete[] VhKV; 
   if (VhCV!=0) delete[] VhCV; 
   VhMV = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];
   VhKV = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];
   VhCV = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];

   for(int i=0;i<nRHS+nOrtho;i++) {
     allOps->K->mult(*(u[i]), *a);
     allOps->M->mult(*(u[i]), *b);
     for(int k=0;k<sol->size();k++) (*c)[k] = 0;
     if (allOps->C_deriv) {
       if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*(u[i]), *c);
     } else (*c).zero(); 
     //else if (allOps->C) allOps->C->mult(*(u[i]), *c);  

     for(int j=0;j<nOrtho+nRHS;j++) {
       VhMV[i*(nOrtho+nRHS)+j] = *b * *u[j];
       VhCV[i*(nOrtho+nRHS)+j] = *c * *u[j];
       VhKV[i*(nOrtho+nRHS)+j] = *a * *u[j] + (w-deltaw)*(w-deltaw) * VhMV[i*(nOrtho+nRHS)+j];
       ScalarTypes::addComplex(VhKV[i*(nOrtho+nRHS)+j], -(w-deltaw)*VhCV[i*(nOrtho+nRHS)+j] );
     }
   }

   delete b;
   delete c;
   delete a;

   for(int i=0;i<nOrtho+nRHS;i++) *v[i] = *u[i];
   if(domain->solInfo().isCoupled) {
     for(int i=0;i<nOrtho+nRHS;i++) 
        scaleDisp(*v[i]);
   }
/*
for(int i=0;i<nOrtho+nRHS;i++) 
fprintf(stderr,"gaga %d %e\n",i,ScalarTypes::Real((*u[i])[0]));
for(int i=0;i<nOrtho+nRHS;i++) 
fprintf(stderr,"gugu %d %e\n",i,ScalarTypes::Real((*v[i])[0]));
*/
 }

 // Project 
 probDesc->getRHS(*f, w,deltaw);

 Scalar *z = new Scalar[nOrtho+nRHS];
 for(int i=0;i<nOrtho+nRHS;i++) {
   z[i] = *f * *u[i];
 }
 delete f;

 Scalar *zz = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];
 for(int i=0;i<nOrtho+nRHS;i++)
   for(int j=0;j<nOrtho+nRHS;j++) {
     zz[i*(nOrtho+nRHS)+j] = VhKV[i*(nOrtho+nRHS)+j] - w*w * VhMV[i*(nOrtho+nRHS)+j];
     ScalarTypes::addComplex(zz[i*(nOrtho+nRHS)+j], 
                      w * VhCV[i*(nOrtho+nRHS)+j]);
   }

//--- Solve the reduced linear system
 std::vector<int> ipiv(nOrtho+nRHS);
 int lwork = 3 * (nOrtho+nRHS);
 std::vector<Scalar> work(lwork);
 int info = 0;
 Tgesv(nOrtho+nRHS, 1, &zz[0], nOrtho+nRHS, &ipiv[0], &z[0], nOrtho+nRHS, info);

/* fprintf(stderr,"Coefs:");
 for(int i=0;i<nOrtho+nRHS;i++)
   fprintf(stderr," %e %e",ScalarTypes::Real(z[i]),
                          ScalarTypes::Imag(z[i]));
 fprintf(stderr,"\n");
*/

//--- Compute the approximant vector
 (*sol).zero();
 for(int i=0;i<nOrtho+nRHS;i++) 
   (*sol).linAdd(z[i],*u[i]);

 if(domain->solInfo().isCoupled) scaleDisp(*sol);

 delete[] z;
 delete[] zz;
}


template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
void
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::qrGalProjection(int nRHS, int nOrtho,
                  VecType *sol, VecType **u, VecType **v,
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  double w, double deltaw)
{
 fprintf(stderr,"w deltaw size: %f %f %d %d %d\n",w,deltaw,nRHS,nOrtho,sol->size());

 VecType *f = new VecType(probDesc->solVecInfo());
 if (nRHS>0) {
  
   VecType *a = new VecType(probDesc->solVecInfo());
   Scalar *y = new Scalar[nRHS];
   Scalar *R = new Scalar[nRHS*nRHS];
   Scalar *S = new Scalar[nRHS*nRHS];
   Scalar *dotp = new Scalar[nRHS];
  
   for(int i=0;i<nRHS*nRHS;i++) R[i] = 0.0;
  
  // u stores q (columns of Q); v stores z (columns of Z)
   for(int p=0;p<nRHS;p++) {
     if (p==0)
       probDesc->getRHS(*f);
     else 
       probDesc->getFreqSweepRHS(f,0, p);
     allOps->sysSolver->solve(*f, *a);
     for(int i=0;i<p;i++) {
       y[i] = 2.0*double(p)*w*R[(p-1)*nRHS+i];
       if (p>1) y[i] += double(p)*double(p-1)*R[(p-2)*nRHS+i];
      }

     for(int i=0;i<p;i++) {
       (*a).linAdd(y[i],*v[i]);
        fprintf(stderr,"y[%d]= %.20e %.20e\n",i, ScalarTypes::Real(y[i]),
           ScalarTypes::Imag(y[i]));
     }
     for(int i=0;i<p;i++) {
       R[p*nRHS+i] = *a * *u[i];
     }
     for(int i=0;i<p;i++) {
        fprintf(stderr,"R[%d]= %.20e %.20e\n",i, ScalarTypes::Real(R[p*nRHS+i]),
           ScalarTypes::Imag(R[p*nRHS+i]));
       (*a).linAdd(-R[p*nRHS+i],*u[i]);
     }

/*
     for(int i=0;i<p;i++) {
       R[p*nRHS+i] = 0;
       R[p*nRHS+i] += *a * *u[i];
       for(int j=0;j<p;j++) {
         R[p*nRHS+i] += S[j*nRHS+i] * y[j];
       }
     }
     for(int i=0;i<p;i++) {
        fprintf(stderr,"y[%d]= %.20e %.20e\n",i, ScalarTypes::Real(y[i]),
           ScalarTypes::Imag(y[i]));
       (*a).linAdd(y[i],*v[i]);
        fprintf(stderr,"R[%d]= %.20e %.20e\n",i, ScalarTypes::Real(R[p*nRHS+i]),
           ScalarTypes::Imag(R[p*nRHS+i]));
       (*a).linAdd(-R[p*nRHS+i],*u[i]);
     }
*/

  // Extra reortho
     int ngs = 2;
     for(int m=0;m<ngs;m++) { 
       for(int i=0;i<p;i++) {
         dotp[i] = *a *  *u[i];
         fprintf(stderr,"%d %d %.20e %.20e\n",m,i,
           ScalarTypes::Real(dotp[i]),
           ScalarTypes::Imag(dotp[i]));
       }
       for(int i=0;i<p;i++) {
         (*a).linAdd(-dotp[i],*u[i]);
         R[p*nRHS+i] += dotp[i];
       }
     }
  
     Scalar nrm = *a * *a;
        fprintf(stderr,"nrm %.20e %.20e\n",
           ScalarTypes::Real(nrm),
           ScalarTypes::Imag(nrm));
     R[p*nRHS+p] = sqrt(ScalarTypes::Real(nrm));
     *a *= 1.0/sqrt(ScalarTypes::Real(nrm));
     *u[p] = *a;
  
     allOps->M->mult(*u[p], *f);
     allOps->sysSolver->solve(*f, *a);
     *v[p] = *a;
     for(int i=0;i<p;i++) {
       S[p*nRHS+i] = *v[p] * *u[i];
       S[i*nRHS+p] = *v[i] * *u[p];
     }
   }

   for(int p=0;p<nRHS;p++) for(int i=0;i<nRHS;i++) {
    Scalar d = *v[p] * *v[i];
        fprintf(stderr,"Z^HZ %d %d  %.20e %.20e\n",p,i,
           ScalarTypes::Real(d),
           ScalarTypes::Imag(d));
   }

   VecType *b = new VecType(probDesc->solVecInfo());
   VecType *c = new VecType(probDesc->solVecInfo());

   if (VhMV!=0) delete[] VhMV; 
   if (VhKV!=0) delete[] VhKV; 
   if (VhCV!=0) delete[] VhCV; 
   VhMV = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];
   VhKV = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];
   VhCV = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];

   for(int i=0;i<nRHS+nOrtho;i++) {
     allOps->K->mult(*(u[i]), *a);
     allOps->M->mult(*(u[i]), *b);
     for(int k=0;k<sol->size();k++) (*c)[k] = 0;
     if (allOps->C_deriv) {
       if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*(u[i]), *c);
     } else (*c).zero(); 
     //else if (allOps->C) allOps->C->mult(*(u[i]), *c);  

     for(int j=0;j<nOrtho+nRHS;j++) {
       VhMV[i*(nOrtho+nRHS)+j] = *b * *u[j];
       VhCV[i*(nOrtho+nRHS)+j] = *c * *u[j];
       VhKV[i*(nOrtho+nRHS)+j] = *a * *u[j] + (w-deltaw)*(w-deltaw) * VhMV[i*(nOrtho+nRHS)+j];
       ScalarTypes::addComplex(VhKV[i*(nOrtho+nRHS)+j], -(w-deltaw)*VhCV[i*(nOrtho+nRHS)+j] );
     }
   }

   delete b;
   delete c;
   delete a;
 }

 // Project 
 probDesc->getRHS(*f, w,deltaw);

 Scalar *z = new Scalar[nOrtho+nRHS];
 for(int i=0;i<nOrtho+nRHS;i++) {
   z[i] = *f * *u[i];
 }
 delete f;

 Scalar *zz = new Scalar[(nOrtho+nRHS)*(nOrtho+nRHS)];
 for(int i=0;i<nOrtho+nRHS;i++)
   for(int j=0;j<nOrtho+nRHS;j++) {
     zz[i*(nOrtho+nRHS)+j] = VhKV[i*(nOrtho+nRHS)+j] - w*w * VhMV[i*(nOrtho+nRHS)+j];
     ScalarTypes::addComplex(zz[i*(nOrtho+nRHS)+j], 
                      w * VhCV[i*(nOrtho+nRHS)+j]);
   }

//--- Solve the reduced linear system
 std::vector<int> ipiv(nOrtho+nRHS);
 int lwork = 3 * (nOrtho+nRHS);
 std::vector<Scalar> work(lwork);
 int info = 0;
 Tgesv(nOrtho+nRHS, 1, &zz[0], nOrtho+nRHS, &ipiv[0], &z[0], nOrtho+nRHS, info);

 fprintf(stderr,"Coefs:");
 for(int i=0;i<nOrtho+nRHS;i++)
   fprintf(stderr," %e %e",ScalarTypes::Real(z[i]),
                          ScalarTypes::Imag(z[i]));
 fprintf(stderr,"\n");

//--- Compute the approximant vector
 (*sol).zero();
 for(int i=0;i<nOrtho+nRHS;i++) 
   (*sol).linAdd(z[i],*u[i]);

 if(domain->solInfo().isCoupled) scaleDisp(*sol);

 delete[] z;
 delete[] zz;
}
