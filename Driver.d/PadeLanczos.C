// The following specialization is for the Pade-Lanczos
// UH 05/20/08 
// revised by PJSA 10/15/08
//

#include <cassert>
#include <cstdlib>
#include <sys/time.h>
#include <Timers.d/GetTime.h>

       #include <sys/types.h>
       #include <sys/stat.h>
       #include <fcntl.h>



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
double
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::galProjection(bool gpReorthoFlag,
                  int nRHS, VecType *sol, VecType **u, VecType **v,
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  double w, double deltaw)
{
 filePrint(stderr,"w deltaw size: %f %f %d\n",w,deltaw,sol->size());


 VecType *f = new VecType(probDesc->solVecInfo());

 if (gpReorthoFlag) {

double time = 0.0;
time -= getTime();
   VecType *a = new VecType(probDesc->solVecInfo());

   for(int i=0;i<nRHS;i++) {
     *u[i] = *v[i];
   }

   if(domain->solInfo().isCoupled) {
     for(int i=0;i<nRHS;i++) {
        scaleInvDisp(*u[i]);
     }
   }

   for(int i=0;i<nRHS;i++) {
filePrint(stderr,"Forcing continuity %d.\n",i);
     forceContinuity(*u[i]);
//fprintf(stderr,"Dumping deriv %d.\n",i);
//       postProcessor->staticOutput(*u[i], *u[i],false);
   }


/*
   double *buf = new double[sol->size()];
   int h = open("vectors",O_CREAT|O_APPEND|O_WRONLY,S_IRUSR|S_IWUSR);
   
     for(int i=0;i<nRHS;i++) {
       for(int j=0;j<sol->size();j++) buf[j] = ScalarTypes::Real((*v[i])[j]);
       write(h,buf,sizeof(double)*sol->size());
    }
   close(h);
*/


// for(int i=0;i<nRHS;i++) ((GenDistrVector<Scalar>*)(u[i]))->printAll();
 // Orthogonalize
//   int ngs = 5;
   int ngs = 2;
   for(int m=0;m<ngs;m++) { 
     for(int i=0;i<nRHS;i++) {
       Scalar nrm0 = *u[i] * *u[i];
       for(int j=0;j<i;j++) {
         Scalar dotp = *u[i] *  *u[j];
         (*u[i]).linAdd(-dotp,*u[j]);
//         if (m==ngs-1)
//             fprintf(stderr,"dot %d %d %.16e %.16e\n",
//                 i,j,ScalarTypes::Real(dotp),ScalarTypes::Imag(dotp));
       }
       Scalar nrm = *u[i] * *u[i];
//       if (m==0)
//          fprintf(stderr,"nrm %d %d %.16e %.16e\n",m,
//                 i,sqrt(ScalarTypes::Real(nrm0)),sqrt(ScalarTypes::Real(nrm)));
       *u[i] *= 1.0/sqrt(ScalarTypes::Real(nrm));
       nrm = *u[i] * *u[i];
//          fprintf(stderr,"nrmnrm %d %d %.16e \n",m,
//                 i,sqrt(ScalarTypes::Real(nrm)));
//       if (m==ngs-1) postProcessor->staticOutput(*u[i], *rhs,false);
     }
   }

/*
   double *buf = new double[sol->size()];
   int h = open("ovectors",O_CREAT|O_APPEND|O_WRONLY,S_IRUSR|S_IWUSR);
   
     for(int i=0;i<nRHS;i++) {
       for(int j=0;j<sol->size();j++) buf[j] = ScalarTypes::Real((*u[i])[j]);
       write(h,buf,sizeof(double)*sol->size());
    }
   close(h);
   delete[] buf;
*/



/*
     for(int i=0;i<nRHS;i++) {
       for(int j=0;j<nRHS;j++) {
         Scalar dotp = *u[i] *  *u[j];
         fprintf(stderr,"check dot %d %d %.16e %.16e\n",
                 i,j,ScalarTypes::Real(dotp),ScalarTypes::Imag(dotp));
       }
     }*/

   VecType *b = new VecType(probDesc->solVecInfo());
   VecType *c = new VecType(probDesc->solVecInfo());
   if (VhMV!=0) delete[] VhMV; 
   if (VhKV!=0) delete[] VhKV; 
   if (VhCV!=0) delete[] VhCV; 
   VhMV = new Scalar[(nRHS)*(nRHS)];
   VhKV = new Scalar[(nRHS)*(nRHS)];
   VhCV = new Scalar[(nRHS)*(nRHS)];
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
       if (allOps->C_deriv[0]) { allOps->C_deriv[0]->mult(*(u[i]), *c);
          fprintf(stderr,"kukuku\n");
       }
     } else (*c).zero(); 
     //else if (allOps->C) allOps->C->mult(*(u[i]), *c);  

     for(int j=0;j<nRHS;j++) {
       VhMV[i*(nRHS)+j] = *b * *u[j];
       VhCV[i*(nRHS)+j] = *c * *u[j];
       VhKV[i*(nRHS)+j] = *a * *u[j];
// + (w-deltaw)*(w-deltaw) * VhMV[i*(nRHS)+j];
//       ScalarTypes::addComplex(VhKV[i*(nRHS)+j], -(w-deltaw)*VhCV[i*(nRHS)+j] );
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
time += getTime();
filePrint(stderr,"Ortho + proj setup time: %e\n",time);
/*
   double *bufm = new double[2*nRHS*nRHS];
   h = open("pmatrices",O_CREAT|O_APPEND|O_WRONLY,S_IRUSR|S_IWUSR);
   
   for(int ii=0;ii<nRHS;ii++) for(int jj=0;jj<nRHS;jj++) {
      bufm[ii*(nRHS)+jj] = ScalarTypes::Real(VhKV[ii*(nRHS)+jj]);
      bufm[nRHS*nRHS+ii*(nRHS)+jj] = ScalarTypes::Real(VhMV[ii*(nRHS)+jj]);
   }
   write(h,bufm,sizeof(double)*2*nRHS*nRHS);
   close(h);
   delete[] bufm;
*/
 }

double time = 0.0;
time -= getTime();

 // Project 
 probDesc->getRHS(*f, w,deltaw);
// for(int k=0;k<sol->size();k++)
//   fprintf(stderr,"f[%d]= %d %e %e\n",k+1,sizeof((*f)[0]),ScalarTypes::Real((*f)[k]),ScalarTypes::Imag((*f)[k]));


 Scalar *z = new Scalar[nRHS];
 for(int i=0;i<nRHS;i++) {
   z[i] = *f * *u[i];
 }

 
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

/*
 fprintf(stderr,"Coefs:");
 for(int i=0;i<nRHS;i++)
   fprintf(stderr," %e %e",ScalarTypes::Real(z[i]),
                          ScalarTypes::Imag(z[i]));
 fprintf(stderr,"\n");
*/

//--- Compute the approximant vector
 (*sol).zero();
 for(int i=0;i<nRHS;i++) 
   (*sol).linAdd(z[i],*u[i]);

time += getTime();
filePrint(stderr,"Projection time: %e\n",time);

time = 0.0;
time -= getTime();
 VecType *a = new VecType(probDesc->solVecInfo());
 VecType *b = new VecType(probDesc->solVecInfo());
 VecType *c = new VecType(probDesc->solVecInfo());
 allOps->K->mult(*sol, *a);
 allOps->M->mult(*sol, *b);
 if (allOps->C_deriv) {
    if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*sol, *c);
 } else (*c).zero(); 
 (*a).linAdd(-w*w,*b);
 for(int k=0;k<sol->size();k++) 
     ScalarTypes::addComplex((*a)[k], 
                      w * (*c)[k]);
// (*a).print();
 (*a).linAdd(-1.0,*f);
 forceAssemble(*a);
 forceAssemble(*f);
 Scalar nrma = *a * *a;
 Scalar nrmf = *f * *f;
 filePrint(stderr,"residual: %e\n", sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf)));
// (*sol).print();
time += getTime();
filePrint(stderr,"Residual compute time: %e\n",time);


 delete a;
 delete b;
 delete c;

 delete f;

 if(domain->solInfo().isCoupled) scaleDisp(*sol);

 delete[] z;
 delete[] zz;

 return sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf));
}


template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
double
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::krylovGalProjection(int nRHS, int nOrtho,
                  VecType *sol, VecType **u, VecType **v,
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  double w, double deltaw)
{
 filePrint(stderr,"w deltaw size: %f %f %d %d %d\n",w,deltaw,nRHS,nOrtho,sol->size());

 VecType *f = new VecType(probDesc->solVecInfo());

 if (nRHS>0) {
double time = 0.0;
double rhstime = 0.0;
time -= getTime();
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
time += getTime();
rhstime -= getTime();
     allOps->sysSolver->solve(*f, *a);
rhstime += getTime();
time -= getTime();
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
       VhKV[i*(nOrtho+nRHS)+j] = *a * *u[j];
// + (w-deltaw)*(w-deltaw) * VhMV[i*(nOrtho+nRHS)+j];
//      ScalarTypes::addComplex(VhKV[i*(nOrtho+nRHS)+j], -(w-deltaw)*VhCV[i*(nOrtho+nRHS)+j] );
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
time += getTime();
filePrint(stderr,"Krylov setup time: %e\n",time);
filePrint(stderr,"RHS solve time: %e\n",rhstime);
 }

double time = 0.0;
time -= getTime();
 // Project 
 probDesc->getRHS(*f, w,deltaw);

 Scalar *z = new Scalar[nOrtho+nRHS];
 for(int i=0;i<nOrtho+nRHS;i++) {
   z[i] = *f * *u[i];
 }

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
time += getTime();
filePrint(stderr,"Projection time: %e\n",time);

time = 0.0;
time -= getTime();
 VecType *a = new VecType(probDesc->solVecInfo());
 VecType *b = new VecType(probDesc->solVecInfo());
 VecType *c = new VecType(probDesc->solVecInfo());
 allOps->K->mult(*sol, *a);
 allOps->M->mult(*sol, *b);
 if (allOps->C_deriv) {
    if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*sol, *c);
 } else (*c).zero(); 
 (*a).linAdd(-w*w,*b);
 for(int k=0;k<sol->size();k++) 
     ScalarTypes::addComplex((*a)[k], 
                      w * (*c)[k]);
// (*a).print();
 (*a).linAdd(-1.0,*f);
 forceAssemble(*a);
 forceAssemble(*f);
 Scalar nrma = *a * *a;
 Scalar nrmf = *f * *f;
 filePrint(stderr,"residual: %e\n", sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf)));
// (*sol).print();
time += getTime();
filePrint(stderr,"Residual compute time: %e\n",time);

 delete a;
 delete b;
 delete c;

 delete f;


 if(domain->solInfo().isCoupled) scaleDisp(*sol);

 delete[] z;
 delete[] zz;

 return sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf));
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
  ::adaptGP(int dgpFlag, int minRHS, int maxRHS, int deltaRHS, int &nOrtho,
                  VecType *sol, VecType **u, VecType **v,
                  VecType **aa, VecType **bb, VecType **cc, 
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  int ncheck, double *wcheck,
                  double alpha, double tol)
{
 filePrint(stderr,"minRHS,maxRHS,deltaRHS,nOrtho...: %d %d %d  %d %d\n",minRHS,maxRHS,deltaRHS,nOrtho,sol->size());

double time = 0.0;
double rhstime = 0.0;
double projmattime = 0.0;
double projmattime2 = 0.0;
double orthotime = 0.0;
double chrestime = 0.0;

time -= getTime();

 for(int i=0;i<nOrtho;i++) {
   *u[i] = *v[i];
   if(domain->solInfo().isCoupled) scaleDisp(*u[i],1.0/alpha);
 }

 VecType *a = new VecType(probDesc->solVecInfo());
 VecType *b = new VecType(probDesc->solVecInfo());
 VecType *c = new VecType(probDesc->solVecInfo());
 VecType *f = new VecType(probDesc->solVecInfo());

 if (VhKV!=0) delete[] VhKV; 
 if (VhMV!=0) delete[] VhMV; 
 if (VhCV!=0) delete[] VhCV; 
 Scalar *tmpVhKV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
 Scalar *tmpVhMV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
 Scalar *tmpVhCV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];

projmattime -= getTime();
 for(int i=0;i<nOrtho;i++) {
//   allOps->K->mult(*(u[i]), *a);
//   allOps->M->mult(*(u[i]), *b);
   allOps->K->mult(*(u[i]), *(aa[i]));
   allOps->M->mult(*(u[i]), *(bb[i]));
   (*c).zero();
   if (allOps->C_deriv) 
//     if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*(u[i]), *c);
     if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*(u[i]), *(cc[i]));
   for(int j=0;j<nOrtho;j++) {
//     tmpVhKV[i*(nOrtho+maxRHS)+j] = *a * *(u[j]);
//     tmpVhMV[i*(nOrtho+maxRHS)+j] = *b * *(u[j]);
//     tmpVhCV[i*(nOrtho+maxRHS)+j] = *c * *(u[j]);
     tmpVhKV[i*(nOrtho+maxRHS)+j] = *(aa[i]) * *(u[j]);
     tmpVhMV[i*(nOrtho+maxRHS)+j] = *(bb[i]) * *(u[j]);
     tmpVhCV[i*(nOrtho+maxRHS)+j] = *(cc[i]) * *(u[j]);
   }
 }
projmattime += getTime();

 Scalar *z = new Scalar[nOrtho+maxRHS];
 double *rescheck = new double[ncheck];
 double *oldrescheck = new double[ncheck];
 for(int i=0;i<ncheck;i++) oldrescheck[i] = 0.0;
 bool done = false;
 int lRHS = 0;
 int uRHS = minRHS;
 double maxres = 0.0;
 double oldmaxres = 0.0;
 while (!done) {
 
// Add new vectors 
orthotime -= getTime();
   int offset = nOrtho+maxRHS;
   for(int i=lRHS;i<uRHS;i++) {
     if (dgpFlag) {
       if (i == 0) {
         u[offset]->zero();
         probDesc->getRHS(*f);
       } else {
         f->zero();
         probDesc->getFreqSweepRHS(f, u+offset, i);
       }
     } else {
       if (i == 0) {
         probDesc->getRHS(*f);
       }
       else {
         allOps->M->mult(*u[nOrtho+i-1], *f);
       }
     }
     filePrint(stderr,"\n ... Solving RHS #%3d               ...\n",i);
rhstime -= getTime();
     allOps->sysSolver->solve(*f, *a);
rhstime += getTime();
     if (dgpFlag) *u[offset+i+1] = *a;
     if(domain->solInfo().isCoupled) scaleDisp(*a,alpha);
     int ngs = 2;
     for(int m=0;m<ngs;m++) { 
        for(int j=0;j<nOrtho+i;j++) {
          Scalar dotp = *a *  *v[j];
          (*a).linAdd(-dotp,*v[j]);
        }
        Scalar nrm = *a * *a;
        *a *= 1.0/sqrt(ScalarTypes::Real(nrm));
     }
     *v[nOrtho+i] = *a;
     *u[nOrtho+i] = *a;
      if(domain->solInfo().isCoupled) scaleDisp(*u[nOrtho+i],1.0/alpha);
   }
orthotime += getTime();

// Update matrices
projmattime2 -= getTime();
   for(int i=nOrtho+lRHS;i<nOrtho+uRHS;i++) {
//     allOps->K->mult(*(u[i]), *a);
//     allOps->M->mult(*(u[i]), *b);
     allOps->K->mult(*(u[i]), *(aa[i]));
     allOps->M->mult(*(u[i]), *(bb[i]));
     (*c).zero();
     if (allOps->C_deriv) 
//       if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*(u[i]), *c);
       if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*(u[i]), *(cc[i]));
     for(int j=0;j<nOrtho+uRHS;j++) {
//       tmpVhKV[i*(nOrtho+maxRHS)+j] = *a * *(u[j]);
//       tmpVhMV[i*(nOrtho+maxRHS)+j] = *b * *(u[j]);
//       tmpVhCV[i*(nOrtho+maxRHS)+j] = *c * *(u[j]);
       tmpVhKV[i*(nOrtho+maxRHS)+j] = *(aa[i]) * *(u[j]);
       tmpVhMV[i*(nOrtho+maxRHS)+j] = *(bb[i]) * *(u[j]);
       tmpVhCV[i*(nOrtho+maxRHS)+j] = *(cc[i]) * *(u[j]);
     }
   }
   for(int i=0;i<nOrtho+uRHS;i++) {
//     allOps->K->mult(*(u[i]), *a);
//     allOps->M->mult(*(u[i]), *b);
//     (*c).zero();
//     if (allOps->C_deriv) 
//       if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*(u[i]), *c);
     for(int j=nOrtho+lRHS;j<nOrtho+uRHS;j++) {
//       tmpVhKV[i*(nOrtho+maxRHS)+j] = *a * *(u[j]);
//       tmpVhMV[i*(nOrtho+maxRHS)+j] = *b * *(u[j]);
//       tmpVhCV[i*(nOrtho+maxRHS)+j] = *c * *(u[j]);
       tmpVhKV[i*(nOrtho+maxRHS)+j] = *(aa[i]) * *(u[j]);
       tmpVhMV[i*(nOrtho+maxRHS)+j] = *(bb[i]) * *(u[j]);
       tmpVhCV[i*(nOrtho+maxRHS)+j] = *(cc[i]) * *(u[j]);
     }
   }
projmattime2 += getTime();

//Check residual at check points

chrestime -= getTime();
   Scalar *zz = new Scalar[(nOrtho+uRHS)*(nOrtho+uRHS)];
   for(int icheck=0;icheck<ncheck;icheck++) {
     probDesc->getRHS(*f, wcheck[icheck],0.0);
     for(int i=0;i<nOrtho+uRHS;i++) {
       z[i] = *f * *u[i];
     }
     for(int i=0;i<nOrtho+uRHS;i++) {
       for(int j=0;j<nOrtho+uRHS;j++) {
         zz[i*(nOrtho+uRHS)+j] = tmpVhKV[i*(nOrtho+maxRHS)+j] -
                  wcheck[icheck]*wcheck[icheck] * tmpVhMV[i*(nOrtho+maxRHS)+j];
         ScalarTypes::addComplex(zz[i*(nOrtho+uRHS)+j], 
                    wcheck[icheck] * tmpVhCV[i*(nOrtho+maxRHS)+j]);
       }
     }

     std::vector<int> ipiv(nOrtho+uRHS);
     int lwork = 3 * (nOrtho+uRHS);
     std::vector<Scalar> work(lwork);
     int info = 0;
     Tgesv(nOrtho+uRHS, 1, &zz[0], nOrtho+uRHS, &ipiv[0], &z[0],
           nOrtho+uRHS, info);

     (*a).zero();
     for(int i=0;i<nOrtho+uRHS;i++) (*a).linAdd(z[i],*aa[i]);
     (*b).zero();
     for(int i=0;i<nOrtho+uRHS;i++) (*b).linAdd(z[i],*bb[i]);
     (*c).zero();
     for(int i=0;i<nOrtho+uRHS;i++) (*c).linAdd(z[i],*cc[i]);
/*
     (*sol).zero();
     for(int i=0;i<nOrtho+uRHS;i++) 
       (*sol).linAdd(z[i],*u[i]);
     allOps->K->mult(*sol, *a);
     allOps->M->mult(*sol, *b);
     if (allOps->C_deriv) {
        if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*sol, *c);
     } else (*c).zero(); 
*/
     (*a).linAdd(-wcheck[icheck]*wcheck[icheck],*b);
     for(int k=0;k<sol->size();k++) 
       ScalarTypes::addComplex((*a)[k], wcheck[icheck] * (*c)[k]);
     (*a).linAdd(-1.0,*f);
     forceAssemble(*a);
     forceAssemble(*f);
     Scalar nrma = *a * *a;
     Scalar nrmf = *f * *f;
     rescheck[icheck] = 
      sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf));
    filePrint(stderr,"rescheck[%d] = %e at %e\n",icheck,rescheck[icheck],wcheck[icheck]);
   }
   delete[] zz;
chrestime += getTime();

   maxres = 0.0;
   for(int i=0;i<ncheck;i++) if (rescheck[i]>maxres) maxres = rescheck[i];
   filePrint(stderr,"maxres = %e, oldmaxres = %e  tol = %e \n",
             maxres,oldmaxres,tol);
   if (maxres<tol) done = true;
   int nsmaller=0;
//   for(int i=0;i<ncheck;i++) if (oldrescheck[i]!=0.0) if (rescheck[i]>tol) if (rescheck[i]<pow(0.8,deltaRHS)*oldrescheck[i]) nsmaller++; 
// For coupled
   for(int i=0;i<ncheck;i++) if (oldrescheck[i]!=0.0) if (rescheck[i]>tol) if (rescheck[i]<pow(0.9,deltaRHS)*oldrescheck[i]) nsmaller++; 
   if (oldmaxres!=0.0) if (nsmaller==0) done = true;
   for(int i=0;i<ncheck;i++) oldrescheck[i] = rescheck[i];
   oldmaxres = maxres;
   lRHS = uRHS;
   uRHS += deltaRHS;
   if (uRHS>maxRHS) done = true; 
 }

 VhKV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
 VhMV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];
 VhCV = new Scalar[(nOrtho+maxRHS)*(nOrtho+maxRHS)];


 for(int i=0;i<nOrtho+lRHS;i++) {
    for(int j=0;j<nOrtho+lRHS;j++) {
      VhKV[i*(nOrtho+lRHS)+j] = tmpVhKV[i*(nOrtho+maxRHS)+j];
      VhMV[i*(nOrtho+lRHS)+j] = tmpVhMV[i*(nOrtho+maxRHS)+j];
      VhCV[i*(nOrtho+lRHS)+j] = tmpVhCV[i*(nOrtho+maxRHS)+j];
    }
 }

 nOrtho += lRHS;
 filePrint(stderr,"nOrtho = %d\n",nOrtho);
 
 delete[] tmpVhKV;
 delete[] tmpVhMV;
 delete[] tmpVhCV;
 delete[] rescheck;
 delete[] oldrescheck;
 delete[] z;

 delete f;
 delete c;
 delete b;
 delete a;

time += getTime();
filePrint(stderr,"Total setup time: %e\n",time);
filePrint(stderr,"Matrix setup time: %e\n",projmattime);
filePrint(stderr,"Matrix setup time2: %e\n",projmattime2);
filePrint(stderr,"Ortho+ time: %e\n",orthotime-rhstime);
filePrint(stderr,"RHS solve time: %e\n",rhstime);
filePrint(stderr,"Check res time: %e\n",chrestime);
}


template < class Scalar,
           class OpSolver,
           class VecType,
           class PostProcessor,
           class ProblemDescriptor,
           class ComplexVecType>
double
StaticSolver< Scalar, OpSolver, VecType,
              PostProcessor, ProblemDescriptor, ComplexVecType >
  ::adaptGPSolRes(int dgpFlag, int nOrtho,
                  VecType *sol, VecType **u, VecType **v,
                  VecType **aa, VecType **bb, VecType **cc,
                  Scalar *&VhKV, Scalar *&VhMV, Scalar *&VhCV,
                  double w, double deltaw)
 //                  ,double alpha)
{
 filePrint(stderr,"w deltaw size: %f %f  %d %d\n",w,deltaw,nOrtho,sol->size());

double time = 0.0;
time -= getTime();

// This should be already done when this function is called
/*
 for(int i=0;i<nOrtho;i++) {
   *u[i] = *v[i];
   if(domain->solInfo().isCoupled) scaleDisp(*u[i],1.0/alpha);
 }
*/

 // Project 
 VecType *f = new VecType(probDesc->solVecInfo());
 probDesc->getRHS(*f, w,deltaw);

 Scalar *z = new Scalar[nOrtho];
 for(int i=0;i<nOrtho;i++) {
   z[i] = *f * *u[i];
 }

 Scalar *zz = new Scalar[(nOrtho)*(nOrtho)];
 for(int i=0;i<nOrtho;i++) {
   for(int j=0;j<nOrtho;j++) {
     zz[i*(nOrtho)+j] = VhKV[i*(nOrtho)+j] - w*w * VhMV[i*(nOrtho)+j];
     ScalarTypes::addComplex(zz[i*(nOrtho)+j], 
                      w * VhCV[i*(nOrtho)+j]);
   }
 }

//--- Solve the reduced linear system
 std::vector<int> ipiv(nOrtho);
 int lwork = 3 * (nOrtho);
 std::vector<Scalar> work(lwork);
 int info = 0;
 Tgesv(nOrtho, 1, &zz[0], nOrtho, &ipiv[0], &z[0], nOrtho, info);

/* fprintf(stderr,"Coefs:");
 for(int i=0;i<nOrtho+nRHS;i++)
   fprintf(stderr," %e %e",ScalarTypes::Real(z[i]),
                          ScalarTypes::Imag(z[i]));
 fprintf(stderr,"\n");
*/

//--- Compute the approximant vector
 (*sol).zero();
 for(int i=0;i<nOrtho;i++) 
   (*sol).linAdd(z[i],*u[i]);
time += getTime();
filePrint(stderr,"Projection time: %e\n",time);

time = 0.0;
time -= getTime();
 VecType *a = new VecType(probDesc->solVecInfo());
 VecType *b = new VecType(probDesc->solVecInfo());
 VecType *c = new VecType(probDesc->solVecInfo());

     (*a).zero();
     for(int i=0;i<nOrtho;i++) (*a).linAdd(z[i],*aa[i]);
     (*b).zero();
     for(int i=0;i<nOrtho;i++) (*b).linAdd(z[i],*bb[i]);
     (*c).zero();
     for(int i=0;i<nOrtho;i++) (*c).linAdd(z[i],*cc[i]);
/*
 allOps->K->mult(*sol, *a);
 allOps->M->mult(*sol, *b);
 if (allOps->C_deriv) {
    if (allOps->C_deriv[0]) allOps->C_deriv[0]->mult(*sol, *c);
 } else (*c).zero(); 
*/
 (*a).linAdd(-w*w,*b);
 for(int k=0;k<sol->size();k++) 
     ScalarTypes::addComplex((*a)[k], w * (*c)[k]);
 (*a).linAdd(-1.0,*f);
 forceAssemble(*a);
 forceAssemble(*f);
 Scalar nrma = *a * *a;
 Scalar nrmf = *f * *f;
 filePrint(stderr,"residual: %e\n", sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf)));
time += getTime();
filePrint(stderr,"Residual compute time: %e\n",time);

 delete a;
 delete b;
 delete c;
 delete f;
 delete[] z;
 delete[] zz;

 if(domain->solInfo().isCoupled) scaleDisp(*sol);

 return sqrt(ScalarTypes::Real(nrma))/sqrt(ScalarTypes::Real(nrmf));
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
 filePrint(stderr,"w deltaw size: %f %f %d %d %d\n",w,deltaw,nRHS,nOrtho,sol->size());

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

