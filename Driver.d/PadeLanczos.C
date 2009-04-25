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
void _FORTRAN(dspev)(const char &JOBZ, const char &UPLO,
                     const int &N, double *AP, double *W,
                     double *Z, const int &LDZ, double *WORK, int &INFO);
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

