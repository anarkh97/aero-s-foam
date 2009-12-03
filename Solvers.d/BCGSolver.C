#include <stdio.h>
#include <Math.d/Vector.h>
#include <Timers.d/GetTime.h>

/*
  Bi-CG is for real nonsymmetric linear systems
*/

template<class Scalar, class AnyVector, class AnyOperator, class AnyPreconditioner>
void
GenBCGSolver<Scalar, AnyVector, AnyOperator, AnyPreconditioner>::solve(AnyVector &rhs, AnyVector &sol)
{
  double resid;
  Scalar rho_1, rho_2, alpha, beta;
  AnyVector z(rhs), ztilde(rhs), p(rhs), ptilde(rhs), q(rhs), qtilde(rhs);

  double normb = rhs.norm();
  sol.zero();
  AnyVector r(rhs); // r = b - A * x  (ASSUMING x^0 = 0)
  AnyVector rtilde(r);

  if (normb == 0.0)
    normb = 1.0;
  
  if ((resid = r.norm() / normb) <= tolerance) {
    return;
  }

  for (int i = 1; i <= maxiter; i++) {
    fprintf(stderr," ... Iteration #%d\tTwo norm = %1.7e\tRelative residual = %1.7e\n",i,r.sqNorm(),resid);
    if(P) P->apply(r, z); else z=r; // z = M^{1}*r
    if(P) P->apply(rtilde, ztilde); else ztilde = rtilde; // ztilde = M^{-T}*rtilde (ASSUMING here that preconditioner is symmetric)
    rho_1 = z*rtilde;
    if (rho_1 == 0.0) { 
      return;
    }
    if (i == 1) {
      p = z;
      ptilde = ztilde;
    } else {
      beta = rho_1 / rho_2;
      p.linC(z, beta, p); // p = z + beta*p
      ptilde.linC(ztilde, beta, ptilde); // ptilde = ztilde + beta*ptilde
    }
    A->mult(p, q); // q = A * p
    A->transposeMult(ptilde,qtilde); //qtilde = A^T*ptilde;
    alpha = rho_1 / (ptilde*q);
    sol.linAdd(alpha, p); // sol += alpha * p;
    r.linAdd(-alpha, q); // r -= alpha*q
    rtilde.linAdd(-alpha, qtilde); // rtilde -= alpha*qtilde

    rho_2 = rho_1;
    if ((resid = r.norm() / normb) < tolerance) {
      return;
    }
  }

  return;
}


