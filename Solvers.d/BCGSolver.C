#include <stdio.h>
#include <Math.d/Vector.h>
#include <Timers.d/GetTime.h>

/*
  All dot products defined as:

                  T
      (x, y) = {x} {y} (operator *)
*/

template<class Scalar, class AnyVector, class AnyOperator>
void
GenBCGSolver<Scalar, AnyVector, AnyOperator>::solve(AnyVector &rhs, AnyVector &sol)
{
 solveTime = -getTime();

 AnyVector r( rhs );

 // ... INITIALIZE SOL TO ZERO
 sol.zero();

 double r0r0 = r.squareNorm();

 //AnyScalar r1r1 = r*r;

 if(r0r0 == 0.0) return;

// ... FIRST ITERATION

 AnyVector p ( r   );
 AnyVector Ap( rhs );

 int niter = 1;

 A->mult(p, Ap);

 while( 1 )
 {

   double r2 = r.squareNorm();
   Scalar rr = r*r;

   // FIRST CHECK CONVERGENCE

  fprintf(stderr," ... Iteration #%d\tTwo norm = %1.7e\n",niter,r2);	

   if( r2 <= sqrt(r0r0*tolerance*tolerance) || niter >= maxiter )  {
      fprintf(stderr,"\n ... Total # Iterations = %d",niter);
      fprintf(stderr,"\n ...     Final Two norm = %1.7e\n",r2);
      solveTime += getTime();
      return;
   }

   A->mult(p, Ap);

   Scalar pAp = p*Ap;

   Scalar alpha = rr/pAp;

   // sol = sol + alpha*p
   sol.linAdd(alpha,p);
  
   // r = r - alpha*Ap 
   r.linAdd(-alpha,Ap);

   Scalar r1r1 = r*r;

   Scalar beta = r1r1/rr;

   // p = r + beta*p
   p.linC(r, beta, p);

   ++niter;

   //For convergence check
   //fprintf(stderr, "\n %d    %1.12e ", niter, r2 );
  
 }

}


