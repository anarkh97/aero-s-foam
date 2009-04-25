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
GenCRSolver<Scalar, AnyVector, AnyOperator>::solve(AnyVector &rhs, AnyVector &sol)
{
 solveTime = -getTime();

 AnyVector r(rhs);

 // ... INITIALIZE SOL TO ZERO
 sol.zero(); 

 double r0r0 = r.squareNorm();

 if(r0r0 == 0.0) return;

// ... FIRST ITERATION

 AnyVector p ( r   );
 AnyVector Ap( rhs );
 AnyVector Ar( rhs );

 int niter = 1;

 A->mult(p, Ap);

// ofstream writefile("rel_residual",ios::out);

 while( 1 )
 {

   // FIRST CHECK CONVERGENCE
   double r2 = r.squareNorm();

   fprintf(stderr," ... Iteration #%d\tTwo norm = %1.7e\t  Rel. residual  %1.7e \n",niter,r2,sqrt(r2/r0r0));
//   writefile << sqrt(r2/r0r0) << endl;
   if( r2 <= r0r0*tolerance*tolerance || niter >= maxiter )  {
      fprintf(stderr,"\n ... Total # Iterations = %d",niter);
      fprintf(stderr,"\n ...     Final Two norm = %1.7e\n",r2);
      solveTime += getTime();
//      writefile.close();
      return;
   }

   A->mult(r, Ar);

   Scalar ApAp = Ap*Ap;

   Scalar Ar_r = Ar*r;

   Scalar alpha = Ar_r/ApAp;

   // sol = sol + alpha*p
   sol.linAdd(alpha,p);
  
   // r = r - alpha*Ap 
   r.linAdd(-alpha,Ap);

   A->mult(r, Ar);

   /* 
      Book's choice. Works with the dot 
      product defined as:

                  T
      (x, y) = {x} {y}
   */
   //Scalar beta = (Ar^r)/(Ar_r);
   Scalar beta = (Ar*r)/(Ar_r);

   /* 
      Michel's choice. Slow convergence and needs the 
      following definition of the dot product:

                  *
      (x, y) = {x} {y}  (note: * is transpose conjugate)
   */
 
   //Scalar beta = -(Ap^Ar)/(ApAp);

   p.linC(r, beta, p);

   // Ap = Ar + beta*Ap
   Ap.linC(Ar, beta, Ap);

   ++niter;

  
 }

}


