#include <stdio.h>

#include <Solvers.d/ComplexSolver.h>
#include <Math.d/ComplexVector.h>

int
ComplexSolver::numRBM()
{
 return 0;
}

void
ComplexSolver::reSolve(FullMC *) {

 fprintf(stderr,"Selected reSolve does not support FullMC matrix input.\n");

}


void
ComplexSolver::reSolve(int nRHS, DComplex **RHS)
{
 // nRHS = number of rhs
 int i;
 for (i = 0; i < nRHS; ++i)
    reSolve(RHS[i]);
}

void
ComplexSolver::reSolve(int nRHS, DComplex *RHS) {

 fprintf(stderr,"Selected reSolve does not support DComplex* input.\n");

}

void
ComplexSolver::reSolve(int nRHS, ComplexVector *RHS)
{
 int i;
 for (i = 0; i < nRHS; ++i)
   reSolve(RHS[i]);
}

void
ComplexSolver::reSolve(ComplexVector &v)
{
 reSolve(v.data());
}

void
ComplexSolver::reSolve(DComplex*)
{
 fprintf(stderr,"Selected reSolve does not support complex.\n");
}

void
ComplexSolver::solve(DComplex *, DComplex *)
{
 fprintf(stderr,"Selected Solver does not support complex.\n");
}

void
ComplexSolver::solve(ComplexVector &, ComplexVector &)
{
 fprintf(stderr,"Selected Solver does not support a Complex Vector.\n");
}

void
ComplexSolver::factor()
{
 
}

double
ComplexSolver::getSolutionTime()
{
 return 0.0; 
}

long
ComplexSolver::size()
{
 return 0;
}
