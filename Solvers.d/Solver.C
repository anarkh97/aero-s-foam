#include <stdio.h>
#include <Utils.d/dbg_alloca.h>

#include <Math.d/Vector.h>

template<class Scalar> 
void
GenSolver<Scalar>::reSolve(int nRHS, Scalar **RHS)
{
 int i;
 for (i = 0; i < nRHS; ++i)
    reSolve(RHS[i]);
}

template<class Scalar> 
void
GenSolver<Scalar>::reSolve(int nRHS, Scalar *rhs)
{

  int numUncon=neqs();
  Scalar **rhsP = (Scalar**) dbg_alloca(sizeof(Scalar*)*nRHS);
  int i;
  for(i=0; i<nRHS; ++i)
    rhsP[i] = rhs+i*numUncon;
  reSolve(nRHS,rhsP);
}


template<class Scalar> 
void
GenSolver<Scalar>::reSolve(int nRHS, GenVector<Scalar> *RHS)
{ 
 int i;
 for (i = 0; i < nRHS; ++i)
   reSolve(RHS[i]);
}

/*
template<class Scalar> 
void
GenSolver<Scalar>::reSolve(GenVector<Scalar> &v)
{
 reSolve(v.data());
}
*/

template<>
void
GenSolver<DComplex>::reSolve(ComplexVector &v);

template<>
void
GenSolver<double>::reSolve(ComplexVector &v);

template<>
void
GenSolver<double>::reSolve(Vector &v);

template<>
void
GenSolver<DComplex>::reSolve(Vector &v);


template<class Scalar> 
void
GenSolver<Scalar>::reSolve(Scalar*)
{
 fprintf(stderr,"Selected reSolve(Scalar*) not supported\n");
}

template<class Scalar>
void
GenSolver<Scalar>::reSolve(GenFullM<Scalar> *)
{
 fprintf(stderr,"Selected reSolve not supported\n");
}

template<class Scalar>
void
GenSolver<Scalar>::forward(GenVector<Scalar>&)
{
 fprintf(stderr,"forward not implemented for selected solver\n");
}

template<class Scalar>
void
GenSolver<Scalar>::backward(GenVector<Scalar>&)
{
 fprintf(stderr,"backward not implemented for selected solver\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::reBuild(FullSquareMatrix *, int, int)
{
 fprintf(stderr,"Selected Solver does not support non-linear analysis.\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::reBuild(FullSquareMatrix *, FullSquareMatrix *, Scalar)
{
 fprintf(stderr,"Selected Solver does not support non-linear analysis.\n");
}

template<class Scalar>
int
GenSolver<Scalar>::dim()
{
 fprintf(stderr,"Selected Solver does not support dim() function\n");
 return 0;
}


template<class Scalar> 
int
GenSolver<Scalar>::numRBM()
{
 fprintf(stderr,"Selected Solver does not support Rigid Body Modes\n");
 return 0;
}

template<class Scalar> 
void
GenSolver<Scalar>::getRBMs(double *)
{
 fprintf(stderr,"Selected Solver does not support Rigid Body Modes\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::getRBMs(Vector *)
{
 fprintf(stderr,"Selected Solver does not support Rigid Body Modes\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::getRBMs(VectorSet &)
{
 fprintf(stderr,"Selected Solver does not support Rigid Body Modes\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::solve(GenVector<Scalar>  &rhs, GenVector<Scalar>  &sol)
{
 solve(rhs.data(), sol.data());
}

template<class Scalar> 
void
GenSolver<Scalar>::solve(Scalar  *, Scalar *)
{
 fprintf(stderr,"Selected Solver does not support a Scalar*.\n");
}

// This function does nothing, it is supposed to do nothing.
// When a solver does not have a factoring step, i.e. pcg, bcg
// or frontal, factoring is not used.

template<class Scalar> 
void
GenSolver<Scalar>::factor()
{
 fprintf(stderr,"Selected Solver does not implement factor() .\n"); 
}

template<class Scalar>
void
GenSolver<Scalar>::parallelFactor()
{
 factor();
}


template<class Scalar> 
void
GenSolver<Scalar>::reBuildGeometricRbms(GeomState *)
{
 fprintf(stderr,"Selected Solver does not support reBuilding Geomtric Rbms\n");
}

template<class Scalar> 
void
GenSolver<Scalar>::clean_up()
{
}

template<class Scalar> 
void 
GenSolver<Scalar>::addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier)
{
 fprintf(stderr,"Selected Solver does not support addBoeing\n");
}

template<class Scalar>
void
GenSolver<Scalar>::addone(Scalar d, int dofi, int dofj)
{
  fprintf(stderr,"Selected Solver does not support addone(...)\n");
}

template<class Scalar>
Scalar
GenSolver<Scalar>::getone(int dofi, int dofj)
{
  fprintf(stderr,"Selected Solver does not support getone(...)\n");
  return 0;
}

template<class Scalar>
void
GenSolver<Scalar>::unify(FSCommunicator *)
{
  fprintf(stderr,"Selected Solver does not support unify()\n");
}

template<class Scalar>
void
GenSolver<Scalar>::add(Scalar *d)
{
  fprintf(stderr,"Selected Solver does not support add(double *)\n");
}

template<class Scalar>
Scalar*
GenSolver<Scalar>::getData()
{
  fprintf(stderr,"Selected Solver does not support getData()\n");
  return 0;
}

