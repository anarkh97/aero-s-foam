#include <Solvers.d/KProject.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Rbm.h>
#include <Solvers.d/Preconditioner.h>

template<class Scalar, class AnyVector, class AnyOperator>
GenNonLinPCGSolver<Scalar,AnyVector,AnyOperator>::GenNonLinPCGSolver(AnyOperator* K, int _precno, int _maxitr,
                 double _tolpcg, int _numele, Connectivity* _dofs, int _kryflg,
                 int _initflg, int _reorthoflg, int _maxVecStorage, Rbm *_rbm) :
   GenPCGSolver<Scalar,AnyVector,AnyOperator>(K, _precno, _maxitr, _tolpcg, 0, _rbm)
{
 numele        = _numele;
 allDofs       = _dofs;
 this->kryflg        = _kryflg;
 this->initflg       = _initflg;
 this->reorthoflg    = _reorthoflg;
 this->maxVecStorage = _maxVecStorage;

 // ... CALL KRYLOV CONSTRUCTOR
 int numDof = K->dim();
 if(this->kryflg)
   this->proj = new KrylovProjector<Scalar,AnyVector>(numDof,this->maxVecStorage);
 else
   this->proj = 0;
}

template<class Scalar, class AnyVector, class AnyOperator>
void
GenNonLinPCGSolver<Scalar,AnyVector,AnyOperator>::reBuild(FullSquareMatrix *kel, int, int)
{
 if(this->proj) this->proj->newKrylov();

 this->A->zeroAll();

 int iele;
 for(iele=0; iele<numele; ++iele)
   this->A->add(kel[iele],(*allDofs)[iele]);
}

template<class Scalar, class AnyVector, class AnyOperator>
void
GenNonLinPCGSolver<Scalar,AnyVector,AnyOperator>::reBuild(FullSquareMatrix *kel,FullSquareMatrix *mel,
                                                          double delta)
{
 if(this->proj) this->proj->newKrylov();

 this->A->zeroAll();

 int iele,i, j;

 double delta2 = delta*delta;

 for(iele=0; iele<numele; ++iele) {

   int dim = kel[iele].dim();

   for(i = 0; i < dim; ++i)
     for(j = 0; j < dim; ++j) {
       double m = mel[iele][i][j];
       double k = kel[iele][i][j];
       kel[iele][i][j] = delta2*k + m;
     }

   this->A->add(kel[iele], (*allDofs)[iele]);
 }

}

