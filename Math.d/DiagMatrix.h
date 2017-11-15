#ifndef _DIAGMATRIX_H_
#define _DIAGMATRIX_H_

#include <Utils.d/dofset.h>
#include <Solvers.d/Solver.h>
#include <Math.d/SparseMatrix.h>

template<class Scalar>
class GenDiagMatrix : public GenSparseMatrix<Scalar>, public GenSolver<Scalar>
{
     int neq;
     Scalar *v;
     DofSetArray *dsa;

   public:
     GenDiagMatrix(DofSetArray *_dsa);
     virtual ~GenDiagMatrix();

     void   add(FullSquareMatrix &, int *) override;
     void   add(GenFullM<Scalar> &knd, int fRow, int fCol) override;
     void   add(GenAssembledFullM<Scalar> &kel, int *dofs) override;
     void   addDiscreteMass(int dof, Scalar mass) override;

     void   zeroAll() override { for(int i=0; i<neq; ++i) v[i] = 0.0; }
     int    dim() const override { return neq; }

     int neqs() const override { return neq; }
     int numCol() const override { return neq; }
     double getSolutionTime() override { return 0;}//JFD to be done
     long size() const override { return neq*sizeof(Scalar); }
     void addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier) override;
     void solve(const Scalar *rhs, Scalar *sol) override;
     Scalar diag(int d) const override;
     Scalar &diag(int d) override;
     void reSolve(Scalar*) override;
     void mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const override;
     void mult(const Scalar *rhs, Scalar *result) const override;
     void squareRootMult(Scalar *result) override;
     void inverseSquareRootMult(Scalar *result) override;
     void factor() override;

};

typedef GenDiagMatrix<double> DiagMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/DiagMatrix.C>
#endif

#endif
