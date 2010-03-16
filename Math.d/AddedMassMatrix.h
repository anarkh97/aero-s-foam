#ifndef __ADDEDMASSMATRIX_H__
#define __ADDEDMASSMATRIX_H__

#include <Math.d/DBSparseMatrix.h>
#include <Solvers.d/Solver.h>

template<class Scalar, class ConstraintOperator>
class AddedMassMatrix : public GenDBSparseMatrix<Scalar>
{
  private: 
    GenSolver<Scalar> *fluidMassSolver;
    ConstraintOperator *op;
    void (ConstraintOperator::*mult)(GenVector<Scalar> &, GenVector<Scalar> &);
    void (ConstraintOperator::*trMult)(GenVector<Scalar> &, GenVector<Scalar> &);
    GenVector<Scalar> &tmp;

  public:
    // Constructor
    AddedMassMatrix(Connectivity* con, DofSetArray* dsa, ConstrainedDSA* c_dsa, GenSolver<Scalar>* sol,
                    ConstraintOperator *_op, void (ConstraintOperator::*_mult)(GenVector<Scalar> &, GenVector<Scalar> &),
                    void (ConstraintOperator::*_trMult)(GenVector<Scalar> &, GenVector<Scalar> &)) 
     : fluidMassSolver(sol), op(_op), mult(_mult), trMult(_trMult), tmp(sol->neqs()), GenDBSparseMatrix<Scalar>(con, dsa, c_dsa) {}

    // Destructor
    ~AddedMassMatrix() { }

    void mult(const GenVector<Scalar> &x, GenVector<Scalar> &y) { 
      (op->*trMult)(x,tmp);
      fluidMassSolver->reSolve(tmp);
      (op->*mult)(tmp,y);
      GenDBSparseMatrix<Scalar>::multAdd(y,y);
    }
};

#endif
