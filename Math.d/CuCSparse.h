#ifndef CUC_SPARSE_H_
#define CUC_SPARSE_H_

#include <Math.d/SparseMatrix.h>
#include <Utils.d/MyComplex.h>

class ControlLawInfo;

class LMPCons;

// Constrained to Unconstrained Sparse Matrix (For a rectangular matrix
// that needs to be stored in sparse format)
template<class Scalar>
class GenCuCSparse : public SparseData, public GenSparseMatrix<Scalar> {
    Scalar *Kuc;
    int myKuc;
  public:
    GenCuCSparse(Connectivity *con, DofSetArray *dsa, int *bc);
    GenCuCSparse(Connectivity *con, DofSetArray *dsa, DofSetArray *c_dsa);
    GenCuCSparse(Connectivity *con, DofSetArray *dsa, int *glbmap, int *glimap);

    GenCuCSparse(LMPCons **mpc, int numMPC, DofSetArray *c_dsa);
    GenCuCSparse(int numInterface, int *glbmap, int numRBM, Scalar *rbms, int lda = -1);
    GenCuCSparse(int, int, int *, int *, Scalar*);
    virtual ~GenCuCSparse();

    double getMemoryUsed();

    void add(FullSquareMatrix &mel, int *dofs);
    void addImaginary(FullSquareMatrix &kel, int *dofs);
    void add(FullSquareMatrixC &mel, int *dofs);
    void addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier);
    void addDiscreteMass(int dof, Scalar dimass);
    void add(int row_dof, int col_dof, Scalar s);

    void mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result);
    void mult(const Scalar *rhs, Scalar *result);

    void multIdentity(Scalar **result);
    void multIdentity(Scalar *result);
    void multIdentity(Scalar **v, int start, int stop);

    void multSubtract(const GenVector<Scalar> &rhs, GenVector<Scalar> &result);
    void multSubtract(const Scalar *rhs, Scalar *result);

    void multSubAdd(Scalar *rhs, Scalar *result);

    Scalar diag(int) const { throw "GenCuCSparse::diag - 1 - should never be called"; }
    Scalar &diag(int) { throw "GenCuCSparse::diag - 2 - should never be called"; }
    void zeroAll();
    int  dim()    { return numConstrained; }
    int  neqs()   { return numConstrained; }
    int  numRow() { return neq-numConstrained; }
    int  numCol() { return numConstrained; }
    void print();
    long size();

    void negate();
    void clean_up();

    void multSub(const Scalar *rhs, Scalar *result);
    void multSub(int numRHS, Scalar **rhs, Scalar **result);
    void transposeMultAdd(const Scalar *rhs,Scalar *result);
    void multAdd(const Scalar *rhs,Scalar *result);
    void transposeMultSubtract(const Scalar *rhs,Scalar *result);
    void transposeMultSubtractClaw(const Scalar *rhs,Scalar *result, int ncl,
                                   int *clawfDofs);
    void transposeMult(const Scalar *rhs,Scalar *result);

    void doWeighting(Scalar *weight);  // PJSA 4-3-03
    void doWeighting(int *weight);

    void multAddNew(const Scalar *rhs, Scalar *result);
    void transposeMultNew(const Scalar *rhs, Scalar *result);
    void transposeMultAddNew(const Scalar *rhs, Scalar *result);
    void transposeMultSubNew(const Scalar *rhs, Scalar *result);

    void  mult(  const Scalar *rhs, Scalar *result, Scalar alpha, Scalar beta); //HB 12/15/04
    void  trMult(const Scalar *rhs, Scalar *result, Scalar alpha, Scalar beta); //HB 12/15/04
};

typedef GenCuCSparse<double> CuCSparse;
typedef GenCuCSparse<DComplex> CuCComplexSparse;

#ifdef _TEMPLATE_FIX_
#include <Math.d/CuCSparse.C>
#endif

#endif
