#ifndef _FETI_MPC_SPARSE_H_
#define _FETI_MPC_SPARSE_H_

#include <cstdio>
#include <Math.d/SparseMatrix.h>
#include <Driver.d/Mpc.h>
#include <Utils.d/MyComplex.h>

// Container class that holds the MPC equation information

template<class Scalar>
class GenMpcSparse : public GenSparseMatrix<Scalar> 
{
    SubLMPCons<Scalar> **mpc;
    int NumRow;
    int NumCol;
    DofSetArray *dsa;
    DofSetArray *DSA; 
    int *wetInterfaceMap;
    int mpcOffset;
   
  public:
    GenMpcSparse(int numMPC, SubLMPCons<Scalar> **_mpc, DofSetArray *_dsa);
    GenMpcSparse(int numMPC, SubLMPCons<Scalar> **_mpc, DofSetArray *_dsa, DofSetArray *_DSA, 
                 int *_wetInterfaceMap, int mpcOffset);
    virtual ~GenMpcSparse();

    double getMemoryUsed();

    void add(const double *const *kel, int *dofs, int kndof);
    void add(FullSquareMatrix &mel, int *dofs);
    void add(FullM&, int, int);
    void add(GenAssembledFullM<Scalar> &, int*);

    void addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier);

    void mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result);
    void mult(const Scalar *rhs, Scalar *result);
    
    void multIdentity(int dof, Scalar *result);
    void multIdentity(Scalar **result);
    void multIdentity(Scalar *result);
    void multSubtract(const GenVector<Scalar> &rhs, GenVector<Scalar> &result);
    void multSubtract(const Scalar *rhs, Scalar *result);

    void zeroAll();
    int  dim()    { return 0; }
    int  neqs()    { return 0; }
    int  numRow() { return NumRow; }
    int  numCol() { return NumCol; }

    Scalar diag(int) const { throw "GenMpcSparse::diag - 1 - should never be called"; }
    Scalar &diag(int) { throw "GenMpcSparse::diag - 2 - should never be called"; }
    long size() const { return 0; }
    void print(FILE *file=stderr, const char* msg="A");
    void negate();

    void multSub(const Scalar *rhs, Scalar *result);
    void multSub(int numRHS, Scalar **rhs, Scalar **result);
    void multSubWI(const Scalar *rhs, Scalar *result);
    void transposeMultAdd(const Scalar *rhs, Scalar *result);
    void multAdd(const Scalar *rhs, Scalar *result);
    void transposeMultSubtract(const Scalar *rhs, Scalar *result);
    void transposeMultSubtractWI(const Scalar *rhs, Scalar *wi_result);
    void transposeMult(const Scalar *rhs, Scalar *result);

  private:
    GenMpcSparse<Scalar>& operator = (const GenMpcSparse<Scalar> &);
    GenMpcSparse<Scalar>(const GenMpcSparse<Scalar> &);
};

typedef GenMpcSparse<double> MpcSparse;
typedef GenMpcSparse<DComplex> MpcSparseC;

#ifdef _TEMPLATE_FIX_
#include <Math.d/MpcSparse.C>
#endif


#endif
