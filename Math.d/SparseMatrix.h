#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#include <Utils.d/dofset.h>

#include <Math.d/matrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/MyComplex.h>

#include <string>

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
typedef GenVector<DComplex> ComplexVector;
class Connectivity;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenAssembledFullM;
typedef GenAssembledFullM<DComplex> AssembledFullMC;
template <class Scalar> class GenFullM;
typedef GenFullM<DComplex> FullMC;
template <class Scalar> class GenSolver;
template <class Scalar> class DistrBlockVector;


template<class Scalar>
class GenSMV {
 GenSparseMatrix<Scalar> &M;
 GenVector<Scalar>       &v;
 public:
   GenSMV(GenSparseMatrix<Scalar> &_M, GenVector<Scalar> &_v) : M(_M), v(_v) {}
};

template<class Scalar>
class GenSparseMatrix {
        GenSolver<Scalar>* meansolver;
        int *firstdof;
        Scalar* scalarfactors;
   public:       
        virtual void zeroAll() = 0;
        virtual int dim() const = 0;
        virtual double norm() const;
        virtual void clean_up();
        virtual double getMemoryUsed() const;
        virtual int  numRow() const;
        virtual int  numCol() const;
        virtual Scalar diag(int dof) const = 0;
        virtual Scalar &diag(int dof) = 0;
        virtual void invertDiag();
        virtual void printSparse(const std::string& filename);

        virtual void add(const FullSquareMatrix &, const int *dofs) = 0;
        virtual void addImaginary(const FullSquareMatrix &, const int *dofs);
        virtual void add(const FullSquareMatrixC &, const int *dofs);
        virtual void add(const GenFullM<Scalar> &knd, int fRow, int fCol);
        virtual void add(const GenAssembledFullM<Scalar> &kel, const int *dofs) ;
        virtual void addDiscreteMass(int dof, Scalar mass);
        virtual void add(int row_dof, int col_dof, Scalar s);

        virtual void mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const;
        virtual void mult(const GenVector<Scalar> &rhs, Scalar *result) const;
        virtual void mult(const Scalar *rhs, Scalar *result) const;
        virtual void multAdd(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const;
        virtual void multAdd(const Scalar *rhs, Scalar *result) const;
        virtual void multSubtract(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) const;
        virtual void multSubtract(const Scalar *rhs, Scalar *result) const;
        virtual void transposeMult(const GenVector<Scalar> & rhs, GenVector<Scalar> & result) const;
        virtual void transposeMult(const Scalar *, Scalar *) const;
        virtual void transposeMultAdd(const Scalar *, Scalar *) const;
        virtual void transposeMultSubtract(const Scalar *, Scalar *) const;
        virtual void transposeMultSubtractClaw(const Scalar *, Scalar *, int, int *) const;
        virtual void multSub(const Scalar *, Scalar *) const;
        virtual void multSub(int, const Scalar **, Scalar **) const;
        virtual void multDiag(const Scalar *, Scalar *) const;
        /** \brief Compute \f$ R = R + M\f$.
         *
         * @param R Pointer to result matrix to modify in dense row-wise matrix storage.
         */
        virtual void multIdentity(Scalar *R) const;
        virtual void multIdentity(Scalar **) const;
        virtual void multIdentity(Scalar **v, int start, int stop) const;
        virtual void squareRootMult(Scalar * result);
        virtual void inverseSquareRootMult(Scalar * result);
        virtual GenFullM<Scalar> * getFullMatrix();
        virtual int* getFirstDof();
        virtual int numNodes() const;
        virtual GenFullM<Scalar>* getDiagMatrix(int i);
        virtual ~GenSparseMatrix() = 0;
        virtual Scalar* getBlockScalarMultipliers();
        virtual void setMeanSolver(GenSolver<Scalar> *prc);
        virtual GenSolver<Scalar>* getMeanSolver();
        virtual int getBlockSize();
        virtual int  neqs() const = 0;
        virtual void print();

        void mult(DistrBlockVector<double>&, DistrBlockVector<double>&) { }; // hack to get code to compile

        virtual void matvec(GenVector<Scalar> &rhs, GenVector<Scalar> &result) { mult(rhs, result); }
};


class LMPCons;

class SparseData {
 protected:
  int *unconstrNum;
  int *constrndNum;
  int *xunonz;
  int *rowu;
  int *colu;
  int numConstrained;
  int numUncon;
  int neq;
  int myMem;
  int myMem_rowu;
 public: 
    SparseData();
    // Constructors for data structures of type CuCSparse and CuCComplexSparse
    SparseData(Connectivity *con, DofSetArray *dsa, const int *bc);
    SparseData(Connectivity *con, DofSetArray *dsa, DofSetArray *c_dsa);
    SparseData(Connectivity *con, DofSetArray *dsa, const int *glbmap, const int *glimap);

    // Constructors for data structures of type DBSparseMatrix 
    // and ComplexDBSparseMatrix and Spooles/Mumps
    SparseData(DofSetArray *dsa, DofSetArray *c_dsa, Connectivity *con,
       	       int expand = 0, int make_colu = 0, bool unsym = false);
    SparseData(EqNumberer *dsa, Connectivity *con, const int* rcn, int expand = 0, int make_colu = 0);

    SparseData(DofSetArray *_dsa, const int *glInternalMap,
               Connectivity *cn, int expand);

    // This constructor is for the Esmond sparse solver (BLKSparseMatrix)
    SparseData(Connectivity *cn, EqNumberer *eqn, double trbm, int expand=1);

    // MLX This constructor is for the Padma sparse solver
    SparseData(EqNumberer *eqn, Connectivity *cn, double trbm);

    // KHP: for storing mpcs
    SparseData(LMPCons **mpc, int numColumns, DofSetArray *cdsa);
    SparseData(int num, const int *xyzCount, int *xyzList);

    // for storing G for use in FETI-DP
    SparseData(int numInterface,
               const int *glbmap, int numModes, int ldm);

    virtual ~SparseData();
    void clean_up();

  private:
    void initialize();

};

typedef GenSMV<double> SMV; 
typedef GenSparseMatrix<double> SparseMatrix; 
typedef GenSparseMatrix<DComplex> ComplexSparseMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/SparseMatrix.C>
#endif

#endif
