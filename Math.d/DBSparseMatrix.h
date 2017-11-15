#ifndef DB_SPARSEMATRIX_H_
#define DB_SPARSEMATRIX_H_

#include <cstdio>
#include <iostream>
#include <Utils.d/MyComplex.h> 
#include <Math.d/SparseMatrix.h>
class FSCommunicator;

// Classical Sparse Matrix

template<class Scalar>
class GenDBSparseMatrix : public SparseData, public GenSparseMatrix<Scalar> {
   Scalar *unonz;
   int     isScaled;
   Scalar *scale;
   int *firstdof;
 public:
   // Constructor
   GenDBSparseMatrix(Connectivity *, DofSetArray *, ConstrainedDSA *c_dsa);
   GenDBSparseMatrix(Connectivity *, DofSetArray *, int* rcn);
   GenDBSparseMatrix(Connectivity *, EqNumberer *);
   virtual ~GenDBSparseMatrix();

   void mult(const GenVector<Scalar> &, GenVector<Scalar> & ) const override; //matrix-vector multiply
   void mult(const GenVector<Scalar> &rhs, Scalar *result) const override;
   void mult(const Scalar *, Scalar *) const override; // matrix-vector multiply
   void multAdd(const Scalar *, Scalar *) const override; // matrix-vector multiply
   void multcomplex(const DComplex *, DComplex *) const;  // matrix-complex vector multiply

   void multDiag(const Scalar *x, Scalar *b) const override;
   void multDiag(int numRHS, const Scalar **x, Scalar **b) const;

   double getMemoryUsed() const override;

   void transposeMult(const GenVector<Scalar> & rhs, GenVector<Scalar> & result) const override;
   void transposeMult(const Scalar *, Scalar *) const override;
   Scalar diag(int dof) const override;      // returns diagonal value of row dof in matrix
   Scalar &diag(int dof) override;
   void add(FullSquareMatrix &, int *dofs) override;
   void add(FullSquareMatrixC &, int *dofs) override;
   void addImaginary(FullSquareMatrix &, int *dofs) override;
   void add(GenFullM<Scalar> &knd, int fRow, int fCol) override;
   void addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier);
   void addDiscreteMass(int dof, Scalar diMass) override;
   void add(int, int, Scalar) override;
   void zeroAll() override;
   void makeIdentity();
   int  dim() const override { return numUncon; }
   int  neqs() const override { return numUncon; }
   int  numRow() const override { return numUncon; }
   long size() const;
   void print(char *fileName);
   void print() override;
   void print1(int dof,FILE *fid);
   void invertDiag() override;
   void deleteMemory() {delete [] unonz; unonz=0; }
   void unify(FSCommunicator *communicator);
   void clean_up() override;
   int  begin(int i)  { return xunonz[i]; }
   int  end(int i)    { return xunonz[i+1]; }
   void symmetricScaling();
   void applyScaling(Scalar *v);
   GenFullM<Scalar> *getFullMatrix() override;
   int*   getFirstDof() override
   { std::cerr << "int*  GenDBSparseMatrix::getFirstDof() called" << std::endl; firstdof = new int[1]; firstdof[0]=0; return firstdof; }
   int getBlockSize() override
   {std::cerr << "dim() = " << dim() << std::endl; return dim();}
};

typedef GenDBSparseMatrix<double> DBSparseMatrix;
typedef GenDBSparseMatrix<DComplex> DBComplexSparseMatrix;

#endif
