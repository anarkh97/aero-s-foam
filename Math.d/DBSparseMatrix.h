#ifndef DB_SPARSEMATRIX_H_
#define DB_SPARSEMATRIX_H_

#include <stdio.h>
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

   void mult(const GenVector<Scalar> &, GenVector<Scalar> & ); //matrix-vector multiply
   void mult(const GenVector<Scalar> &rhs, Scalar *result); 
   void mult(const Scalar *, Scalar *); // matrix-vector multiply
   void multAdd(const Scalar *, Scalar *); // matrix-vector multiply
   void multcomplex(const DComplex *, DComplex *);  // matrix-complex vector multiply

   void multDiag(const Scalar *x, Scalar *b);
   void multDiag(int numRHS, const Scalar **x, Scalar **b);

   double getMemoryUsed();

   void transposeMult(const GenVector<Scalar> & rhs, GenVector<Scalar> & result);
   Scalar diag(int dof) const;      // returns diagonal value of row dof in matrix
   Scalar &diag(int dof);
   void add(FullSquareMatrix &, int *dofs);
   void add(FullSquareMatrixC &, int *dofs);
   void addImaginary(FullSquareMatrix &, int *dofs);
   void add(GenFullM<Scalar> &knd, int fRow, int fCol);
   void addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier);
   void addDiscreteMass(int dof, Scalar diMass);
   void add(int, int, Scalar);
   void zeroAll();
   void makeIdentity();
   int  dim() { return numUncon; }
   int  neqs() { return numUncon; }
   long size();
   void print(char *fileName);
   void print();
   void print1(int dof,FILE *fid);
   void invertDiag();
   void deleteMemory() {delete [] unonz; unonz=0; }
   void unify(FSCommunicator *communicator);
   void clean_up();   
   int  begin(int i)  { return xunonz[i]; }
   int  end(int i)    { return xunonz[i+1]; }
   void symmetricScaling();
   void applyScaling(Scalar *v);
   GenFullM<Scalar> *getFullMatrix();
   int*   getFirstDof() { cerr << "int*  GenDBSparseMatrix::getFirstDof() called" << endl; firstdof = new int[1]; firstdof[0]=0; return firstdof; }
   int getBlockSize() {cerr << "dim() = " << dim() << endl; return dim();}
};

typedef GenDBSparseMatrix<double> DBSparseMatrix;
typedef GenDBSparseMatrix<DComplex> DBComplexSparseMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/DBSparseMatrix.C>
#endif

#endif
