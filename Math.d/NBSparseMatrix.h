#ifndef NB_SPARSEMATRIX_H_
#define NB_SPARSEMATRIX_H_

#include <Utils.d/MyComplex.h>
#include <Math.d/SparseMatrix.h>

template <class Scalar> class GenSparseMatrix;
class Connectivity;
template <class Scalar> class GenFullM;
class ConstrainedDSA;
template <class Scalar> class GenVector;
template <class Scalar> class GenFullSquareMatrix;

// Node Based Sparse Matrix

template<class Scalar>
class GenNBSparseMatrix : public GenSparseMatrix<Scalar>
{
   Connectivity *con;           // Node to node connectivity
   GenFullM<Scalar> *allM;
   int *firstDof;               // Array of the first DOF number for a node
   ConstrainedDSA *dsa;
   int *dofToNode;
   int numnodes;
   int *bc;

 public:
   GenNBSparseMatrix(Connectivity *con, ConstrainedDSA *c_dsa);
   virtual ~GenNBSparseMatrix();   
 
   void mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result ); // Matrix-Vector multiply
   void mult(const Scalar *rhs, Scalar *result );   // Matrix-Vector multiply
   void multAdd(const Scalar *rhs, Scalar *result);
   Scalar diag(int n) const;                 // returns diagonal values of matrix
   Scalar &diag(int n);                 // returns diagonal values of matrix
   void add(FullSquareMatrix &m, int *dofs);
   void add(FullSquareMatrixC &m, int *dofs);
   void addImaginary(FullSquareMatrix &m, int *dofs);
   void zeroAll();

   // returns number of unconstrained dof
   int  dim() { return dsa->size(); }
   int  neqs() { return dsa->size(); }

   // returns a pointer to a diagonal matrix of node i
   GenFullM<Scalar>* getDiagMatrix(int i);
   int*   getFirstDof() { return firstDof; }
   int    numNodes() {return numnodes; }
};

typedef GenNBSparseMatrix<double> NBSparseMatrix;
typedef GenNBSparseMatrix<DComplex> NBComplexSparseMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/NBSparseMatrix.C>
#endif


#endif
