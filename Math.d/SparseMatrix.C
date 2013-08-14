#include <cstdio>

template<class Scalar>
GenSparseMatrix<Scalar>::~GenSparseMatrix() { } // empty destructor

template<class Scalar> 
void
GenSparseMatrix<Scalar>::clean_up()
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::clean_up() not implemented\n");
}

template<class Scalar> 
double
GenSparseMatrix<Scalar>::getMemoryUsed()
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::getMemoryUsed() not implemented\n");
 return 0;
}

template<class Scalar> 
int
GenSparseMatrix<Scalar>::numRow()
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::numRow() not implemented\n");
 return 0;
}

template<class Scalar> 
int
GenSparseMatrix<Scalar>::numCol()
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::numCol() not implemented\n");
 return 0;
}

template<class Scalar>
double
GenSparseMatrix<Scalar>::norm()
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::norm() not implemented\n");
 return 1.0;
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::invertDiag()
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::invertDiag() not implemented\n");
}


template<class Scalar>
void
GenSparseMatrix<Scalar>::printSparse(const std::string& filename)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::printSparse() not implemented\n");
}


template<class Scalar> 
void
GenSparseMatrix<Scalar>::addImaginary(FullSquareMatrix &, int *dofs)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::addImaginary(FullSquareMatrix &, int *dofs) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::add(FullSquareMatrixC &, int *dofs)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::add(FullSquareMatrixC &, int *dofs) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::add(GenFullM<Scalar> &knd, int fRow, int fCol)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::add(GenFullM<Scalar> &knd, int fRow, int fCol) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::add(GenAssembledFullM<Scalar> &knd, int *dofs)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::add(GenAssembledFullM<Scalar> &kel, int *dofs) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::addDiscreteMass(int dof, Scalar s)
{
  int dofs[1] = { dof };
  double reS, imS;
  if((reS = ScalarTypes::Real(s)) != 0) {
    double d[1] = { reS };
    FullSquareMatrix mat(1, d);
    add(mat,dofs);
  }
  if((imS = ScalarTypes::Imag(s)) != 0) {
    double d[1] = { imS };
    FullSquareMatrix mat(1, d);
    addImaginary(mat,dofs);
  }
}

template<class Scalar>
void
GenSparseMatrix<Scalar>::add(int dofi, int dofj, Scalar s)
{
  if(dofi == dofj) {
    addDiscreteMass(dofi, s);
  }
  else {
    double reS, imS;
    int dofs[2] = { dofi, dofj };
    if((reS = ScalarTypes::Real(s)) != 0) {
      double d[4] = { 0, reS, reS, 0 };
      FullSquareMatrix mat(2, d);
      add(mat,dofs);
    }
    if((imS = ScalarTypes::Imag(s)) != 0) {
      double d[4] = { 0, imS, imS, 0 };
      FullSquareMatrix mat(2, d);
      addImaginary(mat,dofs);
    }
  }
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) not implemented\n");
}

template<class Scalar>
void
GenSparseMatrix<Scalar>::mult(const GenVector<Scalar> &rhs, Scalar *result) 
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::mult(const GenVector<Scalar> &rhs, Scalar *result) not implemented\n");
}


template<class Scalar> 
void
GenSparseMatrix<Scalar>::mult(const Scalar *rhs, Scalar *result)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::mult(const Scalar *rhs, Scalar *result) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::multAdd(const GenVector<Scalar> &rhs, GenVector<Scalar> &result)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::multAdd(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::multAdd(const Scalar *rhs, Scalar *result)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::multAdd(const Scalar *rhs, Scalar *result) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::multSubtract(const GenVector<Scalar> &rhs, GenVector<Scalar> &result)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::multSubtract(const GenVector<Scalar> &rhs, GenVector<Scalar> &result) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::multSubtract(const Scalar *rhs, Scalar *result)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::multSubtract(const Scalar *rhs, Scalar *result) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::squareRootMult(Scalar *result)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::squareRootMult(Scalar *result) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::inverseSquareRootMult(Scalar *result)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::inverseSquareRootMult(Scalar *result) not implemented\n");
}

template<class Scalar>
void 
GenSparseMatrix<Scalar>::transposeMult(const GenVector<Scalar> & rhs, GenVector<Scalar> & result)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::transposeMult((const GenVector<Scalar> &, GenVector<Scalar> &) not implemented\n");
}

template<class Scalar>
void
GenSparseMatrix<Scalar>::transposeMult(const Scalar *, Scalar *)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::transposeMult(const Scalar *, Scalar *) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::transposeMultAdd(const Scalar *, Scalar *)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::transposeMultAdd(const Scalar *, Scalar *) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::transposeMultSubtract(const Scalar *, Scalar *)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::transposeMultSubtract(const Scalar *, Scalar *) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::transposeMultSubtractClaw(const Scalar *, Scalar *, int, int *)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::transposeMultSubtractClaw(const Scalar *, Scalar *, int, int *) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::multSub(const Scalar *, Scalar *r)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::multSub(const Scalar *, Scalar *) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::multSub(int, Scalar **, Scalar **)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::multSub(int, Scalar **, Scalar **) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::multDiag(const Scalar *, Scalar *)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::multDiag(const Scalar *, Scalar *) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::multIdentity(Scalar *)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::multIdentity(Scalar *) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::multIdentity(Scalar **)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::multIdentity(Scalar **) not implemented\n");
}

template<class Scalar> 
void
GenSparseMatrix<Scalar>::multIdentity(Scalar **v, int start, int stop)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::multIdentity(Scalar **v, int start, int stop) not implemented\n");
}

template<class Scalar>
GenFullM<Scalar> *
GenSparseMatrix<Scalar>::getFullMatrix()
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::getFullMatrix() not implemented\n");
 return NULL;
}

template<class Scalar>
int* GenSparseMatrix<Scalar>::getFirstDof()
{
 if(!firstdof) {
    firstdof = new int[1];
    firstdof[0]=0;
  }
  return firstdof;
  cerr << "GenSparseMatrix<Scalar>::getFirstDof() called " << endl;
//fprintf(stderr,"GenSparseMatrix<Scalar>::getFirstDof() not implemented\n");
}

template<class Scalar>
int GenSparseMatrix<Scalar>::numNodes()
{
  return 1;
//fprintf(stderr,"GenSparseMatrix<Scalar>::numNodes() not implemented\n");
}

template<class Scalar>
GenFullM<Scalar>* 
GenSparseMatrix<Scalar>::getDiagMatrix(int i)
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::getDiagMatrix(int i) not implemented\n");
 return NULL;
}

template<class Scalar>
Scalar* GenSparseMatrix<Scalar>::getBlockScalarMultipliers()
{
 scalarfactors=new Scalar[1];
 scalarfactors[0] = 1;
 return scalarfactors;
}

template<class Scalar>
void GenSparseMatrix<Scalar>::setMeanSolver(GenSolver<Scalar> *prc)
{
 meansolver= prc;
}

template<class Scalar>
GenSolver<Scalar>* GenSparseMatrix<Scalar>::getMeanSolver()
{
 return meansolver;
}

template<class Scalar>
int GenSparseMatrix<Scalar>::getBlockSize()
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::getBlockSize() not implemented\n");
 return 0;
}

template<class Scalar>
void GenSparseMatrix<Scalar>::print()
{
 fprintf(stderr,"GenSparseMatrix<Scalar>::print() not implemented\n");
}

