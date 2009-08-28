#ifndef _DIAGMATRIX_H_
#define _DIAGMATRIX_H_

#include <Utils.d/dofset.h>
#include <Math.d/SparseMatrix.h>

template<class Scalar>
class GenDiagMatrix: public GenSparseMatrix<Scalar> , public GenSolver<Scalar>
{
     int neq;
     Scalar *v;
     DofSetArray *dsa;

   public:
     GenDiagMatrix(DofSetArray *_dsa);
     virtual ~GenDiagMatrix();

     void   add(FullSquareMatrix &, int *);
     void   add(GenFullM<Scalar> &knd, int fRow, int fCol);
     void   add(GenAssembledFullM<Scalar> &kel, int *dofs);
     void   addDiscreteMass(int dof, Scalar mass);

     void   zeroAll() { for(int i=0; i<neq; ++i) v[i] = 0.0; }
     int    dim() { return neq; }

     int neqs() { return neq; }
     double getSolutionTime(){ return 0;}//JFD to be done
     long size() { return neq*sizeof(Scalar); }
     void addBoeing(int, const int *, const int *, const double *, int *, Scalar multiplier);
     void solve(Scalar *rhs, Scalar *sol);
     Scalar diag(int d) const;
     Scalar &diag(int d);
     void reSolve(Scalar*);
     //void print();
     void mult(const GenVector<Scalar> &rhs, GenVector<Scalar> &result);
     void mult(const Scalar *rhs, Scalar *result);
     void factor() { int count=0; for(int i=0; i<neq; ++i) if(v[i] == 0.0) count++; cerr << count << " zero diagonal terms\n"; }

};

typedef GenDiagMatrix<double> DiagMatrix;

#ifdef _TEMPLATE_FIX_
#include <Math.d/DiagMatrix.C>
#endif

#endif
