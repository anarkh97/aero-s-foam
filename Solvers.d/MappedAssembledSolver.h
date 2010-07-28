/*
 * MappedAssembledSolver.h
 *
 *  Created on: Apr 15, 2009
 *      Author: michel
 */

#ifndef MAPPEDASSEMBLEDSOLVER_H_
#define MAPPEDASSEMBLEDSOLVER_H_
#include <iostream>

#include <Utils.d/resize_array.h>
#include <Math.d/FullMatrix.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Math.d/Vector.h>

struct DOFMap {
  int ndofs;
  int *dofs;
  double *coefs;
};

DOFMap *getDofMaps(int size);

class CoordinateMap
{
    long step;
    ResizeArray<long> flag;
    ResizeArray<int> subScript;

  protected:
    DOFMap *dofMaps;
    int numMappedEqs;
    DOFMap *eqMaps;

    template <class Scalar, class MatrixType >
    MatrixType map(MatrixType &mat, int *dofs, ResizeArray<int> &mappedDofs,
                   ResizeArray<Scalar> &m) {
      // TODO this function is not thread safe
      // Create the set of DOFs that are used.
      step += 1;
      int nMappedDOFs = 0;
      int ndofs = mat.dim();
      // Whether the mapping is a one to one renumbering only
      bool isRenumbering = true;
      for(int i = 0; i < mat.dim(); ++i) {
        if(dofMaps[dofs[i]].ndofs != 1 || flag[dofMaps[dofs[i]].dofs[0]] == step ||
            dofMaps[dofs[i]].coefs[0] != 1.0)
          isRenumbering = false;
        //std::cerr << "mapping " << dofs[i] << " to" ;
        for(int j = 0; j < dofMaps[dofs[i]].ndofs; ++j) {
          //std::cerr << " " << dofMaps[dofs[i]].dofs[j]
          //          << " with flag " << flag[dofMaps[dofs[i]].dofs[j]];
          if(flag[dofMaps[dofs[i]].dofs[j]] != step) {
            flag[dofMaps[dofs[i]].dofs[j]] = step;
            mappedDofs[nMappedDOFs] = dofMaps[dofs[i]].dofs[j];
            subScript[dofMaps[dofs[i]].dofs[j]] = nMappedDOFs;
            nMappedDOFs++;
          }
        }
        //std::cerr << std::endl;
      }

      // If we only have a renumbering, we are done and can return;
      if(isRenumbering)
        return mat;

      // If we have a more complex mapping, we do the transformation of the matrix
      // The intermediate matrix needs to be nMappedDOFs x ndofs. Allocate memory for it.
      Scalar *dataSpace = (Scalar *) dbg_alloca(nMappedDOFs*mat.dim()*sizeof(Scalar));
      GenStackFSFullMatrix<Scalar> mTmp(nMappedDOFs, ndofs, dataSpace);
      // Multiply by the mapping matrix on the left
      mTmp = (Scalar) 0;
      // mTmp = M^T * mat
      // mTmp[i][j] = sum_k ( M[k][i] * mat[k][j] )
      for(int k = 0; k < ndofs; ++k) {
        DOFMap &dMapK = dofMaps[dofs[k]];
        for(int iptr = 0; iptr < dMapK.ndofs; ++iptr) {
          int i = subScript[dMapK.dofs[iptr]];
          for(int j = 0; j < ndofs; ++j)
            mTmp[i][j] += dMapK.coefs[iptr]*mat[k][j];
        }
      }

      // Multiply the mapping on the right.
      // insure the matrix we build is big enough.
      m[nMappedDOFs*nMappedDOFs] = 0;
      MatrixType res(nMappedDOFs, m.data());
      res.zero();
      // res = mTmp*M
      // res[i][j] =  sum_k ( mTmp[i][k] * M[k][j] )
      for(int k = 0; k < ndofs; ++k) {
        DOFMap &dMapK = dofMaps[dofs[k]];
        for(int jptr = 0; jptr < dMapK.ndofs; ++jptr) {
          int j = subScript[dMapK.dofs[jptr]];
          for(int i = 0; i < nMappedDOFs; ++i)
            res[i][j] += mTmp[i][k]*dMapK.coefs[jptr];
        }
      }

      return res;
    }

  public:
    CoordinateMap(int nEq, DOFMap *baseMap, int nMappedEq, DOFMap *eqMap):
      dofMaps(baseMap), numMappedEqs(nMappedEq), eqMaps(eqMap),
      step(0), flag(0), subScript(0) {}
    ~CoordinateMap() {}

   template <class VectorType>
   double norm(VectorType &v) {
     VectorType v2(numMappedEqs);
     v2.zero();
     for(int i = 0; i < numMappedEqs; ++i) {
       for(int j = 0; j < eqMaps[i].ndofs; ++j) {
         v2[eqMaps[i].dofs[j]] += eqMaps[i].coefs[j]*v[i];
       }
     }
     return v2.norm();
   }

};


template<class BaseSolver, class Scalar, class Map = CoordinateMap>
class MappedAssembledSolver : public BaseSolver, public Map
{
    ResizeArray<int> mappedDofs;
    ResizeArray<double> dd;
    ResizeArray<std::complex<double> > dz;
    ResizeArray<Scalar> ds;

  public:
    template <class BaseArgs>
    MappedAssembledSolver(BaseArgs &ba, int nEq, DOFMap *baseMap, int nMappedEq, DOFMap *eqMap):
      BaseSolver(ba), Map(nEq, baseMap, nMappedEq, eqMap), mappedDofs(0), dd(0), dz(0), ds(0) {}
    ~MappedAssembledSolver() {}

    // Note: None of the add routines are thread safe
    void add(FullSquareMatrix &mat, int *dofs) {
      FullSquareMatrix m = Map::map(mat, dofs, mappedDofs, dd);
      m.setMyval(0);
      BaseSolver::add(m, mappedDofs.data());
    }

    void addImaginary(FullSquareMatrix &mat, int *dofs) {
      FullSquareMatrix m = Map::map(mat, dofs, mappedDofs, dd);
      BaseSolver::addImaginary(m, mappedDofs.data());
    }

    void add(FullSquareMatrixC &mat, int *dofs) {
      FullSquareMatrixC m = Map::map(mat, dofs, mappedDofs, dz);
      BaseSolver::add(m, mappedDofs.data());
    }

    void add(GenFullM<Scalar> &mat, int *dofs) {
      GenFullM<Scalar> m = Map::map(mat, dofs, mappedDofs, dz);
      BaseSolver::add(m, mappedDofs.data());
    }

    // TODO revisit
    void add(GenFullM<Scalar> &mat, int fi, int fj) {
      cerr << "MappedAssembledSolver::add(GenFullM<Scalar> &, int, int) is not implemented\n";
    }

    // TODO revisit
    void add(GenAssembledFullM<Scalar> &mat, int *dofs) {
      cerr << "MappedAssembledSolver::add(GenAssembledFullM<Scalar> &, int *) is not implemented\n";
    }

    void addDiscreteMass(int dof, Scalar s) {
      int dofs[1] = { dof };
      double d[1] = { ScalarTypes::Real(s) };
      FullSquareMatrix mat(1, d);
      FullSquareMatrix m = Map::map(mat, dofs, mappedDofs, dd);
      BaseSolver::add(m, mappedDofs.data());
    }

    void add(int dofi, int dofj, Scalar s) {
      if(dofi == dofj) {
        addDiscreteMass(dofi, s);
      }
      else {
        int dofs[2] = { dofi, dofj };
        double d[4] = { 0, ScalarTypes::Real(s), 0, 0 };
        FullSquareMatrix mat(2, d);
        FullSquareMatrix m = Map::map(mat, dofs, mappedDofs, dd);
        BaseSolver::add(m, mappedDofs.data());
      }
    }

    void factor() {
      BaseSolver::factor();
    } 

    void reSolve(Scalar *rhs);

    void solve(GenVector<Scalar> &rhs, GenVector<Scalar> &solution) { 
      solution = rhs;
      reSolve(solution.data());
    }

};

template<class BaseSolver, class Scalar, class Map>
void MappedAssembledSolver<BaseSolver, Scalar, Map>::reSolve(Scalar *s)
{
  GenVector<Scalar> s2(BaseSolver::neqs());
  s2.zero();
  for(int i = 0; i < Map::numMappedEqs; ++i) {
    for(int j = 0; j < Map::eqMaps[i].ndofs; ++j) {
      s2[Map::eqMaps[i].dofs[j]] += Map::eqMaps[i].coefs[j]*s[i];
    }
  }
  BaseSolver::reSolve(s2.data());
  for(int i = 0; i < Map::numMappedEqs; ++i) {
    s[i] = 0;
    for(int j = 0; j < Map::eqMaps[i].ndofs; ++j)
      s[i] += s2[Map::eqMaps[i].dofs[j]]*Map::eqMaps[i].coefs[j];
  }
}

#endif /* MAPPEDASSEMBLEDSOLVER_H_ */
