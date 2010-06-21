#include <Element.d/Rigid.d/RigidBeam.h>
#include <Element.d/Joint.d/ConstantDistanceConstraint.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType1.h>
#include <Element.d/Joint.d/DotConstraintType1.h>
#include <Element.d/Joint.d/ParallelAxesConstraintType2.h>
#include <Element.d/Joint.d/LinearConstraintType1.h>
#include <Element.d/Joint.d/DotConstraintType2.h>

#ifdef USE_EIGEN2
#include <Math.d/rref.h>
#include <Eigen/Core>
using namespace Eigen;
#include <map>
#include <vector>
#endif

RigidBeam::RigidBeam(int* _nn)
{
  nSubElems = 4; 
  subElems = new Element * [nSubElems];
  int indices[2] = { 0, 1 };

  subElems[0] = new ConstantDistanceConstraint(indices);
  subElems[1] = new ParallelAxesConstraintType2(indices);
  subElems[2] = new DotConstraintType1(indices, 2, 1);
  subElems[3] = new ParallelAxesConstraintType1(indices);

  initialize(2, _nn);
}

LMPCons **
RigidBeam::getMPCs()
{
  int numLMPC = getNumMPCs();
  LMPCons **lmpc = SuperElement::getMPCs();
#ifdef USE_EIGEN2
  using std::map;
  using std::pair;
  using std::vector;
  // Create a unique integer ID for every DOF involved in an MPC and count
  // in how many MPC each DOF appears.
  int nID = 0;
  map<pair<int,int>, int> dofID;

  int masterNode = nn[1];
  // first the slave terms
  for(int i = 0; i < numLMPC; ++i) {
    // flush the MPC from any zero terms
    lmpc[i]->removeNullTerms();
    for(int j = 0; j < lmpc[i]->nterms; ++j) {
      if(lmpc[i]->terms[j].nnum == masterNode) continue;
      pair<int,int> p(lmpc[i]->terms[j].nnum, lmpc[i]->terms[j].dofnum);
      map<pair<int,int>, int>::iterator it = dofID.find(p);
      if(it == dofID.end()) {
        dofID[p] = nID++;
      }
    }
  }
  // now the master terms
  for(int i = 0; i < numLMPC; ++i) {
    for(int j = 0; j < lmpc[i]->nterms; ++j) {
      if(lmpc[i]->terms[j].nnum != masterNode) continue;
      pair<int,int> p(lmpc[i]->terms[j].nnum, lmpc[i]->terms[j].dofnum);
      map<pair<int,int>, int>::iterator it = dofID.find(p);
      if(it == dofID.end()) {
        dofID[p] = nID++;
      }
    }
  }

  // Obtain the MPC to DOF connectivity
  SetAccess<LMPCons> lmpcAccess(numLMPC, lmpc, dofID);
  Connectivity lmpcToDof(lmpcAccess);
  Connectivity *dofToLMPC = lmpcToDof.reverse();

  vector<int> *term2col = new vector<int>[numLMPC];
  vector<pair<int,int> > col2pair(dofToLMPC->csize());
  for(int i = 0; i < numLMPC; ++i) {
    for(int j = 0; j < lmpc[i]->nterms; ++j) {
      pair<int, int> p(lmpc[i]->terms[j].nnum, lmpc[i]->terms[j].dofnum);
      map<pair<int,int>, int>::iterator it = dofID.find(p);
      term2col[i].push_back(it->second); col2pair[it->second] = p;
    }
  }

  // copy lmpc coefficients into a dense matrix
  Matrix<double,Dynamic,Dynamic> c(numLMPC, dofToLMPC->csize());
  c.setZero();
  for(int i = 0; i < numLMPC; ++i) {
    for(int j = 0; j < lmpc[i]->nterms; ++j) {
      c(i,term2col[i][j]) = lmpc[i]->terms[j].coef.r_value;
    }
  }
  delete [] term2col;

  // compute the reduced row echelon form, without column pivoting
  double t = -getTime();
  int *rowmap = new int[c.rows()];
  for(int i = 0; i < c.rows(); ++i) rowmap[i] = i;
  ToReducedRowEchelonForm<double, Matrix<double,Dynamic,Dynamic> >(c, rowmap);

  // copy the coefficients of the dense matrix back into the lmpc data structure 
  for(int i = 0; i < numLMPC; ++i) {
    lmpc[rowmap[i]]->terms.clear();
    lmpc[rowmap[i]]->nterms = 0;
    for(int j = i; j < c.cols(); ++j) {
      if(j > i && j < numLMPC) continue;
      if(std::abs<double>(c(i,j)) > 1e-9) { // std::numeric_limits<double>::epsilon()
        LMPCTerm t(col2pair[j].first, col2pair[j].second, c(i,j));
        lmpc[rowmap[i]]->terms.push_back(t);
        lmpc[rowmap[i]]->nterms++;
      }
    }
  }

  delete [] rowmap;
#endif
  return lmpc;
}
