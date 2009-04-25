#ifndef _DYNAMSTATESET_C_
#define _DYNAMSTATESET_C_

#include <cstdio>
#include <Pita.d/DynamStateSet.h>
#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>

using std::fprintf;

// Memory management

template <typename Scalar>
void DynamStateSet<Scalar>::init_(int vectorSize, int maxNumStates)
{
  stateSet_.reserve((maxNumStates > 0) ? maxNumStates : 0);
  vectorSize_ = (vectorSize > 0) ? vectorSize : 0;
}

template <typename Scalar>
void DynamStateSet<Scalar>::merge(const DynamStateSet<Scalar> &dss)
{
  if (dss.vectorSize_ != vectorSize_)
  {
    fprintf(stderr, "Warning -- in DynamStateSet::merge : vectorSize mismatch -- function aborted\n");
    return;
  }

  stateSet_.reserve(stateSet_.size() + dss.stateSet_.size());
  typename vector<State>::iterator it;
  for (it = dss.stateSet_.begin(); it != dss.stateSet_.end(); ++it)
  {
    stateSet_.push_back(*it);
  }
}

template <typename Scalar>
void DynamStateSet<Scalar>::addState(const GenVector<DataType> & pos, const GenVector<DataType> & vel)
{
  stateSet_.push_back(DynamState<Scalar>(vectorSize_));
  stateSet_.back().setState(pos, vel);
}

// Get raw data to contiguous buffer
template <typename Scalar>
void DynamStateSet<Scalar>::getRaw(Scalar * buffer) const
{
  int shift = 2 * vectorSize_;
  int numStates = stateSet_.size();
  for (int i = 0; i < numStates; ++i)
  {
    stateSet_[i].getRaw(buffer + i * shift);
  }
}

// Orthogonalization
template <typename Scalar>
void DynamStateSet<Scalar>::addRawSetAndOG(DynamStateSet<Scalar> & dss, int numStatesToAdd, Scalar * rawDataPtr, const GenSparseMatrix<double> * K, const GenSparseMatrix<double> * M)
{

  numStatesToAdd = (numStatesToAdd > 0) ? numStatesToAdd : 0;
  int finalNumStates = stateSet_.size() + numStatesToAdd;
  stateSet_.reserve(finalNumStates);
  dss.stateSet_.reserve(finalNumStates);

  // Loop on states to add 
  for (int i = 0; i < numStatesToAdd; ++i)
  {
    State tempState(vectorSize_, rawDataPtr + 2 * vectorSize_ * i);
    addThisStateAndOG_(dss, tempState, K, M);  
  }

  // Check ortho
  /*fprintf(stderr, "Check ortho : %d / %d\n", stateSet_.size(), dss.stateSet_.size());
  for (int i = 0; i < stateSet_.size(); ++i)
  {
    for (int j = 0; j <  dss.stateSet_.size(); ++j)
    {
      fprintf(stderr, " %e", dss.stateSet_[i] * stateSet_[j]);  
    }
    fprintf(stderr, "\n");
  }*/
}

template <typename Scalar>
void DynamStateSet<Scalar>::addStateSetAndOG(DynamStateSet<Scalar> & dss, const DynamStateSet<Scalar> & toAdd, const GenSparseMatrix<double> * K, const GenSparseMatrix<double> * M)
{
  int finalNumStates = stateSet_.size() + dss.stateSet_.size();
  stateSet_.reserve(finalNumStates);
  dss.stateSet_.reserve(finalNumStates);

  int imax = toAdd.stateSet_.size();
  for (int i = 0; i < imax; ++i)
  {
    State tempState(toAdd.stateSet_[i]);
    addThisStateAndOG_(dss, tempState, K, M);
  }
}

template <typename Scalar>
void DynamStateSet<Scalar>::addStateAndOG(DynamStateSet<Scalar> & dss, const DynamState<Scalar> & state, const GenSparseMatrix<double> * K, const GenSparseMatrix<double> * M)
{
  int finalNumStates = stateSet_.size() + 1;
  stateSet_.reserve(finalNumStates);
  dss.stateSet_.reserve(finalNumStates);
  State tempState(state);
  addThisStateAndOG_(dss, tempState, K, M);
}

template <typename Scalar>
void DynamStateSet<Scalar>::addThisStateAndOG_(DynamStateSet<Scalar> & dss, DynamState<Scalar> & tempState, const GenSparseMatrix<double> * K, const GenSparseMatrix<double> * M)
{
  int jmax = stateSet_.size();
  for (int j = 0; j < jmax; ++j)
  {
    tempState.linAdd(- (dss.stateSet_[j] * tempState), stateSet_[j]);
  }

  State tempMKState(tempState, K, M);
  double norm_i = tempState * tempMKState;
  if (norm_i < 0.0) // Check for negative tangent stiffness matrix -- Should not happen in a perfect world, but happens anyway with bricks and beams (at least).
  {
    fprintf(stderr, "Warning -- in DynamStateSet<Scalar>::addThisStateAndOG_(...) : negative norm\n");
    return;  
  }

  norm_i = sqrt(norm_i);
  if (norm_i < abstol) return;
  
  norm_i = 1.0 / norm_i;
  tempState *= norm_i;
  tempMKState *= norm_i;
  stateSet_.push_back(tempState);
  dss.stateSet_.push_back(tempMKState);
}

template <typename Scalar>
void DynamStateSet<Scalar>::projectorOG(const DynamStateSet<Scalar> & finalSet, const DynamState<Scalar> & inputState, DynamState<Scalar> & outputState) const
{
  outputState = 0.0; 
 
  if (stateSet_.size() != finalSet.stateSet_.size()) fprintf(stderr, "Size mismatch init != final in DynamStateSet<Scalar>::projectorOG()\n");
 
  typename std::vector<State>::const_iterator init, final;
  for (init = stateSet_.begin(), final = finalSet.stateSet_.begin(); init != stateSet_.end() && final != finalSet.stateSet_.end(); ++init, ++final)
  {
    outputState.linAdd(inputState * (*init), *final);
  }
}

#endif


