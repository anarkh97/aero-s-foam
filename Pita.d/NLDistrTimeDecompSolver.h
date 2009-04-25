#ifndef _NLDISTRTIMEDECOMPSOLVER_H_
#define _NLDISTRTIMEDECOMPSOLVER_H_

#include <Math.d/Vector.h>
#include <Pita.d/NLTimeSlice.h>
#include <Pita.d/TimeSliceMapping.h>
#include <Pita.d/SimpleBuffer.h>
#include <Pita.d/PitaNonLinDynam.h>
#include <Pita.d/NLDynamTimeIntegrator.h>
#include <Pita.d/LinearizedTimeIntegrator.h>
#include <vector>

template <typename Scalar> class DynamState;
template <typename Scalar> class DynamStateSet;
class Communicator;
class Connectivity;

class NLDistrTimeDecompSolver
{
public:
  typedef GenVector<double> VecType;
  typedef VecType::DataType DataType;
  typedef DynamState<DataType> State;
  typedef DynamStateSet<DataType> StateSet;
    
  explicit NLDistrTimeDecompSolver(PitaNonLinDynamic * pbDesc);
  ~NLDistrTimeDecompSolver();
  void solve();

private:
  typedef std::vector<NLTimeSlice>::iterator sliceIterator;

  int myCPU;
  PitaNonLinDynamic * probDesc;
  NLDynamTimeIntegrator seqIntegrator; 
  LinearizedTimeIntegrator linIntegrator;
  Communicator * timeCom;
  TimeSliceMapping * sliceMap;
  std::vector<NLTimeSlice> timeSliceSet;
  sliceIterator firstActive, firstInactive; 
  SimpleBuffer<DataType> baseBuffer, seedBuffer;
 
  // Initialization
  void initialize();
  void buildConnectivities();
  void buildTimeSlices();

  // Algorithm steps
  void getInitialSeeds();
  void fineGridComputation();
  void improveBases();
  void clearBases();
  void computeCorrection();

  // Subroutines
  void improveBasesWithAllSeeds();
  void improveBasesWithLocalIncrements();

  // Timeslice routines 
  void coarseIntegrator(NLTimeSlice &);
  void fineIntegrator(NLTimeSlice &);
  void initializeSeqIntegrator(NLTimeSlice &, double);
  void computeSliceUpdate(NLTimeSlice &);
  void trivialSliceUpdate(NLTimeSlice &);

  // Base Management
  void reBuildLocalK(const NLTimeSlice &);
  void addStateToBase(NLTimeSlice &, State &);
  void addStateSetToBase(NLTimeSlice &, StateSet &);
  void addRawDataToBase(NLTimeSlice &, DataType *, int);

  // Helper functions
  void getIntegratorState(State &);
  void setIntegratorState(const State &);

  // Get number of unconstrained dofs (problem size)
  int getProbSize() const { return probDesc->solVecInfo(); }
};

inline void NLDistrTimeDecompSolver::initializeSeqIntegrator(NLTimeSlice & timeSlice, double timeStep)
{
  setIntegratorState(timeSlice.seedState);
  seqIntegrator.currentTimeIs(timeSlice.initialTime);
  seqIntegrator.timeStepIs(timeStep);
  seqIntegrator.currentTimeStepNumberIs(0);
}

inline void NLDistrTimeDecompSolver::computeSliceUpdate(NLTimeSlice & timeSlice)
{
  timeSlice.jumpState = timeSlice.seedState;
  timeSlice.jumpState -= timeSlice.oldSeedState;
  
  // Convergence criterion to be figured out
  timeSlice.converged = false;
  //fprintf(stderr, "TS # %d : %d vectors in base\n", this->sliceRank, this->seedBase.numStates());
  timeSlice.orthoBase.projectorOG(timeSlice.propBase, timeSlice.jumpState, timeSlice.nextSeedState);
  //fprintf(stderr, "TS #%d, norm U(j) = %e, norm C(j+1) = %e\n", this->getRank(), localNorm(this->jumpState), localNorm(this->nextSeedState));
  timeSlice.nextSeedState += timeSlice.propState;
}

inline void NLDistrTimeDecompSolver::trivialSliceUpdate(NLTimeSlice & timeSlice)
{
  timeSlice.converged = true;
  timeSlice.jumpState = 0.0;
  timeSlice.nextSeedState = timeSlice.propState;
}

inline void NLDistrTimeDecompSolver::addStateToBase(NLTimeSlice & timeSlice, NLDistrTimeDecompSolver::State & state)
{
  timeSlice.seedBase.addStateAndOG(timeSlice.orthoBase, state, probDesc->getStiffMatrix(), probDesc->getMassMatrix());
}

inline void NLDistrTimeDecompSolver::addRawDataToBase(NLTimeSlice & timeSlice, NLDistrTimeDecompSolver::DataType * dataPtr, int numStates)
{
  timeSlice.seedBase.addRawSetAndOG(timeSlice.orthoBase, numStates, dataPtr, probDesc->getStiffMatrix(), probDesc->getMassMatrix());  
}

inline void NLDistrTimeDecompSolver::addStateSetToBase(NLTimeSlice & timeSlice, NLDistrTimeDecompSolver::StateSet & stateSet)
{
  timeSlice.seedBase.addStateSetAndOG(timeSlice.orthoBase, stateSet, probDesc->getStiffMatrix(), probDesc->getMassMatrix());
}

inline void NLDistrTimeDecompSolver::getIntegratorState(NLDistrTimeDecompSolver::State & state)
{
  seqIntegrator.getCurrentVelocity(state.vel());
  seqIntegrator.getCurrentDisplacement(state.disp());
}

inline void NLDistrTimeDecompSolver::setIntegratorState(const NLDistrTimeDecompSolver::State & state)
{
  seqIntegrator.setCurrentVelocity(state.vel());
  seqIntegrator.setCurrentDisplacement(state.disp());
}

#endif
