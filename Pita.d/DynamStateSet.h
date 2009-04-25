#ifndef _DYNAMSTATESET_H_
#define _DYNAMSTATESET_H_

#include <Pita.d/DynamState.h>
#include <vector>

template <typename Scalar> class GenSparseMatrix;

template <typename Scalar>
class DynamStateSet
{

  public:
  
    typedef Scalar DataType; 
    typedef DynamState<DataType> State;
    typedef GenVector<DataType> VectorType;
 
  protected:

    // Data Fields
    int vectorSize_;
    vector<State> stateSet_; 

    // Memory management
    void init_(int vectorSize, int maxNumStates);

    // Tolerance parameters for orthogonalization process
    static const double abstol = 1.0e-10;

  public:

    // Constructors & destructor
    DynamStateSet(int vectorSize = 0, int maxNumStates = 0) { init_(vectorSize, maxNumStates); }
    DynamStateSet(const DynamStateSet<Scalar> & dss) : vectorSize_(dss.vectorSize_), stateSet_(dss.stateSet_) {}
    ~DynamStateSet() {}

    // Accessors
    int numStates() const { return stateSet_.size(); }
    int maxNumStates() const { return stateSet_.capacity(); }
    int vectorSize() const { return vectorSize_; }
    State & operator[](int i) { return stateSet_.at(i); }
    const State & operator[](int i) const { return stateSet_.at(i); }

    // Global operations
    DynamStateSet<Scalar> & operator=(const DynamStateSet<Scalar> & dss) { copy(dss); return *this;}
    void reset(int vectorSize, int maxNumStates) { stateSet_.clear(); init_(vectorSize, maxNumStates); }
    void copy(const DynamStateSet<Scalar> &dss) { vectorSize_ = dss.vectorSize_; stateSet_ = dss.stateSet_; }
    void merge(const DynamStateSet<Scalar> &dss);
    void resize(int newMaxNumStates) { stateSet_.reserve(newMaxNumStates); }
    void clear() { stateSet_.clear(); }

    // Add states
    void addState(const DynamState<Scalar> & newState) { stateSet_.push_back(newState); }
    void addState(const VectorType & pos, const VectorType & vel);
    void addState() { stateSet_.push_back(DynamState<Scalar>(vectorSize_)); }
    
    // Orthogonalization
    void addRawSetAndOG(DynamStateSet<Scalar> &dss, int numStates, Scalar * rawDataPtr, const GenSparseMatrix<double> * K, const GenSparseMatrix<double> * M); 
    void addStateAndOG(DynamStateSet<Scalar> &dss, const DynamState<Scalar> &state, const GenSparseMatrix<double> * K, const GenSparseMatrix<double> * M);
    void addStateSetAndOG(DynamStateSet<Scalar> &dss, const DynamStateSet<Scalar> &toAdd, const GenSparseMatrix<double> * K, const GenSparseMatrix<double> * M);
    void projectorOG(const DynamStateSet<Scalar> &finalSet, const DynamState<Scalar> & inputState, DynamState<Scalar> & outputState) const;

    // Extract to raw data array
    void getRaw(Scalar * buffer) const;

  protected:
 
    void addThisStateAndOG_(DynamStateSet<Scalar> &dss, DynamState<Scalar> &tempState, const GenSparseMatrix<double> * K, const GenSparseMatrix<double> * M);
 
};

#ifdef _TEMPLATE_FIX_
#include <Pita.d/DynamStateSet.C>
#endif

#endif


