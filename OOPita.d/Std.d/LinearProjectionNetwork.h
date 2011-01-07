#ifndef PITA_STD_LINEARPROJECTIONNETWORK_H
#define PITA_STD_LINEARPROJECTIONNETWORK_H

#include "Fwk.h"
#include "Types.h"

#include "../DynamOps.h"

class Communicator;
#include "SliceMapping.h"

#include "../DynamStatePlainBasis.h"
#include "../RankDeficientSolver.h"
#include <Math.d/FullSquareMatrix.h>

#include "AffineBasisCollector.h"

#include "../NamedTask.h"

#include <list>

namespace Pita { namespace Std {

class LinearProjectionNetwork : public Fwk::PtrInterface<LinearProjectionNetwork> {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearProjectionNetwork);

  // Data collection
  AffineBasisCollector * collector() const { return collector_.ptr(); }

  // Projection operators
  size_t reducedBasisSize() const { return projectionBasis_->stateCount(); }
  const DynamOps * metric() const { return metric_.ptr(); }
  const DynamStateBasis * projectionBasis() const { return projectionBasis_.ptr(); }
  const DynamStateBasis * propagatedBasis() const { return propagatedBasis_.ptr(); }
  const RankDeficientSolver * normalMatrixSolver() const { return normalMatrixSolver_.ptr(); }
  const FullSquareMatrix * reprojectionMatrix() const { return &reprojectionMatrix_; }

  // Execution
  NamedTask::Ptr projectionTaskNew();

  enum Kind {
    INITIAL = 0,
    FINAL = 1
  };

  typedef GenId<Kind> IterStateId;
  
  class StateId {
  public:
    IterationRank iteration() const { return iteration_; }
    SliceRank slice() const { return inIterId_.rank(); }
    Kind type() const { return inIterId_.type(); }

    StateId(IterationRank iter, IterStateId inIterId) :
      iteration_(iter), inIterId_(inIterId)
    {}

    bool operator==(const StateId & other) const {
      return iteration() == other.iteration() && inIterId_ == other.inIterId_;
    }

    bool operator<(const StateId & other) const {
      return iteration() == other.iteration() ? inIterId_ < other.inIterId_ : iteration() < other.iteration();
    }

  private:
    IterationRank iteration_;
    IterStateId inIterId_;
  };
  
  void print_debug();

  static Ptr New(const SliceMapping * mapping, Communicator * timeComm, const DynamOps * metric, size_t vecSize, RankDeficientSolver * solver) {
    return new LinearProjectionNetwork(mapping, timeComm, metric, vecSize, solver);
  }

protected:
  LinearProjectionNetwork(const SliceMapping * mapping, Communicator * timeComm, const DynamOps * metric, size_t vecSize, RankDeficientSolver * solver);
 
  // Numbering for Allgather 
  class StateExchgNumbering;
  class MatrixExchgNumbering;

  // Execution
  class Task;
  friend class Task;
  void buildProjection();

private:
  // Problem characteristics 
  size_t vectorSize_;
  DynamOps::PtrConst metric_;

  // Network & Communication
  SliceMapping::PtrConst mapping_;
  Fwk::Ptr<MatrixExchgNumbering> numbering_;
  Communicator * timeCommunicator_;

  // Projection operators
  DynamStatePlainBasis::Ptr projectionBasis_; 
  DynamStatePlainBasis::Ptr propagatedBasis_;
  RankDeficientSolver::Ptr normalMatrixSolver_;
  FullSquareMatrix reprojectionMatrix_;

  // Rank-deficient operators
  DynamStatePlainBasis::Ptr originalProjectionBasis_; 
  DynamStatePlainBasis::Ptr originalPropagatedBasis_; 
  FullSquareMatrix normalMatrix_;
  FullSquareMatrix transmissionMatrix_;

  // Local data collection
  AffineBasisCollector::Ptr collector_;
  typedef std::map<StateId, DynamState> LocalStateMap;
  LocalStateMap localState_; // Accumulated index to local state

  DISALLOW_COPY_AND_ASSIGN(LinearProjectionNetwork);
};


inline
OStream & operator<<(OStream & out, LinearProjectionNetwork::Kind k) {
  char c = (k == LinearProjectionNetwork::INITIAL) ? 'i' : 'f';
  return out << c;
}

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_LINEARPROJECTIONNETWORK_H */
