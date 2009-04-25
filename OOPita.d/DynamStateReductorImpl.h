#ifndef PITA_DYNAMSTATEREDUCTORIMPL_H
#define PITA_DYNAMSTATEREDUCTORIMPL_H

#include "DynamStateReductor.h"

#include "DynamStateBasis.h"
#include "PivotedCholeskySolver.h"

#include <map>

namespace Pita {

class DynamStateReductorImpl : public DynamStateReductor {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamStateReductorImpl);

  virtual void initialStateIs(const DynamState & is); // Overriden
  
  const DynamStateBasis * reductionBasis() const { return reductionBasis_.ptr(); }
  const PivotedCholeskySolver * solver() const { return solver_.ptr(); }

  class Manager;

protected:
  DynamStateReductorImpl(const DynamStateBasis * reductionBasis,
                         const PivotedCholeskySolver * solver);

  void setReductionBasis(const DynamStateBasis * mb) { reductionBasis_ = mb; }
  void setSolver(const PivotedCholeskySolver * s) { solver_ = s; }

  void reset();

private:
  DynamStateBasis::PtrConst reductionBasis_;
  PivotedCholeskySolver::PtrConst solver_;

  friend class Manager;
};

class DynamStateReductorImpl::Manager : public DynamStateReductor::Manager {
public:
  EXPORT_PTRINTERFACE_TYPES(Manager);
  
  // Overriden members
  virtual DynamStateReductorImpl * instance(const String & key) const;
  virtual size_t instanceCount() const;
  virtual DynamStateReductorImpl * instanceNew(const String & key);
  virtual void instanceDel(const String & key);

  // Added members
  const DynamStateBasis * defaultReductionBasis() const { return defaultReductionBasis_.ptr(); }
  const PivotedCholeskySolver * defaultSolver() const { return defaultSolver_.ptr(); }

  void defaultSolverIs(const PivotedCholeskySolver * s);
  void defaultReductionBasisIs(const DynamStateBasis * mb);

  static Ptr New(const DynamStateBasis * defaultReductionBasis, const PivotedCholeskySolver * defaultSolver) {
    return new Manager(defaultReductionBasis, defaultSolver);
  }

protected:
  Manager(const DynamStateBasis * defaultReductionBasis, const PivotedCholeskySolver * defaultSolver);

  void resetInstances();

private:
  typedef std::map<String, DynamStateReductorImpl::Ptr> InstanceContainer;
  InstanceContainer instance_;

  DynamStateBasis::PtrConst defaultReductionBasis_;
  PivotedCholeskySolver::PtrConst defaultSolver_;
  
  //const SymFullMatrix & normalMatrix() const { return *normalMatrix_; }
  //void defaultNormalMatrixIs(const SymFullMatrix & nm);
  //const SymFullMatrix * normalMatrix_; // Only the upper triangle is used
  //static const SymFullMatrix * defaultNormalMatrix();
};

} // end namespace Pita

#endif /* PITA_DYNAMSTATEREDUCTORIMPL_H */
