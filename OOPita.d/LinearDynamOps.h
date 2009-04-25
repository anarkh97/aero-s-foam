#ifndef PITA_LINEARDYNAMOPS_H
#define PITA_LINEARDYNAMOPS_H

#include "DynamOps.h"
#include "Types.h"
#include <Driver.d/Dynam.h>
#include <map>

template<typename Scalar> class SingleDomainDynamic;;

namespace Pita {
 
class LinearDynamOps : public DynamOps {
public:
  typedef Fwk::Ptr<LinearDynamOps> Ptr;
  typedef Fwk::Ptr<const LinearDynamOps> PtrConst;

  typedef double Scalar;
  typedef GenDynamMat<Scalar> DynOpsType;
  typedef GenSparseMatrix<Scalar> MatrixType;
  typedef GenSolver<Scalar> SolverType;

  /* Not implemented // virtual size_t unconCount() const { return unconCount_; } */

  virtual const MatrixType * massMatrix()        const;
  virtual const MatrixType * stiffnessMatrix()   const;
  virtual const MatrixType * dampingMatrix()     const;
  virtual const SolverType * dynamicMassSolver() const;
  
  MatrixType * massMatrix();
  MatrixType * stiffnessMatrix();
  MatrixType * dampingMatrix();
  SolverType * dynamicMassSolver();

  const DynOpsType * dynamMat() const { return dynamMat_; }

  class Manager;
  friend class Manager;

protected:
  explicit LinearDynamOps(DynOpsType * dMat);

private:
  DynOpsType * dynamMat_;
  /* Not implemented // size_t unconCount_; */
};


class GeneralizedAlphaParameter {
public:
  double alpham() const { return alpham_; }
  double alphaf() const { return alphaf_; }
  double beta()   const { return beta_;   }
  double gamma()  const { return gamma_;  }

  // rhoInfinity = dissipation coefficient of forward time-integration when frequency -> infinity
  double rhoInfinity() const { return rhoInfinity_; }
  void rhoInfinityIs(double rhoInf); // 0.0 <= rhoInf <= 1.0
  
  Seconds timeStepSize() const { return timeStepSize_; }
  void timeStepSizeIs(Seconds dt);

  explicit GeneralizedAlphaParameter(Seconds stepSize, double rhoInf = 0.0);
  // Default GeneralizedAlphaParameter(const GeneralizedAlphaParameter &);
  // Default GeneralizedAlphaParameter & operator=(const GeneralizedAlphaParameter &);

  bool operator==(const GeneralizedAlphaParameter &) const;
  bool operator<(const GeneralizedAlphaParameter &) const; // (Arbitrary) total ordering defined for std::map

private:
  double alpham_, alphaf_, beta_, gamma_;
  double rhoInfinity_;
  Seconds timeStepSize_;
};


class LinearDynamOps::Manager : public PtrInterface<LinearDynamOps::Manager> {
public:
  typedef Fwk::Ptr<LinearDynamOps::Manager> Ptr;
  typedef Fwk::Ptr<const LinearDynamOps::Manager> PtrConst;

  typedef double Scalar;
  typedef SingleDomainDynamic<Scalar> ProblemDescriptor;

  LinearDynamOps::Ptr dynOpsNew(const GeneralizedAlphaParameter & param);
  LinearDynamOps::Ptr dynOps(const GeneralizedAlphaParameter & param) const;
  void dynOpsDel(const GeneralizedAlphaParameter & param);
  size_t dynOpsCount() const;

  const ProblemDescriptor * probDesc() const { return probDesc_; }
  ProblemDescriptor * probDesc() { return probDesc_; }

  static Manager::Ptr New(ProblemDescriptor * pbDesc) {
    return new Manager(pbDesc);
  }

protected:
  explicit Manager(ProblemDescriptor * pbDesc);

private:
  ProblemDescriptor * probDesc_;

  typedef std::map<GeneralizedAlphaParameter, LinearDynamOps::Ptr> DynOpsMap;
  DynOpsMap dynOpsMap_;

  // Optimization possible when no damping
  bool noPhysicalDamping_;
};

} // end namespace Pita

#endif /* PITA_LINEARDYNAMOPS_H */
