#ifndef ROM_ELEMENTSAMPLINGDRIVER_H
#define ROM_ELEMENTSAMPLINGDRIVER_H

#include "DriverInterface.h"

#include "SparseNonNegativeLeastSquaresSolver.h"
#include <Math.d/Vector.h>
#include "VecBasis.h"

#include <vector>

class Domain;
class Corotator;
class GeomState; 
class StaticTimers;
template <typename Scalar> class GenFullSquareMatrix; 
template <typename Scalar> class GenVector;
typedef GenVector<double> Vector;

namespace Rom {

class ElementSamplingDriver : public DriverInterface {
public:
  virtual void solve(); // overriden
  void computeSolution(Vector &solution, bool verboseFlag = true);
  void postProcess(Vector &solution, bool firstTime = true, bool verboseFlag = true);
  VecBasis& podBasis() { return podBasis_; }
  VecBasis& displac() { return displac_; }
  int vectorSize() const;
  void timeStampsIs(const std::vector<double> &tst) { timeStamps_ = tst; }

  explicit ElementSamplingDriver(Domain *);
  ~ElementSamplingDriver();

protected:
  virtual void preProcess();
  void buildDomainCdsa();  
  template <typename DblFwdIt>
  void assembleTrainingData(const VecBasis &displac, DblFwdIt timeStampFirst, const VecBasis &podBasis,
                            typename SparseNonNegativeLeastSquaresSolver::MatrixBufferType::iterator elemContributions,
                            Vector &trainingTarget, VecBasis *veloc, VecBasis *accel);
  int elementCount() const;

  Domain *domain_;

  Corotator **corotators_;
  GeomState *geomState_;
  GenFullSquareMatrix<double> *kelArray_, *melArray_;

  VecBasis podBasis_;
  VecBasis displac_;
  VecBasis *veloc_;
  VecBasis *accel_;
  std::vector<double> timeStamps_;
  
  SparseNonNegativeLeastSquaresSolver solver_;
};

} // end namespace Rom

Rom::DriverInterface *elementSamplingDriverNew(Domain *);

#endif /* ROM_ELEMENTSAMPLINGDRIVER_H */
