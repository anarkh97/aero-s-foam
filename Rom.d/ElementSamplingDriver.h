#ifndef ROM_ELEMENTSAMPLINGDRIVER_H
#define ROM_ELEMENTSAMPLINGDRIVER_H

#include "DriverInterface.h"

#include "SparseNonNegativeLeastSquaresSolver.h"
#include "VecBasis.h"

class Domain;
class Corotator;
class GeomState; 
class StaticTimers;
template <typename Scalar> class GenFullSquareMatrix; 

namespace Rom {

class ElementSamplingDriver : public DriverInterface {
public:
  virtual void solve(); // overriden

  explicit ElementSamplingDriver(Domain *);
  ~ElementSamplingDriver();

private:
  void preProcess();
  void buildDomainCdsa();  

  template <typename DblFwdIt>
  void assembleTrainingData(const VecBasis &snapshots, DblFwdIt timeStampFirst, const VecBasis &podBasis, VecBasis &elemContributions);

  int elementCount() const;
  int vectorSize() const;

  Domain *domain_;

  Corotator **corotators_;
  GeomState *geomState_;
  GenFullSquareMatrix<double> *kelArray_;
  
  SparseNonNegativeLeastSquaresSolver solver_;
};

} // end namespace Rom

Rom::DriverInterface *elementSamplingDriverNew(Domain *);

#endif /* ROM_ELEMENTSAMPLINGDRIVER_H */
