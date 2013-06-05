#ifndef ROM_SUBELEMENTSAMPLINGDRIVER_H
#define ROM_SUBELEMENTSAMPLINGDRIVER_H

#include "ElementSamplingDriver.h"

namespace Rom {

class SubElementSamplingDriver : public ElementSamplingDriver {
public:
  explicit SubElementSamplingDriver(Domain *);
  //~SubElementSamplingDriver();

private:
  void preProcess();
/*
  void buildDomainCdsa();  

  template <typename DblFwdIt>
  void assembleTrainingData(const VecBasis &snapshots, DblFwdIt timeStampFirst, const VecBasis &podBasis,
                            typename SparseNonNegativeLeastSquaresSolver::MatrixBufferType::iterator elemContributions, 
                            Vector &trainingTarget);

  int elementCount() const;
  int vectorSize() const;

  Domain *domain_;

  Corotator **corotators_;
  GeomState *geomState_;
  GenFullSquareMatrix<double> *kelArray_;
  
  SparseNonNegativeLeastSquaresSolver solver_;
*/
};

} // end namespace Rom

Rom::DriverInterface *subElementSamplingDriverNew(Domain *);

#endif /* ROM_SUBELEMENTSAMPLINGDRIVER_H */
