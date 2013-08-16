#ifndef ROM_ELEMENTSAMPLINGDRIVER_H
#define ROM_ELEMENTSAMPLINGDRIVER_H

#include "DriverInterface.h"

#include "Problems.d/DynamDescr.h"

#include "SparseNonNegativeLeastSquaresSolver.h"
#include <Math.d/Vector.h>
#include "VecBasis.h"
#include <vector>
#include "MeshDesc.h"
#include "FileNameInfo.h"

class Domain;
class Corotator;
class GeomState; 
class StaticTimers;
class FileNameInfo;
class MeshDesc;

template <typename Scalar> class GenFullSquareMatrix; 
template <typename Scalar> class GenVector;
typedef GenVector<double> Vector;

namespace Rom {
void outputMeshFile(const FileNameInfo &fileInfo, const MeshDesc &mesh, const int podVectorCount);

template<typename MatrixBufferType = std::vector<double>, typename SizeType = size_t>
class ElementSamplingDriver : public SingleDomainDynamic, public DriverInterface {
public:
  virtual void solve(); // overriden
  void computeSolution(Vector &trainingTarget, Vector &solution, double relativeTolerance, bool verboseFlag = true);
  void postProcess(Vector &solution, bool verboseFlag = true);
  VecBasis& podBasis() { return podBasis_; }
  VecBasis& displac() { return displac_; }
  VecBasis* veloc() { if(!veloc_) veloc_ = new VecBasis; return veloc_; }
  VecBasis* accel() { if(!accel_) accel_ = new VecBasis; return accel_; }
  int vectorSize() const;
  void timeStampsIs(const std::vector<double> &tst) { timeStamps_ = tst; }

  explicit ElementSamplingDriver(Domain *);
  ~ElementSamplingDriver();

  virtual void preProcess();
  void assembleTrainingData(Vector &trainingTarget);
  void clean();

protected:
  int elementCount() const;

  void buildDomainCdsa();
  Domain *domain_;

  Corotator **corotators_;
  GeomState *geomState_;
  GenFullSquareMatrix<double> *kelArray_, *melArray_;

  VecBasis podBasis_;
  VecBasis displac_;
  VecBasis *veloc_;
  VecBasis *accel_;
  std::vector<double> timeStamps_;
  SparseNonNegativeLeastSquaresSolver<MatrixBufferType,SizeType> solver_;
};

} // end namespace Rom

Rom::DriverInterface *elementSamplingDriverNew(Domain *);

#endif /* ROM_ELEMENTSAMPLINGDRIVER_H */
