#include "ElementSamplingDriver.C"

namespace Rom {

template
void
outputFullWeights<Vector,std::vector<int> >(const Vector &weights, const std::vector<int> &elemIds, int j);

template
void
outputFullWeights<std::vector<double>,std::vector<int> >(const std::vector<double> &weights, const std::vector<int> &elemIds, int j);

template
int
ElementSamplingDriver<std::vector<double>,size_t>
::elementCount() const;

template
int
ElementSamplingDriver<std::vector<double>,size_t>
::vectorSize() const;

template
ElementSamplingDriver<std::vector<double>,size_t>
::ElementSamplingDriver(Domain *d);

template
ElementSamplingDriver<std::vector<double>,size_t>
::~ElementSamplingDriver();

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::assembleTrainingData(const VecBasis &podBasis, const int podVectorCount, const VecBasis &displac,
                       const VecBasis *veloc, const VecBasis *accel, int j);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::assembleTrainingData(const std::vector<StackVector> &podBasis, const int podVectorCount, const std::vector<StackVector> &displac,
                       const std::vector<StackVector> *veloc, const std::vector<StackVector> *accel, int j);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::solve();

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::computeSolution(Vector &solution, double tol, bool verboseFlag);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::postProcess(Vector &solution, bool verboseFlag);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::preProcess();

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::buildDomainCdsa();

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::clean();

} // end namespace Rom
