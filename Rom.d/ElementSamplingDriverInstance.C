#include "ElementSamplingDriver.C"

namespace Rom {

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
                       const VecBasis *veloc, const VecBasis *accel);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::assembleTrainingData(const std::vector<StackVector> &podBasis, const int podVectorCount, const std::vector<StackVector> &displac,
                       const std::vector<StackVector> *veloc, const std::vector<StackVector> *accel);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::assembleWeightedTrainingData(const std::vector<StackVector> &podBasis, const int podVectorCount, const std::vector<StackVector> &displac,
                               const std::vector<StackVector> *veloc, const std::vector<StackVector> *accel, const Vector &weights);

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

#ifdef USE_STXXL
template
int
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::elementCount() const;

template
int
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::vectorSize() const;

template
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::ElementSamplingDriver(Domain *d);

template
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::~ElementSamplingDriver();

template
void
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::assembleTrainingData(const VecBasis &podBasis, const int podVectorCount, const VecBasis &displac,
                       const VecBasis *veloc, const VecBasis *accel);

template
void
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::solve();

template
void
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::computeSolution(Vector &solution, double tol, bool verboseFlag);

template
void
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::postProcess(Vector &solution, bool verboseFlag);

template
void
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::preProcess();

template
void
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::buildDomainCdsa();

template
void
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::clean();
#endif

} // end namespace Rom
