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
::assembleTrainingData(Vector &trainingTarget);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::solve();

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::computeSolution(Vector &trainingTarget, Vector &solution, double tol, bool verboseFlag);

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
::assembleTrainingData(Vector &trainingTarget);

template
void
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::solve();

template
void
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::computeSolution(Vector &trainingTarget, Vector &solution, double tol, bool verboseFlag);

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
