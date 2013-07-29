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
::assembleTrainingData(const VecBasis &displac, std::vector<double>::iterator timeStampFirst, const VecBasis &podBasis,
                       std::vector<double>::iterator elemContributions, Vector &trainingTarget, VecBasis *veloc, VecBasis *accel);

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::solve();

template
void
ElementSamplingDriver<std::vector<double>,size_t>
::computeSolution(Vector &solution, bool verboseFlag);

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
::assembleTrainingData(const VecBasis &displac, std::vector<double>::iterator timeStampFirst, const VecBasis &podBasis,
                       stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result::iterator elemContributions,
                       Vector &trainingTarget, VecBasis *veloc, VecBasis *accel);

template
void
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::solve();

template
void
ElementSamplingDriver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::computeSolution(Vector &solution, bool verboseFlag);

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
#endif

} // end namespace Rom
