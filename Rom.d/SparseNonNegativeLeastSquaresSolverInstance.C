#include "SparseNonNegativeLeastSquaresSolver.C"

namespace Rom {

template
SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t>
::SparseNonNegativeLeastSquaresSolver();

template
void
SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t>
::problemSizeIs(long eqnCount, long unkCount);

template
void
SparseNonNegativeLeastSquaresSolver<std::vector<double>,size_t>
::solve();

#ifdef USE_STXXL
template
SparseNonNegativeLeastSquaresSolver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::SparseNonNegativeLeastSquaresSolver();

template
void
SparseNonNegativeLeastSquaresSolver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::problemSizeIs(long eqnCount, long unkCount);

template
void
SparseNonNegativeLeastSquaresSolver<stxxl::VECTOR_GENERATOR<double,16,32,8388608,stxxl::RC,stxxl::random>::result,stxxl::uint64>
::solve();
#endif

} // end namespace Rom
