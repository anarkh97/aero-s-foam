#ifndef ROM_DISTRSNAPSHOTCLUSTERINGSOLVER_H
#define ROM_DISTRSNAPSHOTCLUSTERINGSOLVER_H

#if defined(USE_SCALAPACK) && defined(USE_EIGEN3)
class Communicator;
class SCDoubleMatrix;

#include <Eigen/Core>
#include <vector>

namespace Rom {

class DistrSnapshotClusteringSolver {
public:
  // Local data distribution
  int localRows() const { return localRows_; }

  // Local buffers: Internal column-major ordering, zero-based indexing
  // Local matrix buffer: [localRows by colCount]
  double *matrixColBuffer(int col);
  // Local cluster buffer: [localRows by clusterDimension]
  const double *clusterColBuffer(int i, int col) const;
  const int clusterColCount(int i) const;

  void solve();

  DistrSnapshotClusteringSolver(Communicator * comm, int rowCount, int colCount, int localRows, int numClusters,
                                int blockSize);

private:
  // Disallow copy & assignment
  DistrSnapshotClusteringSolver(const DistrSnapshotClusteringSolver &);
  DistrSnapshotClusteringSolver & operator=(const DistrSnapshotClusteringSolver &);

  Communicator * communicator_;
  int rowCount_, colCount_, localRows_, numClusters_, blockSize_;

  Eigen::MatrixXd matrixBuffer_;
  std::vector<std::vector<int> > clusterCols_;
};

/* Helper functions for buffer access */
inline
double *
DistrSnapshotClusteringSolver::matrixColBuffer(int col) {
  return matrixBuffer_.data() + col*localRows_;
}

inline
const double *
DistrSnapshotClusteringSolver::clusterColBuffer(int i, int col) const {
  return matrixBuffer_.data() + clusterCols_[i][col]*localRows_;
}

inline
const int
DistrSnapshotClusteringSolver::clusterColCount(int i) const {
  return clusterCols_[i].size();
}

} // end namespace Rom
#endif

#endif /* ROM_DISTRSNAPSHOTCLUSTERINGSOLVER_H */
