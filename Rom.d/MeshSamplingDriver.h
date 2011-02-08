#ifndef ROM_MESHSAMPLINGDRIVER_H
#define ROM_MESHSAMPLINGDRIVER_H

#include <vector>

class Domain;

class MeshSamplingDriver {
public:
  explicit MeshSamplingDriver(Domain *);

  void solve();

private:
  void preProcess();
  std::vector<int> determineSampleNodes();
  void buildReducedMesh(const std::vector<int> &);

  Domain *domain_;

  // Disallow copy and assignment
  MeshSamplingDriver(const MeshSamplingDriver &);
  MeshSamplingDriver &operator=(const MeshSamplingDriver &);
};

#endif /* ROM_MESHSAMPLINGDRIVER_H */
