#ifndef ROM_BASISINPUTFILE_H
#define ROM_BASISINPUTFILE_H

#include <string>
#include <cstdio>

#include <cassert>

namespace Rom {

class BasisInputFile {
public:
  // Static parameters
  const std::string &fileName() const { return fileName_; }
  
  int stateCount() const { return stateCount_; }
  int nodeCount() const { return nodeCount_; }

  // Iteration and retrieval
  int currentStateIndex() const { return currentStateIndex_; }
  double currentStateHeaderValue() const { return currentStateHeaderValue_; }
  template <typename NodeDof6Type>
  const NodeDof6Type &currentStateBuffer(NodeDof6Type &target) const;
  void currentStateIndexInc(); // Must have called currentStateBuffer() at least once before calling currentStateIndexInc()

  bool validCurrentState() const { return currentStateIndex() < stateCount(); }

  // Ctor & dtor
  explicit BasisInputFile(const std::string &fileName);

  ~BasisInputFile();

private:
  const std::string fileName_;
  
  int nodeCount_;
  int stateCount_;
  
  FILE *stream_;
  
  int currentStateIndex_;
  double currentStateHeaderValue_;
  mutable bool currentStateRead_;
  std::fpos_t currentStatePosition_;

  void readCurrentStateHeader();
  void positionAtStateStart() const;

  // Disallow copy & assignment
  BasisInputFile(const BasisInputFile&);
  BasisInputFile& operator=(const BasisInputFile&);
};

template <typename NodeDof6Type>
const NodeDof6Type &
BasisInputFile::currentStateBuffer(NodeDof6Type &target) const {
  positionAtStateStart();

  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    const int info = std::fscanf(stream_, "%le %le %le %le %le %le",
                                 &target[iNode][0], &target[iNode][1], &target[iNode][2], 
                                 &target[iNode][3], &target[iNode][4], &target[iNode][5]);
    assert(info == 6); 
  }

  currentStateRead_ = true;

  return target;
}

} /* end namespace Rom */

#endif /* ROM_BASISINPUTFILE_H */
