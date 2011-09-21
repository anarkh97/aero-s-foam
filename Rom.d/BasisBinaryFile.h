#ifndef ROM_BASISBINARYFILE_H
#define ROM_BASISBINARYFILE_H

#include <Utils.d/BinFileHandler.h>

#include <string>

namespace Rom {

class BasisBinaryOutputFile {
public:
  const std::string &fileName() const { return fileName_; }

  int nodeCount() const  { return nodeCount_;  }
  int stateCount() const { return stateCount_; }

  enum StateCountStatus { UP_TO_DATE, OUTDATED };
  StateCountStatus stateCountStatus() const;
  void updateStateCountStatus();

  // NodeBufferType must implement double indexation ([i][j]) to yield a type convertible to (double &)
  // First dimension is nodeCount
  // Second dimension is 6 (#dofs per node)
  // Note: (double *)[6] is a valid type
  template <typename NodeBufferType>
  void stateAdd(const NodeBufferType &data);
  template <typename NodeBufferType>
  void stateAdd(const NodeBufferType &data, double headValue);
  
  BasisBinaryOutputFile(const std::string &fileName, int nodeCount);

private:
  const std::string fileName_;
  const int nodeCount_;
  int stateCount_;
 
  BinFileHandler binHandler_;

  void writeStateHeader(double headValue);
  
  template <typename Dof6Type>
  void writeCurrentNode(const Dof6Type &buffer);

  void writeCurrentNode(const double *buffer); 
  
  void incrementStateCount();

  // Disallow copy & assignment
  BasisBinaryOutputFile(const BasisBinaryOutputFile&);
  BasisBinaryOutputFile& operator=(const BasisBinaryOutputFile&);
};

template <typename Dof6Type>
void
BasisBinaryOutputFile::writeCurrentNode(const Dof6Type &buffer) {
  for (int iDof = 0; iDof < 6; ++iDof) {
    binHandler_.write(&buffer[iDof], 1);
  }
}

inline
void
BasisBinaryOutputFile::writeCurrentNode(const double *buffer) {
  binHandler_.write(buffer, 6);
}

template <typename NodeBufferType>
void
BasisBinaryOutputFile::stateAdd(const NodeBufferType &data, double headValue) {
  writeStateHeader(headValue);

  // Write values
  for (int iNode = 0; iNode < nodeCount(); ++iNode) {
    writeCurrentNode(data[iNode]);
  }

  incrementStateCount();
}

template <typename NodeBufferType>
void
BasisBinaryOutputFile::stateAdd(const NodeBufferType &data) {
  stateAdd(data, static_cast<double>(stateCount() + 1));
}


class BasisBinaryInputFile {
public:
  // Static parameters
  const std::string &fileName() const { return fileName_; }
  
  int stateCount() const { return stateCount_; }
  int nodeCount() const { return nodeCount_; }

  // Iteration and retrieval
  int currentStateIndex() const { return currentStateIndex_; }
  double currentStateHeaderValue() const { return currentStateHeaderValue_; }
  template <typename NodeDof6Type>
  const NodeDof6Type &currentStateBuffer(NodeDof6Type &target);
  // Must have called currentStateBuffer() at least once before calling currentStateIndexInc()
  // NOTE: This restriction could be lifted if necessary
  void currentStateIndexInc();

  bool validCurrentState() const { return currentStateIndex() < stateCount(); }

  // Ctor & dtor
  explicit BasisBinaryInputFile(const std::string &fileName);

protected:
  template <typename Dof6Type>
  void readCurrentNode(Dof6Type &buffer);

  void readCurrentNode(double *buffer); 
  
  void seekNextState();
  void seekCurrentState();
  void seekNode(int iNode);
  void validateNextOffset();

private:
  const std::string fileName_;
  
  int nodeCount_;
  int stateCount_;
  
  int currentStateIndex_;
  double currentStateHeaderValue_;

  BinFileHandler binHandler_;
  BinFileHandler::OffType savedOffset_;
  mutable bool validNextOffset_;

  // Disallow copy & assignment
  BasisBinaryInputFile(const BasisBinaryInputFile&);
  BasisBinaryInputFile& operator=(const BasisBinaryInputFile&);
};

template <typename Dof6Type>
void
BasisBinaryInputFile::readCurrentNode(Dof6Type &target) {
  for (int iDof = 0; iDof < 6; ++iDof) {
    binHandler_.read(&target[iDof], 1);
  }
}

template <typename NodeDof6Type>
const NodeDof6Type &
BasisBinaryInputFile::currentStateBuffer(NodeDof6Type &target) {
  if (validCurrentState()) {
    seekCurrentState();

    for (int iNode = 0; iNode < nodeCount(); ++iNode) {
      readCurrentNode(target[iNode]);
    }

    validateNextOffset();
  }

  return target;
}

} /* end namespace Rom */

#endif /* ROM_BASISBINARYFILE_H */
