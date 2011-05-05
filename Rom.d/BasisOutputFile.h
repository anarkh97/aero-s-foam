#ifndef ROM_BASISOUTPUTFILE_H
#define ROM_BASISOUTPUTFILE_H

#include <string>
#include <cstdio>

namespace Rom {

class BasisOutputFile {
public:
  const std::string &fileName() const { return fileName_; }

  int nodeCount() const  { return nodeCount_;  }
  int stateCount() const { return stateCount_; }

  enum StateCountStatus { UP_TO_DATE, OUTDATED };
  StateCountStatus stateCountStatus() const;
  void updateStateCountStatus();

  // NodeBufferType must implement double indexation ([i][j]) to yield doubles
  // First dimension is nodeCount
  // Second dimension is 6 (dofs per node)
  // Note: (double *)[6] is a valid type
  template <typename NodeBufferType>
  void stateAdd(const NodeBufferType &data);
  template <typename NodeBufferType>
  void stateAdd(const NodeBufferType &data, double headValue);
  
  BasisOutputFile(const std::string &fileName, int nodeCount);
 
  ~BasisOutputFile();

private:
  const std::string fileName_;
  const int nodeCount_;
  const int width_;
  const int precision_;

  int stateCount_;
  
  FILE *stream_;

  static const int STATE_COUNT_LENGTH;
  int stateCountOnFile_;
  void writeStateCount();
  void rewindAndWriteStateCount(); 

  void writeNodeCount();
  void writeStateHeader(double value);

  // Disallow copy & assignment
  BasisOutputFile(const BasisOutputFile&);
  BasisOutputFile& operator=(const BasisOutputFile&);
};

inline
BasisOutputFile::StateCountStatus
BasisOutputFile::stateCountStatus() const {
  return (stateCount() == stateCountOnFile_) ? UP_TO_DATE : OUTDATED;
}

template <typename NodeBufferType>
inline
void
BasisOutputFile::stateAdd(const NodeBufferType &data) {
  stateAdd(data, static_cast<double>(stateCount() + 1));
}

template <typename NodeBufferType>
void
BasisOutputFile::stateAdd(const NodeBufferType &data, double headValue) {
  writeStateHeader(headValue);
  
  for (int iNode = 0; iNode < nodeCount(); iNode++)  {
    std::fprintf(stream_, " % *.*E % *.*E % *.*E % *.*E % *.*E % *.*E\n",
                 width_, precision_, data[iNode][0], width_, precision_, data[iNode][1],
                 width_, precision_, data[iNode][2], width_, precision_, data[iNode][3],
                 width_, precision_, data[iNode][4], width_, precision_, data[iNode][5]);
  }

  stateCount_++;
}

} /* end namespace Rom */

#endif /* ROM_BASISOUTPUTFILE_H */
