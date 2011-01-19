#ifndef ROM_BASISOUTPUTFILE_H
#define ROM_BASISOUTPUTFILE_H

#include <string>
#include <cstdio>

class BasisOutputFile {
public:
  const std::string &fileName() const { return fileName_; }

  int nodeCount() const  { return nodeCount_;  }
  int stateCount() const { return stateCount_; }

  enum StateCountStatus { UP_TO_DATE, OUTDATED };
  StateCountStatus stateCountStatus() const;
  void updateStateCountStatus();

  void stateAdd(const double(*data)[6]);
  void stateAdd(const double(*data)[6], double headValue);
  
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
  
  // Disallow copy & assignment
  BasisOutputFile(const BasisOutputFile&);
  BasisOutputFile& operator=(const BasisOutputFile&);
};

inline
BasisOutputFile::StateCountStatus
BasisOutputFile::stateCountStatus() const {
  return (stateCount_ == stateCountOnFile_) ? UP_TO_DATE : OUTDATED;
}

#endif /* ROM_BASISOUTPUTFILE_H */
