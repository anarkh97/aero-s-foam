#ifndef ROM_BASISOUTPUTFILE_H
#define ROM_BASISOUTPUTFILE_H

#include <string>
#include <cstdio>

class BasisOutputFile {
public:
  const std::string &fileName() const { return fileName_; }

  BasisOutputFile(const std::string &fileName, int nodeCount);
  ~BasisOutputFile();

  int stateCount() const { return stateCount_; }

  void stateAdd(const double(*data)[6]);

private:
  const std::string fileName_;
  const int nodeCount_;
  const int width_;
  const int precision_;

  FILE *stream_;
  int stateCount_;

  // Disallow copy & assignment
  BasisOutputFile(const BasisOutputFile&);
  BasisOutputFile& operator=(const BasisOutputFile&);
};

#endif /* ROM_BASISOUTPUTFILE_H */
