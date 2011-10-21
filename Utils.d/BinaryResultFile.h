#ifndef UTILS_BINARYRESULTFILE_H
#define UTILS_BINARYRESULTFILE_H

#include "BinFileHandler.h"

#include <string>
#include <vector>

class BinaryResultOutputFile {
public:
  const std::string &pathName() const { return pathName_; }

  double version() const { return binHandler_.getVersion(); }
  
  // Prelude
  const std::string &description() const { return description_; }
  int dataType() const { return dataType_; }
  int itemCount() const { return itemCount_; }
  int itemDimension() const { return itemDimension_; }
  int stateCount() const { return stateCount_; }
  
  typedef std::vector<int>::const_iterator ItemIdIterator;
  ItemIdIterator itemIdBegin() const { return itemIds_.begin(); }
  ItemIdIterator itemIdEnd() const { return itemIds_.end(); }

  // Data
  int localItemCount() const { return itemIds_.size(); }
  int localDataSize() const { return localItemCount() * itemDimension(); }
  // Writes (localDataSize) x scalars
  void stateAdd(double stamp, const double *data);

  // Constructors
  BinaryResultOutputFile(const std::string &pathName, int dataType, const std::string &description,
                         int itemCount, int itemDimension,
                         double version);
  
  template <typename IdxInpIt>
  BinaryResultOutputFile(const std::string &pathName, int dataType, const std::string &description,
                         int itemCount, int itemDimension,
                         int localOffset, IdxInpIt localIdBegin, IdxInpIt localIdEnd,
                         double version);

private:
  int stateSize() const { return itemCount() * itemDimension(); }
  bool isMaster() const { return localOffset_ == 0; }

  void writePrelude();
  void updateStateCount();

  std::vector<int>::size_type descriptionSize() const { return description_.size(); }

  std::string pathName_;

  int dataType_;
  std::string description_;
  int itemCount_;
  int itemDimension_;
  int stateCount_;
  
  std::vector<int> itemIds_;
  int localOffset_;

  BinFileHandler binHandler_;

  // Disallow copy and assignment
  BinaryResultOutputFile(const BinaryResultOutputFile &);
  BinaryResultOutputFile &operator=(const BinaryResultOutputFile &);
};

template <typename IdxInpIt>
BinaryResultOutputFile::BinaryResultOutputFile(const std::string &pathName, int dataType, const std::string &description,
                                               int itemCount, int itemDimension,
                                               int localOffset, IdxInpIt localIdBegin, IdxInpIt localIdEnd,
                                               double version) :
  pathName_(pathName),
  dataType_(dataType),
  description_(description),
  itemCount_(itemCount),
  itemDimension_(itemDimension),
  stateCount_(0),
  itemIds_(localIdBegin, localIdEnd),
  localOffset_(localOffset),
  binHandler_(pathName.c_str(), isMaster() ? "ws" : "ws+", version)
{
  writePrelude();
}


class BinaryResultInputFile {
public:
  const std::string &pathName() const { return pathName_; }
 
  double version() const { return binHandler_.getVersion(); }

  // Prelude
  const std::string &description() const { return description_; }
  int dataType() const { return dataType_; }
  int itemCount() const { return itemIds_.size(); }
  int itemDimension() const { return itemDimension_; }
  int stateCount() const { return stateCount_; }
  
  typedef std::vector<int>::const_iterator ItemIdIterator;
  ItemIdIterator itemIdBegin() const { return itemIds_.begin(); }
  ItemIdIterator itemIdEnd() const { return itemIds_.end(); } 

  // Data
  int stateRank() const { return stateRank_; }
  double stateStamp();
  const double* state(double *buffer);

  void stateRankInc();

  // Constructor
  explicit BinaryResultInputFile(const std::string &pathName);

private:
  int stateSize() const { return itemCount() * itemDimension(); }
  
  std::vector<int>::size_type descriptionSize() const { return description_.size(); }
  
  std::string pathName_;

  int dataType_;
  std::string description_;
  int itemDimension_;
  int stateCount_;
  
  std::vector<int> itemIds_;

  BinFileHandler binHandler_;

  int stateRank_;

  // Disallow copy and assignment
  BinaryResultInputFile(const BinaryResultInputFile &);
  BinaryResultInputFile &operator=(const BinaryResultInputFile &);
};

#endif /*UTILS_BINARYRESULTFILE_H */
