#ifndef PITA_POSTPROCESSOR_H
#define PITA_POSTPROCESSOR_H

#include "Fwk.h"

#include <map>

class GeoSource;

namespace Pita {

// Constructor: Duplicate the underlying outputfiles
// Destructor: Close the output files
// Invariants: Track the open files, check that only "legal" files are opened
// Treat directly with GeoSource (global geoSource instance to be passed as a constructor parameter)

class PostProcessor : public Fwk::PtrInterface<PostProcessor> {
public:
  EXPORT_PTRINTERFACE_TYPES(PostProcessor);

  enum FileStatus {
    NO_FILE = 0,
    CLOSED,
    OPEN
  };

private:
  class FileSet; 

public:
  class FileSetId : public Fwk::Ordinal<FileSet, int> {
  public:
    explicit FileSetId(int v = -1) :
      Fwk::Ordinal<FileSet, int>(v)
    {}
  }; 
  
  GeoSource * geoSource() const { return geoSource_; }
 
  FileStatus fileStatus(FileSetId fileSetId) const;
  int fileSetCount() const; 
  
  virtual void fileStatusIs(FileSetId fileSet, FileStatus status);

protected:
  PostProcessor(GeoSource * geoSource,
                int localFileCount,
                const int * localFileId); // Pointer to array of size localFileCount
  
  virtual ~PostProcessor();

private:
  GeoSource * geoSource_;
  std::map<FileSetId, FileStatus> fileStatus_;

  DISALLOW_COPY_AND_ASSIGN(PostProcessor);
};
  
} // end namespace Pita

#endif /* PITA_POSTPROCESSOR_H */
