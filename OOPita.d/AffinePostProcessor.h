#ifndef PITA_AFFINEPOSTPROCESSOR_H
#define PITA_AFFINEPOSTPROCESSOR_H

#include "Fwk.h"
#include "Types.h"
#include "PostProcessor.h"

#include "DynamState.h"
#include <map>

class SDDynamPostProcessor;

// Mostly a hack, to output the updated state, not only the update

namespace Pita {

class LinearGenAlphaIntegrator;

class AffinePostProcessor : public GenPostProcessor<LinearGenAlphaIntegrator> {
public:
  EXPORT_PTRINTERFACE_TYPES(AffinePostProcessor);

  virtual void outputNew(FileSetId fileSetId, const LinearGenAlphaIntegrator * oi);

  SDDynamPostProcessor * basePostProcessor() const { return basePostProcessor_; }

  static Ptr New(GeoSource * geoSource,
                 int localFileCount,
                 const int * localFileId,
                 SDDynamPostProcessor * basePostProcessor) {
    return new AffinePostProcessor(geoSource, localFileCount, localFileId, basePostProcessor);
  }

protected:
  AffinePostProcessor(GeoSource * gs, int lfc, const int * lfi, SDDynamPostProcessor * bpp);

private:
  SDDynamPostProcessor * basePostProcessor_;


  typedef std::map<TimeStepCount, DynamState> ConstantTermMap;
  typedef std::map<bool, ConstantTermMap> DirectionMap; // true = forward, false = backward
  typedef std::map<FileSetId, DirectionMap> FileSetMap;

  FileSetMap constantTermMap_;
};

} // namespace Pita

#endif /* PITA_AFFINEPOSTPROCESSOR_H */