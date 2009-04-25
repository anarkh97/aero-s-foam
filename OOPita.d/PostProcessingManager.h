#ifndef PITA_POSTPROCESSINGMANAGER_H
#define PITA_POSTPROCESSINGMANAGER_H

#include "Fwk.h"
#include "LinearPostProcessor.h"
#include "IntegratorPropagator.h"
#include "LinearGenAlphaIntegrator.h"

#include <map>

namespace Pita {

// TODO: Handle other types of postprocessors

class PostProcessingReactor : public LinearGenAlphaIntegrator::NotifieeConst {
public:
  EXPORT_PTRINTERFACE_TYPES(PostProcessingReactor);

  const LinearPostProcessor * target() { return target_; }

  PostProcessor::FileSetId outputFileSet() const { return outputFileSet_; }

  virtual void onInitialCondition();
  virtual void onCurrentCondition();

  virtual void notifierIs(const DynamTimeIntegrator * notifier);

  static Ptr New(const LinearGenAlphaIntegrator * notifier, 
                 LinearPostProcessor * target,
                 PostProcessor::FileSetId outputFileSet) {
    return new PostProcessingReactor(notifier, target, outputFileSet);
  }

protected:
  PostProcessingReactor(const LinearGenAlphaIntegrator * n, 
                        LinearPostProcessor * t, 
                        PostProcessor::FileSetId ofs);

  void performOutput() const;

private:
  const LinearGenAlphaIntegrator * notifier_; // Need a specialized DynamTimeIntegrator
  LinearPostProcessor * target_;
  PostProcessor::FileSetId outputFileSet_;
};


class PostProcessingManager : public Fwk::PtrInterface<PostProcessingManager> {
public:
  EXPORT_PTRINTERFACE_TYPES(PostProcessingManager);

  const LinearPostProcessor * postProcessor() const { return postProcessor_.ptr(); }

  PostProcessor::FileSetId outputFileSet(const IntegratorPropagator * observedPropagator) const;
  void outputFileSetIs(const IntegratorPropagator * observedPropagator, PostProcessor::FileSetId fileSet);

  static Ptr New(SDDynamPostProcessor * basePostProcessor,
                 int localFileCount,
                 const int * localFileId) // Pointer to array of size localFileCount
  {     
    return new PostProcessingManager(basePostProcessor, localFileCount, localFileId);
  }

protected:
  class PropagatorReactor : public IntegratorPropagator::Notifiee {
  public:
    EXPORT_PTRINTERFACE_TYPES(PropagatorReactor);

    virtual void onInitialState();
    virtual void onFinalState();

    PostProcessor::FileSetId fileSet() const { return fileSet_; }
    void fileSetIs(PostProcessor::FileSetId fs) { fileSet_ = fs; }

    PropagatorReactor(const IntegratorPropagator * notifier,
                      LinearPostProcessor * target,
                      PostProcessor::FileSetId fileSet);

  private:
    const IntegratorPropagator * notifier_; // Need a specialized DynamPropagator
    LinearPostProcessor * target_;          // Smart pointer to be kept by creator of PropagatorReactor
    PostProcessor::FileSetId fileSet_;

    PostProcessingReactor::Ptr integratorReactor_; // Collects the output data
  };

  PostProcessingManager(SDDynamPostProcessor * basePostProcessor,
                        int localFileCount,
                        const int * localFileId);

private:
  LinearPostProcessor::Ptr postProcessor_;

  typedef std::map<const IntegratorPropagator *, PropagatorReactor::Ptr> PropagatorReactorContainer;
  PropagatorReactorContainer propagatorReactor_;
  
  DISALLOW_COPY_AND_ASSIGN(PostProcessingManager);
};

} // end namespace Pita

#endif /* PITA_POSTPROCESSINGMANAGER_H */
