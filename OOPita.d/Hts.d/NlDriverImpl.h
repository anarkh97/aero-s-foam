#ifndef PITA_HTS_NLDRIVERIMPL_H
#define PITA_HTS_NLDRIVERIMPL_H

#include "Fwk.h"
#include "Types.h"

#include "../NlDriver.h"

#include "../PitaNonLinDynam.h"
#include <Utils.d/SolverInfo.h>
#include <Driver.d/Domain.h>
#include <Comm.d/Communicator.h>

#include "SliceMapping.h"

namespace Pita { namespace Hts {

class NlDriverImpl : public NlDriver {
public:
  EXPORT_PTRINTERFACE_TYPES(NlDriverImpl);

  virtual void solve();

  GeoSource * geoSource() const { return geoSource_; }
  Domain * domain() const { return domain_; }
  SolverInfo * solverInfo() const { return solverInfo_; }
  Communicator * baseComm() const { return baseComm_; }
  
  static Ptr New(PitaNonLinDynamic * pbDesc,
                 GeoSource * geoSource,
                 Domain * domain, 
                 SolverInfo * solverInfo,
                 Communicator * baseComm) {
    return new NlDriverImpl(pbDesc, geoSource, domain, solverInfo, baseComm);
  }

protected:
  explicit NlDriverImpl(PitaNonLinDynamic *, GeoSource *, Domain *, SolverInfo *, Communicator *);

  void preprocess();

private:
  /* Primary sources */
  GeoSource * geoSource_;
  Domain * domain_;
  SolverInfo * solverInfo_;
  Communicator * baseComm_;
  
  /* Space-domain */
  size_t vectorSize_;
  
  /* Time-domain */ 
  Seconds fineTimeStep_;
  TimeStepCount halfSliceRatio_;
  TimeStepCount sliceRatio_;
  Seconds coarseTimeStep_;
  Seconds initialTime_;
  Seconds finalTime_;

  /* Load balancing */ 
  SliceMapping::Ptr mapping_;

  /* Other parameters */
  IterationRank lastIteration_;
  double projectorTolerance_;
};

} /* end namespace Hts */ } /* end namespac Pta */

Pita::NlDriver::Ptr nlPitaDriverNew(Pita::PitaNonLinDynamic * problemDescriptor);

#endif /* PITA_HTS_NLDRIVERIMPL_H */
