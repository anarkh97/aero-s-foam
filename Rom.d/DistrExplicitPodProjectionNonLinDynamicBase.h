#ifndef ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMICBASE_H
#define ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMICBASE_H

#include <Paral.d/MDDynam.h>
#include <Driver.d/GeoSource.h>
#include <cstdio>
#include "DistrVecBasis.h"

namespace Rom {

class DistrExplicitPodPostProcessor : public MultiDomDynPostProcessor {

public:
 
  DistrExplicitPodPostProcessor(DecDomain *, StaticTimers* , DistrGeomState *, Corotator ***);
  ~DistrExplicitPodPostProcessor();

  void dynamOutput(int, double, MDDynamMat &, DistrVector &, DistrVector *aeroF, SysState<DistrVector> &);
  void printPODSize(int);

  
private:

  OutputInfo *oinfo;
  DistrGeomState *geomState;
  DecDomain *decDomain;
  StaticTimers *times;
  int podSize;
  FILE *fn;

  int numOutInfo;    
};

class DistrExplicitPodProjectionNonLinDynamicBase : public MultiDomainDynam {
public:
  explicit DistrExplicitPodProjectionNonLinDynamicBase(Domain *);

  // Overriding via hiding
  void preProcess(); // Additional pre-processing
  const DistrInfo &solVecInfo();
  DistrInfo &reducedVecInfo();
  void getInitState(SysState<DistrVector> &);
  void updateState(double, DistrVector&, DistrVector&);
  void computeExtForce2(SysState<DistrVector> &distState,
                        DistrVector &f, DistrVector &cnst_f, int tIndex,
                        double t, DistrVector *aero_f=0,
                        double gamma=0.5, double alphaf=0.0);
  void getConstForce(DistrVector& v);
  void getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex);
  void printFullNorm(DistrVector &);
  DistrExplicitPodPostProcessor *getPostProcessor();
  MDDynamMat * buildOps(double, double, double); // Reduced-order matrix solver

protected:
  Domain        *domain;
  DistrInfo      reducedInfo;
  DistrVecBasis  projectionBasis_;
  GenVecBasis<double, GenDistrVector> normalizedBasis_;
  //dummy Variables to fascillitate computation on Reduced Coordinates
  StaticTimers * times;
  DistrVector  * fExt;
  DistrVector  * fInt;
  DistrVector  * cnst_fBig;
  DistrVector  * aero_fBig;
  DistrVector  * d_n;
  DistrVector  * v_n;
  DistrVector  * a_n;
  DistrVector  * v_p;
  DistrVector  * tempVec;
  SysState<DistrVector> *dummyState;
  bool haveRot;
  GenParallelSolver<double> * fullMassSolver;

private:
  DistrExplicitPodPostProcessor *mddPostPro;
  // Disallow copy and assignment
  DistrExplicitPodProjectionNonLinDynamicBase(const DistrExplicitPodProjectionNonLinDynamicBase &);
  DistrExplicitPodProjectionNonLinDynamicBase &operator=(const DistrExplicitPodProjectionNonLinDynamicBase &);

};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMICBASE_H */
