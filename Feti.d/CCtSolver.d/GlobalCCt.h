#ifndef _GLOBALCCT_H_
#define _GLOBALCCT_H_

#include <Feti.d/CCtSolver.d/CCtSolver.h>

template<class Scalar>
class GlobalCCtSolver : public CCtSolver<Scalar>
{
  public:
    GlobalCCtSolver(Connectivity *_mpcToMpc, Connectivity *_mpcToCpu, int _numSubsWithMpcs,
                    GenSubDomain<Scalar> **_subsWithMpcs, FetiInfo *finfo, FSCommunicator *_fetiCom);
    ~GlobalCCtSolver();
    void reSolve(GenDistrVector<Scalar> &v);
    void zeroAll();
    void assemble();
    void factor();
  private:
    SimpleNumberer *mpcEqNums;
    GenSolver<Scalar> *CCtsolver;
    void computeSubContributionToGlobalCCt(int i, SimpleNumberer *mpcEqNums);
    void extractMpcResidual(int iSub, GenDistrVector<Scalar> &v, GenVector<Scalar> &mpcv1);
    void insertMpcResidual(int iSub, GenDistrVector<Scalar> &v, GenVector<Scalar> &mpcv1);
};

#ifdef _TEMPLATE_FIX_
  #include<Feti.d/CCtSolver.d/GlobalCCt.C>
#endif
                                                                                                                                              
#endif

