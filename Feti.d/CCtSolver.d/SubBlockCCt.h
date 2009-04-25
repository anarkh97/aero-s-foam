#ifndef _SUBBLOCKCCT_H_
#define _SUBBLOCKCCT_H_

#include <Feti.d/CCtSolver.d/CCtSolver.h>

template<class Scalar>
class SubBlockCCtSolver : public CCtSolver<Scalar>
{
  public:
    SubBlockCCtSolver(Connectivity *mpcToMpc, Connectivity *mpcToSub, 
                      int numSubsWithMpcs, GenSubDomain<Scalar> **subsWithMpcs,
                      FSCommunicator *fetiCom, Connectivity *cpuToSub);
    ~SubBlockCCtSolver();
    void reSolve(GenDistrVector<Scalar> &v);
    void zeroAll();
    void assemble();
    void factor();
  private:
    Connectivity *mpcToSub;
    FSCommPattern<Scalar> *mpcvPat;
    FSCommPattern<Scalar> *cctPat;
    Connectivity *cpuToSub;
    int myCPU;
    void solveLocalCCt(int iSub, GenDistrVector<Scalar> &v);
    void sendMpcInterfaceVec(int iSub, GenDistrVector<Scalar> &v);
    void combineMpcInterfaceVec(int iSub, GenDistrVector<Scalar> &v);
};
                                          
#ifdef _TEMPLATE_FIX_
  #include<Feti.d/CCtSolver.d/SubBlockCCt.C>
#endif
                                          
#endif

