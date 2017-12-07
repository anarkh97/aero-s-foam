#ifndef _CCTSOLVER_H_
#define _CCTSOLVER_H_

template <typename Scalar>
class GenDistrVector;
template <typename Scalar>
class GenSubDomain;

class FSCommunicator;
class Connectivity;

template<class Scalar>
class CCtSolver
{
  public:
    virtual void reSolve(GenDistrVector<Scalar> &v) = 0;
    virtual void zeroAll() = 0;
    virtual void assemble() = 0;
    virtual void factor() = 0;
    virtual ~CCtSolver() { };
  protected:
    FSCommunicator *fetiCom;
    Connectivity *mpcToCpu;
    int numSubsWithMpcs;
    GenSubDomain<Scalar> **subsWithMpcs;
    int glNumMpc;
};

#endif
