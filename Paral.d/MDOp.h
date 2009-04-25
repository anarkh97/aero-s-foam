#ifndef _MD_OP_H_
#define _MD_OP_H_

//#include <Driver.d/SubDomain.h>
//#include <Feti.d/DistrVector.h>

template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;
template <class Scalar> class GenCuCSparse;
typedef GenCuCSparse<double> CuCSparse;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
class DistrGeomState;

class MultiDomainOp : public TaskDescr {
    SubDomain **sd;
    CuCSparse **Kuc;
    DistrVector *v1, *v2, *v3, *v4;
    double c1;    // obviously c1 refers to the time, duh
    double *userDefDisps;
    DistrGeomState *geomState;

    void  (MultiDomainOp::*f)(int);
 public:
    MultiDomainOp(void (MultiDomainOp::*_f)(int),  SubDomain **_sd,
                  DistrVector *_v1, DistrVector*_v2, double c, CuCSparse **Kuc);
    MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **,
         DistrVector *, DistrVector*, double, double *userDefDisps = 0, DistrGeomState *geomState = 0);
    MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **,
                  DistrVector *, DistrVector*, DistrVector*, DistrVector*);
    MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **,
                  DistrVector *, DistrVector*, DistrVector*);
    MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **,
                  DistrVector *);
    MultiDomainOp(void (MultiDomainOp::*_f)(int), SubDomain **);

    void computeExtForce(int);
    void getConstForce(int);
    void getInitState(int);
    void runFor(int);
    void makeAllDOFs(int);

};

#endif
