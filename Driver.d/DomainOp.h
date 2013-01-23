#ifndef _DOMAIN_OP_H_
#define _DOMAIN_OP_H_

#include <Threads.d/Paral.h>

template <class Scalar> class GenSubDomain;
typedef GenSubDomain<double> SubDomain;

class DomainTask : public TaskDescr {
    SubDomain *sd;
    Domain *g;
    int subn;
    Connectivity *connect, *subToN;
    void (SubDomain::*op)();
  public:
    DomainTask() {}
    DomainTask(Domain *, int i, Connectivity *, Connectivity *);
    DomainTask(SubDomain *_sd, void (SubDomain::*_op)()) { sd = _sd; op = _op; }
    SubDomain *getSD() { return sd; }
    void run();
};

#endif
