#ifndef _GCR_C_H_
#define _GCR_C_H_

#include <Utils.d/MyComplex.h>

class DistTimer;
class TaskDescr;
class ThreadLock;
class GCRCOp;

class GCRC {
    DComplex *allFPFiP;
    int len;
    int numP;
    int maxP;

    DComplex *op1;
    DComplex *op2, *op3;

    void (GCRCOp::*operation)();
    int numTasks;
    TaskDescr **oos;
    ThreadLock lock; // need a lock
 public:
    GCRC(int _len, int maxsize);
    void orthoAdd(DComplex *, DComplex *, DComplex);
    void orthoAddTimed(DistTimer &, DComplex *, DComplex *, DComplex);
    void orthogonalize(DComplex *, DComplex *, DComplex *, DComplex *);
    void orthogonalizeTimed(DistTimer &, DComplex *, DComplex *, 
                            DComplex *, DComplex *);
    void predict(DComplex *, DComplex *);
    int numDir() { return numP; }

    friend class GCRCOp;
};

class GCRCOp : public TaskDescr {
    GCRC *os;
    DComplex *locAllP, *locAllFiP;
    int loclen;
    int numP;
    int index;
  public:
    GCRCOp(GCRC *, int, int);
    void addVec();
    void dot();
    void Fdot();
    void mult();
    void multAdd();
    void multFAdd();

    void run();
	void runFor(int) override { throw "Illegal operation called on GCRCOp"; }
};
#endif
