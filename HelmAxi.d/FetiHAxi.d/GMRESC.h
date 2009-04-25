#ifndef _GMRES_C_H_
#define _GMRES_C_H_

#include <Utils.d/MyComplex.h>


class TaskDescr;
class ThreadLock;
class GMRESCOp;


class GMRESC {

    int len;
    int numV;
    int maxV;
    DComplex *givensC;
    DComplex *givensS;
    DComplex *g;
    DComplex *y;
    DComplex *matrixH;

    DComplex *op1;
    DComplex *op2, *op3;

    void (GMRESCOp::*operation)();
    int numTasks;
    TaskDescr **oos;
    ThreadLock lock; // need a lock

 public:

    GMRESC(int _len, int maxsize);
    ~GMRESC();
    void generateRotation(DComplex a, DComplex b, DComplex &cs, DComplex &ss);
    void applyRotation(DComplex &a, DComplex &b, DComplex cs, DComplex ss);
    void reInit();
    void init(DComplex *v0, double beta);
    double orthoAdd(DComplex *Fv, DComplex *v);
    void solution(DComplex *u);
    int numDir() { return numV; }

    friend class GMRESCOp;
};


class GMRESCOp : public TaskDescr {

    GMRESC *os;
    DComplex *locAllV;
    int loclen;
    int numV;
    int index;

  public:

    GMRESCOp(GMRESC*, int, int);
    void reInit();
    void run();
    void addVec();
    void addVecAndNorm();
    void dot();
    void mult();
    void multAdd();
    void test();
};
#endif
