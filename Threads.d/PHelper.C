#include <Timers.d/GetTime.h>

template <class TA, class TB>
class ZeroArgExecuter : public TaskDescr {
    TA *target;
    void (TB::*f)(int);
    void (TB::*g)();
  public:
    ZeroArgExecuter(TA *t, void (TB::*_f)(int)) : target(t) { f = _f;}
    ZeroArgExecuter(TA *t, void (TB::*_g)() )   : target(t) { g = _g;}
    void run();
    void runFor(int);
};

template <class TA, class TB>
void
ZeroArgExecuter<TA, TB>::run()
{
 TB *tg = static_cast<TB *>(target);
 (tg->*g)();
}

template <class TA, class TB>
void
ZeroArgExecuter<TA, TB>::runFor(int i)
{
 TB *tg = static_cast<TB *>(target);
 (tg->*f)(i);
}

template <class TA, class TB>
void execParal(int n, TA *target, void (TB::*f)(int) )
{
 ZeroArgExecuter<TA, TB> ze(target,f);
 threadManager->execParal(n, &ze);
}

template <class TA, class TB>
void timedParal(DistTimer &timer, int n, TA *target, void (TB::*f)(int) )
{
 ZeroArgExecuter<TA, TB> ze(target, f);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &ze);
 long finalMem = threadManager->memoryUsed();
 timer.addOverAll( finalMem - initMem, getTime()-initTime );
}


template <class TA, class TB, class TC>
class OneArgExecuter : public TaskDescr {
    TA *target;
    void (TB::*f)(int, TC);
    TC c;
  public:
    OneArgExecuter(TA *t, void (TB::*_f)(int, TC), TC _c) : target(t), c(_c) { f = _f;}
    void runFor(int);
};

template <class TA, class TB, class TC>
void
OneArgExecuter<TA, TB, TC>::runFor(int i)
{
 TB *tg = static_cast<TB *>(target);
 (tg->*f)(i, c);
}

template <class TA, class TB, class TC>
void execParal(int n, TA *target, void (TB::*f)(int, TC), TC c)
{
 OneArgExecuter<TA,TB,TC> oe(target,f,c);
 threadManager->execParal(n, &oe);
}


// one pointer argument execute in parallel
template <class TA, class TB, class TC>
void execParal(int n, TA *target, void (TB::*f)(int, TC*), TC *c)
{
 OneArgExecuter<TA,TB,TC*> oe(target,f,c);
 threadManager->execParal(n, &oe);
}

//HB: one reference argument execute in parallel
template <class TA, class TB, class TC>
void execParal(int n, TA *target, void (TB::*f)(int, TC&), TC &c)
{
 OneArgExecuter<TA,TB,TC&> oe(target,f,c);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC>
void timedParal(DistTimer &timer, int n, 
                TA *target, void (TB::*f)(int, TC), TC c)
{
 OneArgExecuter<TA,TB,TC> oe(target,f,c);

 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &oe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC>
void execParal1R(int n, TA *target, void (TB::*f)(int, TC&), TC &c)
{
 OneArgExecuter<TA,TB,TC&> oe(target,f,c);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC>
void execParal1R(int n, TA *target, void (TB::*f)(int, TC*), TC *c)
{
 OneArgExecuter<TA,TB,TC*> oe(target,f,c);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC>
void timedParal1R(DistTimer &timer, int n,
                  TA *target,void (TB::*f)(int, TC&), TC &c)
{
 OneArgExecuter<TA,TB,TC&> oe(target,f,c);

 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();

 threadManager->execTimedParal(timer, n, &oe);

 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}


template <class TA, class TB, class TC, class TD>
class TwoArgExecuter : public TaskDescr {
    TA *target;
    void (TB::*f)(int, TC, TD);
    TC c;
    TD d;
  public:
    TwoArgExecuter(TA *t, void (TB::*_f)(int, TC, TD), TC _c, TD _d) : 
       target(t), c(_c), d(_d) { f = _f;}
    void runFor(int);
};

template <class TA, class TB, class TC, class TD>
void
TwoArgExecuter<TA, TB, TC, TD>::runFor(int i)
{
 TB *tg = static_cast<TB *>(target);
 (tg->*f)(i, c, d);
}

template <class TA, class TB, class TC, class TD>
void execParal(int n, TA *target, void (TB::*f)(int, TC, TD), TC c, TD d)
{
 TwoArgExecuter<TA,TB,TC,TD> oe(target,f,c,d);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD>
void timedParal(DistTimer &timer, int n,
                TA *target, void (TB::*f)(int, TC, TD), TC c, TD d)
{
 TwoArgExecuter<TA,TB,TC,TD> oe(target,f,c,d);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &oe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD>
void execParal2R(int n, TA *target, void (TB::*f)(int, TC&, TD&), TC &c, TD &d)
{
 TwoArgExecuter<TA,TB,TC&,TD&> oe(target,f,c,d);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD>
void execParal2R(int n, TA *target, void (TB::*f)(int, TC&, TD*), TC &c, TD *d)
{
 TwoArgExecuter<TA,TB,TC&,TD*> oe(target,f,c,d);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD>
void execParal2R(int n, TA *target, void (TB::*f)(int, TC&, TD), TC &c, TD d)
{
 TwoArgExecuter<TA,TB,TC&,TD> oe(target,f,c,d);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD>
void execParal2R(int n, TA *target, void (TB::*f)(int, TC*, TD*), TC *c, TD *d)
{
 TwoArgExecuter<TA,TB,TC*,TD*> oe(target,f,c,d);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD>
void timedParal2R(DistTimer &timer, int n,
                  TA *target, void (TB::*f)(int, TC&, TD&), TC &c, TD &d)
{
 TwoArgExecuter<TA,TB,TC&,TD&> oe(target,f,c,d);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();

 threadManager->execTimedParal(timer, n, &oe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD>
void timedParal2R(DistTimer &timer, int n,
                  TA *target, void (TB::*f)(int, TC*, TD*), TC *c, TD *d)
{
 TwoArgExecuter<TA,TB,TC*,TD*> oe(target,f,c,d);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();

 threadManager->execTimedParal(timer, n, &oe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

//HB
template <class TA, class TB, class TC, class TD>
void timedParal2R(DistTimer &timer, int n,
                  TA *target, void (TB::*f)(int, TC&, TD*), TC &c, TD *d)
{
 TwoArgExecuter<TA,TB,TC&,TD*> oe(target,f,c,d);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();

 threadManager->execTimedParal(timer, n, &oe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}


template <class TA, class TB, class TC, class TD, class TE>
class ThreeArgExecuter : public TaskDescr {
    TA *target;
    void (TB::*f)(int, TC, TD, TE);
    TC c;
    TD d;
    TE e;
 public:
   ThreeArgExecuter(TA *t, void (TB::*_f)(int, TC, TD, TE), TC _c, TD _d, TE _e) : target(t), c(_c), d(_d), e(_e)
                            { f = _f;}
    void runFor(int);
};

template <class TA, class TB, class TC, class TD, class TE>
void
ThreeArgExecuter<TA,TB,TC,TD, TE>::runFor(int i)
{
 TB *tg = static_cast<TB *>(target);
 (tg->*f)(i, c, d, e);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal(int n, TA *target, void (TB::*f)(int, TC, TD, TE), TC c, TD d, TE e)
{
 ThreeArgExecuter<TA,TB,TC,TD,TE> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void timedParal(DistTimer & timer, int n, TA *target,
                void (TB::*f)(int, TC, TD, TE), TC c, TD d, TE e)
{
 ThreeArgExecuter<TA,TB,TC,TD,TE> oe(target,f,c,d,e);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &oe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC&, TD&, TE*), TC &c, TD &d, TE *e)
{
 ThreeArgExecuter<TA,TB,TC&,TD&,TE*> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC&, TD&, TE&), TC &c, TD &d, TE &e)
{
 ThreeArgExecuter<TA,TB,TC&,TD&,TE&> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC, TD&, TE*), TC c, TD &d, TE *e)
{
 ThreeArgExecuter<TA,TB,TC,TD&,TE*> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC&, TD&, TE), TC &c, TD &d, TE e)
{
 ThreeArgExecuter<TA,TB,TC&,TD&,TE> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC&, TD, TE), TC &c, TD d, TE e)
{
 ThreeArgExecuter<TA,TB,TC&,TD,TE> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC&, TD, TE*), TC &c, TD d, TE *e)
{
 ThreeArgExecuter<TA,TB,TC&,TD,TE*> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC&, TD*, TE), TC &c, TD *d, TE e)
{
 ThreeArgExecuter<TA,TB,TC&,TD*,TE> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC&, TD, TE&), TC &c, TD d, TE &e)
{
 ThreeArgExecuter<TA,TB,TC&,TD,TE&> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC&, TD*, TE*), TC &c, TD *d, TE *e)
{
 ThreeArgExecuter<TA,TB,TC&,TD*,TE*> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC, TD*, TE*), TC c, TD *d, TE *e)
{
 ThreeArgExecuter<TA,TB,TC,TD*,TE*> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC*, TD, TE*), TC* c, TD d, TE *e)
{
 ThreeArgExecuter<TA,TB,TC*,TD,TE*> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC*, TD*, TE*), TC *c, TD *d, TE *e)
{
 ThreeArgExecuter<TA,TB,TC*,TD*,TE*> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *target, void (TB::*f)(int, TC*, TD*, TE), TC *c, TD *d, TE e)
{
 ThreeArgExecuter<TA,TB,TC*,TD*,TE> oe(target,f,c,d,e);
 threadManager->execParal(n, &oe);
}

template <class TA, class TB, class TC, class TD, class TE>
void timedParal3(DistTimer &timer, int n, TA *target, void (TB::*f)(int, TC*, TD*, TE), TC *c, TD *d, TE e)
{
 ThreeArgExecuter<TA,TB,TC*,TD*, TE> fe(target,f,c,d,e);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &fe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD, class TE>
void timedParal3P(DistTimer &timer, int n, TA *target, void (TB::*f)(int, TC*, TD*, TE*), TC *c, TD *d, TE *e)
{
 ThreeArgExecuter<TA,TB,TC*,TD*,TE*> fe(target,f,c,d,e);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &fe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD, class TE>
void timedParal(DistTimer &timer, int n, TA *target, 
                 void (TB::*f)(int, TC&, TD&, TE&), TC &c, TD &d, TE &e)
{
 ThreeArgExecuter<TA,TB,TC&,TD&,TE&> oe(target,f,c,d,e);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &oe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD, class TE>
class FourArgExecuter: public TaskDescr {
    TA *target;
    void (TA::*f)(int, TB, TC, TD, TE);
    TB b;
    TC c;
    TD d;
    TE e;
  public:
    FourArgExecuter(TA *t, void (TA::*_f)(int, TB, TC, TD, TE), TB _b, TC _c, TD _d, TE _e)
     : target(t), b(_b), c(_c), d(_d), e(_e) { f = _f; }
    void runFor(int);
};

template <class TA, class TB, class TC, class TD, class TE>
void
FourArgExecuter<TA,TB,TC,TD,TE>::runFor(int i)
{
  (target->*f)(i, b, c,d,e);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal(int n, TA *target, void (TA::*f)(int, TB, TC, TD, TE), TB b, TC c, TD d, TE e)
{
 FourArgExecuter<TA,TB,TC,TD,TE> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void timedParal(DistTimer & timer,int n, TA *target,
                void (TA::*f)(int, TB, TC, TD, TE), TB b, TC c,TD d, TE e)
{
 FourArgExecuter<TA,TB,TC,TD,TE> fe(target,f,b,c,d,e);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &fe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *target, void (TA::*f)(int, TB&, TC&, TD&, TE&), TB &b, TC &c, TD &d, TE &e)
{
 FourArgExecuter<TA,TB&,TC&,TD&,TE&> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *target, void (TA::*f)(int, TB&, TC, TD&, TE&), TB &b, TC c, TD &d, TE &e)
{
 FourArgExecuter<TA,TB&,TC,TD&,TE&> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *target, void (TA::*f)(int, TB&, TC&, TD&, TE), TB &b, TC &c, TD &d, TE e)
{
 FourArgExecuter<TA,TB&,TC&,TD&,TE> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *target, void (TA::*f)(int, TB&, TC&, TD&, TE*), TB &b, TC &c, TD &d, TE *e)
{
 FourArgExecuter<TA,TB&,TC&,TD&,TE*> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *target, void (TA::*f)(int, TB&, TC*, TD, TE*), TB &b, TC *c, TD d, TE *e)
{
 FourArgExecuter<TA,TB&,TC*,TD,TE*> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *target, void (TA::*f)(int, TB&, TC&, TD, TE), TB &b, TC &c, TD d, TE e)
{
 FourArgExecuter<TA,TB&,TC&,TD,TE> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *target, void (TA::*f)(int, TB&, TC, TD&, TE), TB &b, TC c, TD &d, TE e)
{
 FourArgExecuter<TA,TB&,TC,TD&,TE> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *target, void (TA::*f)(int, TB&, TC, TD, TE&), TB &b, TC c, TD d, TE &e)
{
 FourArgExecuter<TA,TB&,TC,TD,TE&> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *target, void (TA::*f)(int, TB&, TC, TD*, TE), TB &b, TC c, TD *d, TE e)
{
 FourArgExecuter<TA,TB&,TC,TD*,TE> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *target, void (TA::*f)(int, TB&, TC, TD, TE), TB &b, TC c, TD d, TE e)
{
 FourArgExecuter<TA,TB&,TC,TD,TE> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *target, void (TA::*f)(int, TB*, TC*, TD*, TE*), TB *b, TC *c, TD *d, TE *e)
{
 FourArgExecuter<TA,TB*,TC*,TD*,TE*> fe(target,f,b,c,d,e);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE>
void timedParal4R(DistTimer &timer, int n, TA *target, void (TA::*f)(int, TB&, TC&, TD&, TE&), TB &b, TC &c, TD &d, TE &e)
{
 FourArgExecuter<TA,TB&,TC&,TD&,TE&> fe(target,f,b,c,d,e);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &fe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD, class TE>
void timedParal4(DistTimer &timer, int n, TA *target, void (TA::*f)(int, TB*, TC*, TD*, TE), TB *b, TC *c, TD *d, TE e)
{
 FourArgExecuter<TA,TB*,TC*,TD*,TE> fe(target,f,b,c,d,e);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &fe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD, class TE, class TF>
class FiveArgExecuter: public TaskDescr {
    TA *target;
    void (TA::*fct)(int, TB, TC, TD, TE, TF);
    TB b;
    TC c;
    TD d;
    TE e;
    TF f;
  public:
    FiveArgExecuter(TA *t, void (TA::*_fc)(int, TB, TC, TD, TE, TF), 
                                          TB _b, TC _c, TD _d, TE _e, TF _f)
     : target(t), b(_b), c(_c), d(_d), e(_e), f(_f) { fct = _fc; }
    void runFor(int);
};

template <class TA, class TB, class TC, class TD, class TE, class TF>
void
FiveArgExecuter<TA,TB,TC,TD,TE,TF>::runFor(int i)
{
  (target->*fct)(i, b, c,d,e,f);
}

template <class TA, class TB, class TC, class TD, class TE, class TF> // CKT
void timedParal5R(DistTimer &timer, int n, TA *target, void (TA::*ff)(int, TB&, TC&, TD&, TE&, TF&), TB &b, TC &c, TD &d, TE &e, TF &f)
{
 FiveArgExecuter<TA,TB&,TC&,TD&,TE&,TF&> fe(target,ff,b,c,d,e,f);
 //double initTime = getTime();
 //long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *target, void (TA::*f)(int, TB&, TC&, TD&, TE&, TG&), TB &b, TC &c, TD &d, TE &e, TG &g)
{
 FiveArgExecuter<TA,TB&,TC&,TD&,TE&,TG&> fe(target,f,b,c,d,e,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *target, void (TA::*f)(int, TB&, TC, TD, TE, TG**), TB &b, TC c, TD d, TE e, TG **g)
{
 FiveArgExecuter<TA,TB&,TC,TD,TE,TG**> fe(target,f,b,c,d,e,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *target, void (TA::*f)(int, TB&, TC*, TD&, TE&, TG*), TB &b, TC *c, TD &d, TE &e, TG *g)
{
 FiveArgExecuter<TA,TB&,TC*,TD&,TE&,TG*> fe(target,f,b,c,d,e,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *target, void (TA::*f)(int, TB&, TC*, TD&, TE, TG*), TB &b, TC *c, TD &d, TE e, TG *g)
{
 FiveArgExecuter<TA,TB&,TC*,TD&,TE,TG*> fe(target,f,b,c,d,e,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *target, void (TA::*f)(int, TB&, TC&, TD&, TE, TG), TB &b, TC &c, TD &d, TE e, TG g)
{
 FiveArgExecuter<TA,TB&,TC&,TD&,TE,TG> fe(target,f,b,c,d,e,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF>
void execParal(int n, TA *target, void (TA::*fct)(int, TB, TC, TD, TE, TF), TB b, TC c,TD d, TE e, TF f)
{
 FiveArgExecuter<TA,TB,TC,TD,TE,TF> fe(target,fct, b,c,d,e,f);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *target, void (TA::*f)(int, TB*, TC*, TD*, TE, TG), TB *b, TC *c, TD *d, TE e, TG g)
{
 FiveArgExecuter<TA,TB*,TC*,TD*,TE,TG> fe(target,f,b,c,d,e,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *target, void (TA::*f)(int, TB&, TC*, TD*, TE, TG), TB &b, TC *c, TD *d, TE e, TG g)
{
 FiveArgExecuter<TA,TB&,TC*,TD*,TE,TG> fe(target,f,b,c,d,e,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF>
void timedParal(DistTimer &timer, int n, TA *target, 
     void (TA::*fct)(int, TB, TC, TD, TE, TF), TB b, TC c,TD d, TE e, TF f)
{
 FiveArgExecuter<TA,TB,TC,TD,TE,TF> fe(target,fct, b,c,d,e,f);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &fe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
class SixArgExecuter: public TaskDescr {
    TA *target;
    void (TA::*fct)(int, TB, TC, TD, TE, TF, TG);
    TB b;
    TC c;
    TD d;
    TE e;
    TF f;
    TG g;
  public:
    SixArgExecuter(TA *t, void (TA::*_fc)(int, TB, TC, TD, TE,TF,TG),TB _b, TC _c,
TD _d, TE _e, TF _f, TG _g)
     : target(t), b(_b), c(_c), d(_d), e(_e), f(_f), g(_g) { fct = _fc; }
    void runFor(int);
};

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void
SixArgExecuter<TA,TB,TC,TD,TE,TF,TG>::runFor(int i)
{
  (target->*fct)(i, b, c,d,e,f,g);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void execParal(int n, TA *target, void (TA::*fct)(int, TB, TC, TD, TE, TF,TG), TB b, TC c, TD d, TE e, TF f, TG g)
{
 SixArgExecuter<TA,TB,TC,TD,TE,TF,TG> fe(target,fct, b,c,d,e,f,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void execParal6R(int n, TA *target, void (TA::*fct)(int, TB&, TC&, TD&, TE&, TF&, TG), TB &b, TC &c, TD &d, TE &e, TF &f, TG g)
{
 SixArgExecuter<TA,TB&,TC&,TD&,TE&,TF&,TG> fe(target,fct,b,c,d,e,f,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void execParal6R(int n, TA *target, void (TA::*fct)(int, TB&, TC&, TD&, TE&, TF, TG), TB &b, TC &c, TD &d, TE &e, TF f, TG g)
{
 SixArgExecuter<TA,TB&,TC&,TD&,TE&,TF,TG> fe(target,fct,b,c,d,e,f,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void execParal6R(int n, TA *target, void (TA::*fct)(int, TB&, TC&, TD, TE, TF&, TG&), TB &b, TC &c, TD d, TE e, TF &f, TG &g)
{
 SixArgExecuter<TA,TB&,TC&,TD,TE,TF&,TG&> fe(target,fct,b,c,d,e,f,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TG, class TH>
void execParal6R(int n, TA *target, void (TA::*f)(int, TB&, TC, TD, TE, TG**, TH**), TB &b, TC c, TD d, TE e, TG **g, TH **h)
{
 SixArgExecuter<TA,TB&,TC,TD,TE,TG**,TH**> fe(target,f,b,c,d,e,g,h);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void execParal6R(int n, TA *target, void (TA::*fct)(int, TB&, TC&, TD&, TE, TF*, TG), TB &b, TC &c, TD &d, TE e, TF *f, TG g)
{
 SixArgExecuter<TA,TB&,TC&,TD&,TE,TF*,TG> fe(target,fct,b,c,d,e,f,g);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG> 
void timedParal6R(DistTimer &timer, int n, TA *target, void (TA::*ff)(int, TB&, TC&, TD&, TE&, TF&, TG*), TB &b, TC &c, TD &d, TE &e, TF &f, TG *g)
{
 SixArgExecuter<TA,TB&,TC&,TD&,TE&,TF&,TG*> fe(target,ff,b,c,d,e,f,g);
 //double initTime = getTime();
 //long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void timedParal6R(DistTimer &timer, int n, TA *target, void (TA::*ff)(int, TB&, TC&, TD&, TE&, TF&, TG&), TB &b, TC &c, TD &d, TE &e, TF &f, TG &g)
{
 SixArgExecuter<TA,TB&,TC&,TD&,TE&,TF&,TG&> fe(target,ff,b,c,d,e,f,g);
 //double initTime = getTime();
 //long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void timedParal(DistTimer &timer, int n, TA *target, void (TA::*fct)(int, TB, TC, TD, TE, TF,TG), TB b, TC c, TD d, TE e, TF f, TG g)
{
 SixArgExecuter<TA,TB,TC,TD,TE,TF,TG> fe(target,fct, b,c,d,e,f,g);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &fe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH>
class SixArgExecuterF: public TaskDescr {
    TA *target;
    void (TB::*fct)(int, TC, TD, TE, TF, TG, TH);
    TC c;
    TD d;
    TE e;
    TF f;
    TG g;
    TH h;
  public:
    SixArgExecuterF(TA *t, void (TB::*_fc)(int, TC, TD, TE,TF,TG,TH),TC _c,
TD _d, TE _e, TF _f, TG _g, TH _h)
     : target(t), c(_c), d(_d), e(_e), f(_f), g(_g), h(_h) { fct = _fc; }
    void runFor(int);
};


template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH>
void
SixArgExecuterF<TA,TB,TC,TD, TE,TF,TG,TH>::runFor(int i)
{
 TB *tg = static_cast<TB *>(target);
 (tg->*fct)(i, c, d, e,f,g,h);
}

template <class TA, class TB,  class TC, class TD, class TE, class TF, class TG, class TH>
void timedParal4P2(DistTimer &timer, int n, TA *target, void (TB::*f)(int, TC*, TD*, TE*,TF*,TG,TH), TC *c, TD *d, TE* e, TF *ff, TG g, TH h)
{
 SixArgExecuterF<TA,TB,TC*,TD*, TE*,TF*,TG,TH> fe(target,f,c,d,e,ff,g,h);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &fe);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH>
class SevenArgExecuter: public TaskDescr {
    TA *target;
    void (TA::*fct)(int, TB, TC, TD, TE, TF, TG, TH);
    TB b;
    TC c;
    TD d;
    TE e;
    TF f;
    TG g;
    TH h;
  public:
    SevenArgExecuter(TA *t, void (TA::*_fc)(int, TB, TC, TD, TE,TF,TG,TH),TB _b, TC _c,
TD _d, TE _e, TF _f, TG _g, TH _h)
     : target(t), b(_b), c(_c), d(_d), e(_e), f(_f), g(_g), h(_h) { fct = _fc; }
    void runFor(int);
};


template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH>
void
SevenArgExecuter<TA,TB,TC,TD,TE,TF,TG,TH>::runFor(int i)
{
  (target->*fct)(i, b, c,d,e,f,g,h);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH>
void execParal7R(int n, TA *target, void (TA::*fct)(int, TB&, TC, TD, TE, TF**, TG**, TH**), TB &b, TC c, TD d, TE e, TF **f, TG **g, TH **h)
{
 SevenArgExecuter<TA,TB&,TC,TD,TE,TF**,TG**,TH**> fe(target,fct, b,c,d,e,f,g,h);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH>
void execParal7R(int n, TA *target, void (TA::*fct)(int, TB&, TC&, TD&, TE&, TF&, TG&, TH&), TB &b, TC &c, TD &d, TE &e, TF &f, TG &g, TH &h)
{
 SevenArgExecuter<TA,TB&,TC&,TD&,TE&,TF&,TG&,TH&> fe(target,fct, b,c,d,e,f,g,h);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH>
void execParal7R(int n, TA *target, void (TA::*fct)(int, TB&, TC&, TD&, TE&, TF&, TG, TH), TB &b, TC &c, TD &d, TE &e, TF &f, TG g, TH h)
{
 SevenArgExecuter<TA,TB&,TC&,TD&,TE&,TF&,TG,TH> fe(target,fct, b,c,d,e,f,g,h);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH>
void execParal7R(int n, TA *target, void (TA::*fct)(int, TB&, TC&, TD, TE*, TF*, TG*, TH*), TB &b, TC &c, TD d, TE *e, TF *f, TG *g, TH *h)
{
 SevenArgExecuter<TA,TB&,TC&,TD,TE*,TF*,TG*,TH*> fe(target,fct, b,c,d,e,f,g,h);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG,
          class TH>
void execParal(int n, TA *target, void (TA::*fct)(int, TB, TC, TD, TE, TF,TG,TH),
     TB b, TC c, TD d, TE e, TF f, TG g, TH h)
{
 SevenArgExecuter<TA,TB,TC,TD,TE,TF,TG,TH> fe(target,fct, b,c,d,e,f,g,h);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH, class TI>
class EightArgExecuter: public TaskDescr {
    TA *target;
    void (TA::*fct)(int, TB, TC, TD, TE, TF, TG, TH, TI);
    TB b;
    TC c;
    TD d;
    TE e;
    TF f;
    TG g;
    TH h;
    TI i;
  public:
    EightArgExecuter(TA *t, void (TA::*_fc)(int, TB, TC, TD, TE, TF, TG, TH, TI), TB _b, TC _c,
TD _d, TE _e, TF _f, TG _g, TH _h, TI _i)
     : target(t), b(_b), c(_c), d(_d), e(_e), f(_f), g(_g), h(_h), i(_i) { fct = _fc; }
    void runFor(int);
};

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH, class TI>
void
EightArgExecuter<TA,TB,TC,TD,TE,TF,TG,TH,TI>::runFor(int _i)
{
  (target->*fct)(_i, b, c, d, e, f, g, h, i);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH, class TI>
void execParal8R(int n, TA *target, void (TA::*fct)(int, TB&, TC&, TD&, TE&, TF&, TG, TH, TI), TB &b, TC &c, TD &d, TE &e, TF &f, TG g, TH h, TI i)
{
 EightArgExecuter<TA,TB&,TC&,TD&,TE&,TF&,TG,TH,TI> fe(target,fct, b, c, d, e, f, g, h, i);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH, class TI>
void execParal8R(int n, TA *target, void (TA::*fct)(int, TB*, TC&, TD*, TE, TF*, TG*, TH*, TI*), TB *b, TC &c, TD *d, TE e, TF *f, TG *g, TH *h, TI *i)
{
 EightArgExecuter<TA,TB*,TC&,TD*,TE,TF*,TG*,TH*,TI*> fe(target,fct, b, c, d, e, f, g, h, i);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH, class TI, class TJ>
class NineArgExecuter: public TaskDescr {
    TA *target;
    void (TA::*fct)(int, TB, TC, TD, TE, TF, TG, TH, TI, TJ);
    TB b;
    TC c;
    TD d;
    TE e;
    TF f;
    TG g;
    TH h;
    TI i;
    TJ j;
  public:
    NineArgExecuter(TA *t, void (TA::*_fc)(int, TB, TC, TD, TE, TF, TG, TH, TI, TJ), TB _b, TC _c,
TD _d, TE _e, TF _f, TG _g, TH _h, TI _i, TJ _j)
     : target(t), b(_b), c(_c), d(_d), e(_e), f(_f), g(_g), h(_h), i(_i), j(_j) { fct = _fc; }
    void runFor(int);
};

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH, class TI, class TJ>
void
NineArgExecuter<TA,TB,TC,TD,TE,TF,TG,TH,TI,TJ>::runFor(int _i)
{
  (target->*fct)(_i, b, c, d, e, f, g, h, i, j);
}

template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH, class TI, class TJ>
void execParal9R(int n, TA *target, void (TA::*fct)(int, TB&, TC&, TD&, TE&, TF, TG, TH, TI, TJ), TB &b, TC &c, TD &d, TE &e, TF f, TG g, TH h, TI i, TJ j)
{
 NineArgExecuter<TA,TB&,TC&,TD&,TE&,TF,TG,TH,TI,TJ> fe(target,fct, b, c, d, e, f, g, h, i, j);
 threadManager->execParal(n, &fe);
}

template <class TA, class TB>
class NoArgApplier : public TaskDescr {
    TA **target;
    void (TB::*f)();
 public:
    NoArgApplier(TA **t, void (TB::*_f)()) : target(t) { f = _f;}
    void runFor(int);
};

template <class TA, class TB>
void
NoArgApplier<TA, TB>::runFor(int i)
{
 TB *tg = static_cast<TB *>(target[i]);
 (tg->*f)();
}

template <class TA, class TB>
void paralApplyToAll(int n, TA **a, void (TB::*f)())
{
 NoArgApplier<TA, TB> ap(a,f);
 threadManager->execParal(n, &ap);
}

template <class TA>
void timedParalApplyToAll(DistTimer &timer, int n, TA **a, void (TA::*f)())
{
 NoArgApplier<TA,TA> ap(a,f);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &ap);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB>
class OneArgApplier : public TaskDescr {
    TA **target;
    TB arg;
    void (TA::*f)(TB);
 public:
    OneArgApplier(TA **t, void (TA::*_f)(TB), TB b) : target(t) , 
          arg(b) { f = _f;}
    void runFor(int);
};

template <class TA, class TB>
void
OneArgApplier<TA,TB>::runFor(int i)
{
 (target[i]->*f)(arg);
}

template <class TA, class TB>
void paralApplyToAll(int n, TA **a, void (TA::*f)(TB), TB b)
{
 OneArgApplier<TA,TB> ap(a,f,b);
 threadManager->execParal(n, &ap);
}

template <class TA, class TB>
void timedParalApplyToAll(DistTimer &timer,int n,TA **a,void (TA::*f)(TB), TB b)
{
 OneArgApplier<TA,TB> ap(a,f,b);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &ap);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}

template <class TA, class TB>
class OneArgArrayApplier : public TaskDescr {
    TA *target;
    TB arg;
    void (TA::*f)(TB);
 public:
    OneArgArrayApplier(TA *t, void (TA::*_f)(TB), TB b) : target(t) ,
          arg(b) { f = _f;}
    void runFor(int);
};

template <class TA, class TB>
void
OneArgArrayApplier<TA,TB>::runFor(int i)
{
 (target[i].*f)(arg);
}

template <class TA, class TB>
void paralApplyToAll(int n, TA *a, void (TA::*f)(TB), TB b)
{
 OneArgArrayApplier<TA,TB> ap(a,f,b);
 threadManager->execParal(n, &ap);
}

template <class TA, class TB>
void timedParalApplyToAll(DistTimer &timer,int n,TA *a, void (TA::*f)(TB), TB b)
{
 OneArgArrayApplier<TA,TB> ap(a,f,b);
 double initTime = getTime();
 long initMem  = threadManager->memoryUsed();
 threadManager->execTimedParal(timer, n, &ap);
 timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
}


// KHP: 3-20-98
// Two argument applier template class
// Where TB and TC are also arrays of objects
// i.e. when a subdomain function wants to be run in parallel
// and the arguments are per subdomain also

template <class TA, class TB, class TC>
class TwoArgApplier : public TaskDescr {
    TA **target;
    TB argB;
    TC argC;
    void (TA::*f)(TB);
 public:
    TwoArgApplier(TA **t, void (TA::*_f)(TB), TB b, TC c) : target(t) ,
          argB(b), argC(c) { f = _f;}
    void runFor(int);
};

template <class TA, class TB, class TC>
void
TwoArgApplier<TA,TB,TC>::runFor(int i)
{
 (target[i]->*f)(argB, argC);
}

template <class TA, class TB, class TC>
void paralApplyToAll(int n, TA **a, void (TA::*f)(TB, TC), TB b, TC c)
{
 TwoArgApplier<TA,TB,TC> ap(a,f,b,c);
 threadManager->execParal(n, &ap);
}


// KHP 3-20-98 Two argument array applier template class

template <class TA, class TB, class TC>
class TwoArgArrayApplier : public TaskDescr {
    TA *target;
    TB argB;
    TC argC;
    void (TA::*f)(TB,TC);
 public:
    TwoArgArrayApplier(TA *t, void (TA::*_f)(TB,TC), TB b, TC c) : target(t), argB(b), argC(c) { f = _f;}
    void runFor(int);
};

template <class TA, class TB, class TC>
void
TwoArgArrayApplier<TA,TB,TC>::runFor(int i)
{
 (target[i].*f)(argB, argC);
}

template <class TA, class TB, class TC>
void paralApplyToAll(int n, TA *a, void (TA::*f)(TB,TC), TB b, TC c)
{
 TwoArgArrayApplier<TA,TB,TC> ap(a,f,b,c);
 threadManager->execParal(n, &ap);
}

// PJSA: 2-24-02 rewrite of paralApply functions to remove propriatory code

// no argument ParalApply

template <class A, class B>
class ParalApplyNoArgApplier : public TaskDescr {
    A **target;
    void (B::*fct)();
  public:
    ParalApplyNoArgApplier(A **a, void (B::*_f)()) : target(a), fct(_f) { }
    void runFor(int);
};

template <class A, class B>
void ParalApplyNoArgApplier<A, B>::runFor(int i)
{
  B *tg = static_cast<B *>(target[i]);
  (tg->*fct)();
}

template <class A, class B>
void paralApply(int n, A **a, void (B::*fct)())
{
  ParalApplyNoArgApplier<A,B> ap(a, fct);
  threadManager->execParal(n, &ap);
}

// one argument ParalApply

template <class A, class B, class Arg1, class PArg1>
class ParalApplyOneArgApplier : public TaskDescr {
    A **target;
    void (B::*fct)(Arg1);
    PArg1 arg1;
  public:
    ParalApplyOneArgApplier(A **a, void (B::*_f)(Arg1), PArg1 _arg1) :
              target(a), fct(_f), arg1(_arg1) { }
    void runFor(int);
};

template <class A, class B, class Arg1, class PArg1>
void ParalApplyOneArgApplier<A, B, Arg1, PArg1>::runFor(int i)
{
  B *tg = static_cast<B *>(target[i]);
  (tg->*fct)(Temp<Arg1, PArg1>::subEval(arg1,i));
}

template <class A, class B, class Arg1, class PArg1>
void paralApply(int n, A **a, void(B::*fct)(Arg1), PArg1 arg1)
{
 ParalApplyOneArgApplier<A, B, Arg1, PArg1> ap(a, fct, arg1);
 threadManager->execParal(n, &ap);
}

// two argument ParalApply

template <class A, class B, class Arg1, class Arg2, class PArg1, class PArg2>
class ParalApplyTwoArgApplier : public TaskDescr {
    A **target;
    void (B::*fct)(Arg1, Arg2);
    PArg1 arg1;
    PArg2 arg2;
  public:
    ParalApplyTwoArgApplier(A **a, void (B::*_f)(Arg1, Arg2), PArg1 _arg1, PArg2 _arg2) :
              target(a), fct(_f), arg1(_arg1), arg2(_arg2) { }
    virtual ~ParalApplyTwoArgApplier() { }
    void runFor(int);
};

template <class A, class B, class Arg1, class Arg2, class PArg1, class PArg2>
void ParalApplyTwoArgApplier<A, B, Arg1, Arg2, PArg1, PArg2>::runFor(int i)
{
  B *tg = static_cast<B *>(target[i]);
  (tg->*fct)(Temp<Arg1, PArg1>::subEval(arg1,i), Temp<Arg2, PArg2>::subEval(arg2,i));
}

template <class A, class B, class Arg1, class Arg2, class PArg1, class PArg2>
void paralApply(int n, A **a, void(B::*fct)(Arg1, Arg2), PArg1 arg1, PArg2 arg2)
{
  ParalApplyTwoArgApplier<A, B, Arg1, Arg2, PArg1, PArg2> ap(a, fct, arg1, arg2);
  threadManager->execParal(n, &ap);
}
 
// three argument ParalApply 
 
template <class A, class B, class Arg1, class Arg2, class Arg3, class PArg1, class PArg2, class PArg3>
class ParalApplyThreeArgApplier : public TaskDescr {
    A **target;
    void (B::*fct)(Arg1, Arg2, Arg3);
    PArg1 arg1;
    PArg2 arg2;
    PArg3 arg3;
  public:
    ParalApplyThreeArgApplier(A **a, void(B::*_f)(Arg1, Arg2, Arg3),
              PArg1 _arg1, PArg2 _arg2, PArg3 _arg3) :
              target(a), fct(_f), arg1(_arg1), arg2(_arg2), arg3(_arg3) { }
    void runFor(int);
};
 
template <class A, class B, class Arg1, class Arg2, class Arg3, class PArg1, class PArg2, class PArg3>
void ParalApplyThreeArgApplier<A, B, Arg1, Arg2, Arg3, PArg1, PArg2, PArg3>::runFor(int i)
{
  B *tg = static_cast<B *>(target[i]);
  (tg->*fct)(Temp<Arg1, PArg1>::subEval(arg1,i), Temp<Arg2, PArg2>::subEval(arg2,i),
                    Temp<Arg3, PArg3>::subEval(arg3,i));
}
 
template <class A, class B, class Arg1, class Arg2, class Arg3, class PArg1, class PArg2, class PArg3>
void paralApply(int n, A **a, void(B::*fct)(Arg1, Arg2, Arg3), PArg1 arg1, PArg2 arg2, PArg3 arg3)
{
  ParalApplyThreeArgApplier<A, B, Arg1, Arg2, Arg3, PArg1, PArg2, PArg3> ap(a, fct, arg1, arg2, arg3);
  threadManager->execParal(n, &ap);
}
 
// four argument ParalApply
 
template <class A, class B, class Arg1,  class Arg2,  class Arg3,  class Arg4,
                   class PArg1, class PArg2, class PArg3, class PArg4>
class ParalApplyFourArgApplier : public TaskDescr {
    A **target;
    void (B::*fct)(Arg1, Arg2, Arg3, Arg4);
    PArg1 arg1;
    PArg2 arg2;
    PArg3 arg3;
    PArg4 arg4;
  public:
    ParalApplyFourArgApplier(A **a, void(B::*_f)(Arg1, Arg2, Arg3, Arg4),
              PArg1 _arg1, PArg2 _arg2, PArg3 _arg3, PArg4 _arg4) :
              target(a), fct(_f), arg1(_arg1), arg2(_arg2), arg3(_arg3), arg4(_arg4) { }
    void runFor(int);
};
 
template <class A, class B, class Arg1,  class Arg2,  class Arg3,  class Arg4,
                   class PArg1, class PArg2, class PArg3, class PArg4>
void ParalApplyFourArgApplier<A, B, Arg1,  Arg2,  Arg3,  Arg4,
                   PArg1, PArg2, PArg3, PArg4>::runFor(int i)
{
  B *tg = static_cast<B *>(target[i]);
  (tg->*fct)(Temp<Arg1, PArg1>::subEval(arg1,i), Temp<Arg2, PArg2>::subEval(arg2,i),
                    Temp<Arg3, PArg3>::subEval(arg3,i), Temp<Arg4, PArg4>::subEval(arg4,i));
}
 
template <class A, class B, class Arg1,  class Arg2,  class Arg3,  class Arg4,
                   class PArg1, class PArg2, class PArg3, class PArg4>
void paralApply(int n, A **a, void(B::*fct)(Arg1, Arg2, Arg3, Arg4),
                PArg1 arg1, PArg2 arg2, PArg3 arg3, PArg4 arg4)
{
  ParalApplyFourArgApplier<A, B, Arg1, Arg2, Arg3, Arg4,
            PArg1, PArg2, PArg3, PArg4> ap(a, fct, arg1, arg2, arg3, arg4);
  threadManager->execParal(n, &ap);
}
 
// five argument ParalApply
 
template <class A, class B, class Arg1, class Arg2,  class Arg3,  class Arg4,  class Arg5,
                   class PArg1,class PArg2, class PArg3, class PArg4, class PArg5>
class ParalApplyFiveArgApplier : public TaskDescr {
    A **target;
    void (B::*fct)(Arg1, Arg2, Arg3, Arg4, Arg5);
    PArg1 arg1;
    PArg2 arg2;
    PArg3 arg3;
    PArg4 arg4;
    PArg5 arg5;
  public:
    ParalApplyFiveArgApplier(A **a, void(B::*_f)(Arg1, Arg2, Arg3, Arg4, Arg5),
              PArg1 _arg1, PArg2 _arg2, PArg3 _arg3, PArg4 _arg4, PArg5 _arg5) :
              target(a), fct(_f), arg1(_arg1), arg2(_arg2), arg3(_arg3), arg4(_arg4), arg5(_arg5) { }
    void runFor(int);
};
 
template <class A, class B, class Arg1,  class Arg2,  class Arg3,  class Arg4,  class Arg5,
                   class PArg1, class PArg2, class PArg3, class PArg4, class PArg5>
void ParalApplyFiveArgApplier<A, B, Arg1, Arg2, Arg3, Arg4, Arg5,
                   PArg1, PArg2, PArg3, PArg4, PArg5>::runFor(int i)
{
  B *tg = static_cast<B *>(target[i]);
  (tg->*fct)(Temp<Arg1, PArg1>::subEval(arg1,i), Temp<Arg2, PArg2>::subEval(arg2,i),
                    Temp<Arg3, PArg3>::subEval(arg3,i), Temp<Arg4, PArg4>::subEval(arg4,i),
                    Temp<Arg5, PArg5>::subEval(arg5,i));
}
 
template <class A, class B,
          class Arg1,  class Arg2,  class Arg3,  class Arg4,  class Arg5,
          class PArg1, class PArg2, class PArg3, class PArg4, class PArg5>
void paralApply(int n, A **a, void(B::*fct)(Arg1, Arg2, Arg3, Arg4, Arg5),
                PArg1 arg1, PArg2 arg2, PArg3 arg3, PArg4 arg4, PArg5 arg5)
{
  ParalApplyFiveArgApplier<A, B, Arg1, Arg2, Arg3, Arg4, Arg5,
            PArg1, PArg2, PArg3, PArg4, PArg5> ap(a, fct, arg1, arg2, arg3, arg4, arg5);
  threadManager->execParal(n, &ap);
}

// six argument ParalApply

template <class A, class B, class Arg1,  class Arg2,  class Arg3,  class Arg4,  class Arg5,  class Arg6,
                   class PArg1, class PArg2, class PArg3, class PArg4, class PArg5, class PArg6>
class ParalApplySixArgApplier : public TaskDescr {
    A **target;
    void (B::*fct)(Arg1, Arg2, Arg3, Arg4, Arg5, Arg6);
    PArg1 arg1;
    PArg2 arg2;
    PArg3 arg3;
    PArg4 arg4;
    PArg5 arg5;
    PArg6 arg6;
  public:
    ParalApplySixArgApplier(A **a, void(B::*_f)(Arg1, Arg2, Arg3, Arg4, Arg5, Arg6),
              PArg1 _arg1, PArg2 _arg2, PArg3 _arg3, PArg4 _arg4, PArg5 _arg5, PArg6 _arg6) :
              target(a), fct(_f), arg1(_arg1), arg2(_arg2), arg3(_arg3), arg4(_arg4), arg5(_arg5), arg6(_arg6) { }
    void runFor(int);
};

template <class A, class B, class Arg1,  class Arg2,  class Arg3,  class Arg4,  class Arg5,  class Arg6,
                   class PArg1, class PArg2, class PArg3, class PArg4, class PArg5, class PArg6>
void ParalApplySixArgApplier<A, B, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6,
                   PArg1, PArg2, PArg3, PArg4, PArg5, PArg6>::runFor(int i)
{
  B *tg = static_cast<B *>(target[i]);
  (tg->*fct)(Temp<Arg1, PArg1>::subEval(arg1,i), Temp<Arg2, PArg2>::subEval(arg2,i), 
                    Temp<Arg3, PArg3>::subEval(arg3,i), Temp<Arg4, PArg4>::subEval(arg4,i),
                    Temp<Arg5, PArg5>::subEval(arg5,i), Temp<Arg6, PArg6>::subEval(arg6,i));
}

template <class A, class B,
          class Arg1,  class Arg2,  class Arg3,  class Arg4,  class Arg5,  class Arg6,
          class PArg1, class PArg2, class PArg3, class PArg4, class PArg5, class PArg6>
void paralApply(int n, A **a, void(B::*fct)(Arg1, Arg2, Arg3, Arg4, Arg5, Arg6),
                PArg1 arg1, PArg2 arg2, PArg3 arg3, PArg4 arg4, PArg5 arg5, PArg6 arg6)
{
  ParalApplySixArgApplier<A, B, Arg1, Arg2, Arg3, Arg4, Arg5, Arg6,
            PArg1, PArg2, PArg3, PArg4, PArg5, PArg6> ap(a, fct, arg1, arg2, arg3, arg4, arg5, arg6);
  threadManager->execParal(n, &ap);
}


