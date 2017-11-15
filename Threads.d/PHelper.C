#include <Timers.d/GetTime.h>

template <typename FType>
class FunctorExecuter : public TaskDescr {
public:
	FunctorExecuter(const FType &ft) : fctor(ft) {}

	FunctorExecuter(const FunctorExecuter &) = default;

	void runFor(int) override;
private:
	const FType &fctor; //!< The functor to execute.
};

template <typename FType>
void
FunctorExecuter<FType>::runFor(int i) {
	fctor(i);
}

template <typename FType>
auto makeExecuter(const FType &ftor) {
	return FunctorExecuter<FType>{ftor};
}

template <typename TA, typename TB, typename ... FArgs, typename ... Args>
void execParal(int n, TA *target, void (TB::*f)(int, FArgs ... fa) const, Args &&...args)
{
	auto call =[&](int i) { (static_cast<const TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};

template <typename TA, typename TB, typename ... FArgs, typename ... Args>
void execParal(int n, TA *target, void (TB::*f)(int, FArgs ... fa), Args &&...args)
{
	auto call =[&](int i) { (static_cast<TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	threadManager->execParal(n, &fe);
};

template <typename TA, typename TB, typename ... FArgs, typename ... Args>
void timedParal(DistTimer &timer, int n, TA *target, void (TB::*f)(int, FArgs ... fa), Args &&...args)
{
	double initTime = getTime();
	long initMem  = threadManager->memoryUsed();
	auto call =[&](int i) { (static_cast<TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	threadManager->execTimedParal(timer, n, &fe);
	timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
};


template <typename TA, typename TB, typename ... FArgs, typename ... Args, typename X = typename std::enable_if<sizeof...(FArgs) == sizeof...(Args)>::type>
void timedParal(DistTimer &timer, int n, TA *target, void (TB::*f)(int, FArgs ... fa) const, Args &&...args)
{
	auto call =[&](int i) { (static_cast<const TB *>(target)->*f)(i, std::forward<Args>(args)...); };
	auto fe = makeExecuter(call);
	double initTime = getTime();
	long initMem  = threadManager->memoryUsed();
	threadManager->execTimedParal(timer, n, &fe);
	timer.addOverAll( threadManager->memoryUsed()-initMem, getTime()-initTime );
};


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


