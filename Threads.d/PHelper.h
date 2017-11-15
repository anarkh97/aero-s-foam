#ifndef _PHELPER_H_
#define _PHELPER_H_

#include <Threads.d/Paral.h>
#include <Timers.d/DistTimer.h>

// 4 pointer and two object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH>
void timedParal4P2(DistTimer &timer, int n, TA *, void (TB::*f)(int, TC*, TD*,TE*,TF*, TG,TH), TC*, TD*, TE*, TF*, TG,TH);

// 3 pointer and one object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void timedParal4(DistTimer &timer, int n, TA *, void (TA::*f)(int, TB*, TC*, TD*, TE), TB*, TC*, TD*, TE);

// six pointer arguments execute in parallel

template <class TA, class TB>
void paralApplyToAll(int n, TA **, void (TB::*f)());

template <class TA>
void timedParalApplyToAll(DistTimer &timer, int n, TA **, void (TA::*f)());

template <class TA, class TB>
void paralApplyToAll(int n, TA *,  void (TA::*f)(TB), TB);

template <class TA, class TB>
void timedParalApplyToAll(DistTimer &timer, int n, TA *, void (TA::*f)(TB), TB);

template <class TA, class TB>
void paralApplyToAll(int n, TA **, void (TA::*f)(TB), TB);

template <class TA, class TB>
void timedParalApplyToAll(DistTimer &timer,int n, TA **, void (TA::*f)(TB), TB);

template <class TA, class TB, class TC>
void paralApplyToAll(int n, TA  *, void (TA::*f)(TB, TC), TB, TC);

template <class TA, class TB, class TC>
void paralApplyToAll(int n, TA **, void (TA::*f)(TB, TC), TB, TC);

// PJSA 2-24-03 
// subEval is a new helper function used in rewrite of paralApply to eliminate propriatory code
template <class A, class B> class Temp { };

template <class T> class Temp<T, T> { 
  public:
    inline static T subEval(T z, int) { return z; }
};
template <class T> class Temp<T, T*> {
  public:
    inline static T subEval(T *z, int i) { return z[i]; }
};
template <class T> class Temp<T&, T&> { 
  public:
    inline static T& subEval(T &z, int) { return z; }
};
template <class T> class Temp<T&, T*> {
  public:
    inline static T& subEval(T *z, int) { return *z; }
};
template <class T> class Temp<T&, T> {
  public:
    inline static T& subEval(T &z, int) { return z; }
};

/*
template <class T> 
inline T subEval(T z, int) { return z; }

template <class T> 
inline T subEval(T *z, int i) { return z[i]; }
*/

template <class A, class B>
void paralApply(int n, A **a, void(B::*fct)());

template <class A, class B, class Arg1, class PArg1>
void paralApply(int n, A **a, void(B::*fct)(Arg1), PArg1 arg1);

template <class A, class B,
          class Arg1,  class Arg2,
          class PArg1, class PArg2>
void paralApply(int n, A **a, void(B::*fct)(Arg1, Arg2), PArg1 arg1, PArg2 arg2);

template <class A, class B,
          class Arg1,  class Arg2,  class Arg3,
          class PArg1, class PArg2, class PArg3>
void paralApply(int n, A **a, void(B::*fct)(Arg1, Arg2, Arg3),
                PArg1 arg1, PArg2 arg2, PArg3 arg3);

template <class A, class B,
          class Arg1,  class Arg2,  class Arg3,  class Arg4,
          class PArg1, class PArg2, class PArg3, class PArg4>
void paralApply(int n, A **a, void(B::*fct)(Arg1, Arg2, Arg3, Arg4),
                     PArg1 arg1, PArg2 arg2, PArg3 arg3, PArg4 arg4);

template <class A, class B,
          class Arg1,  class Arg2,  class Arg3,  class Arg4,  class Arg5,
          class PArg1, class PArg2, class PArg3, class PArg4, class PArg5>
void paralApply(int n, A **a, void(B::*fct)(Arg1, Arg2, Arg3, Arg4, Arg5),
                PArg1 arg1, PArg2 arg2, PArg3 arg3, PArg4 arg4, PArg5 arg5);
 
template <class A, class B,
          class Arg1,  class Arg2,  class Arg3,  class Arg4,  class Arg5,  class Arg6,
          class PArg1, class PArg2, class PArg3, class PArg4, class PArg5, class PArg6>
void paralApply(int n, A **a, void(B::*fct)(Arg1, Arg2, Arg3, Arg4, Arg5, Arg6),
             PArg1 arg1, PArg2 arg2, PArg3 arg3, PArg4 arg4, PArg5 arg5, PArg6 arg6);



#ifdef _TEMPLATE_FIX_
#include <Threads.d/PHelper.C>
#endif

#endif
