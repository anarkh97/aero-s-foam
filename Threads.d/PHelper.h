#ifndef _PHELPER_H_
#define _PHELPER_H_

#include <Threads.d/Paral.h>
#include <Timers.d/DistTimer.h>

// zero arguments execute in parallel
template <class TA, class TB>
void execParal(int n, TA *, void (TB::*f)(int));

// zero arguments timed execute in parallel 
template <class TA, class TB>
void timedParal(DistTimer &, int n, TA *, void (TB::*f)(int));

// one argument execute in parallel
template <class TA, class TB, class TC>
void execParal(int n, TA *, void (TB::*f)(int, TC), TC);

// one pointer argument execute in parallel
template <class TA, class TB, class TC>
void execParal(int n, TA *, void (TB::*f)(int, TC*), TC*);

//HB: one reference argument execute in parallel
template <class TA, class TB, class TC>
void execParal(int n, TA *, void (TB::*f)(int, TC&), TC&);

// one pointer argument timed execute in parallel
template <class TA, class TB, class TC>
void timedParal(DistTimer &, int n, TA *, void (TB::*f)(int, TC), TC);

// one reference argument execute in parallel
template <class TA, class TB, class TC>
void execParal1R(int n, TA *, void (TB::*f)(int, TC&), TC&);

// one pointer argument execute in parallel
template <class TA, class TB, class TC>
void execParal1R(int n, TA *, void (TB::*f)(int, TC*), TC*);

// one reference argument timed execute in parallel
template <class TA, class TB, class TC>
void timedParal1R(DistTimer &, int n, TA *, void (TB::*f)(int, TC&), TC&);

// two pointer arguments execute in parallel 
template <class TA, class TB, class TC, class TD>
void execParal(int n, TA *, void (TB::*f)(int, TC, TD), TC, TD);

// two pointer arguments timed execute in parallel 
template <class TA, class TB, class TC, class TD>
void timedParal(DistTimer &, int n, TA *, void (TB::*f)(int, TC, TD), TC, TD);

// one reference argument and one pointer argument execute in parallel
template <class TA, class TB, class TC, class TD>
void execParal2R(int n, TA *, void (TB::*f)(int, TC&, TD*), TC&, TD*);

// one reference argument and one object argument execute in parallel
template <class TA, class TB, class TC, class TD>
void execParal2R(int n, TA *, void (TB::*f)(int, TC&, TD), TC&, TD);

// two reference arguments execute in parallel
template <class TA, class TB, class TC, class TD>
void execParal2R(int n, TA *, void (TB::*f)(int, TC&, TD&), TC&, TD&);

// two pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD>
void execParal2R(int n, TA *, void (TB::*f)(int, TC*, TD*), TC*, TD*);

// two reference arguments timed execute in parallel
template <class TA, class TB, class TC, class TD>
void timedParal2R(DistTimer &,int n,TA *,void (TB::*f)(int, TC&, TD&),TC&,TD&);

// two pointer argument execute in parallel
template <class TA, class TB, class TC, class TD>
void timedParal2R(DistTimer &,int n,TA *,void (TB::*f)(int, TC*, TD*),TC*,TD*);

// one reference + one pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD>
void timedParal2R(DistTimer &,int n,TA *,void (TB::*f)(int, TC&, TD*),TC&,TD*);

// one reference + one object arguments execute in parallel
template <class TA, class TB, class TC, class TD>
void timedParal2R(DistTimer &,int n,TA *,void (TB::*f)(int, TC&, TD),TC&,TD);

// three pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE> 
void execParal(int n, TA *, void (TB::*f)(int, TC, TD, TE), TC, TD, TD, TE);

// three pointer arguments timed execute in parallel
template <class TA, class TB, class TC, class TD, class TE> 
void timedParal(DistTimer &, int n, TA *, 
                void (TB::*f)(int, TC, TD, TE), TC, TD, TE);

// two reference arguments and one pointer argument execute in parallel
template <class TA, class TB, class TC, class TD, class TE> 
void execParal3R(int n, TA *, void (TB::*f)(int, TC&, TD&, TE*), TC&, TD&, TE*);

// three reference arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE> 
void execParal3R(int n, TA *, void (TB::*f)(int, TC&, TD&, TE&), TC&, TD&, TE&);

// 1 obj, 1 reference & 1 pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *, void (TB::*f)(int, TC, TD&, TE*), TC, TD&, TE*);

// two reference arguments and one object argument execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *, void (TB::*f)(int, TC&, TD&, TE), TC&, TD&, TE);

// one reference argument and two object argument execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *, void (TB::*f)(int, TC&, TD, TE), TC&, TD, TE);

// 1 reference, 1 obj & 1 pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *, void (TB::*f)(int, TC&, TD, TE*), TC&, TD, TE*);

// 1 reference, 1 pointer & 1 obj arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *, void (TB::*f)(int, TC&, TD*, TE), TC&, TD*, TE);

// 1 reference, 1 object & 1 reference arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *, void (TB::*f)(int, TC&, TD, TE&), TC&, TD, TE&);

// 2 pointer & 1 object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *, void (TB::*f)(int, TC*, TD*, TE), TC*, TD*, TE);

// 1 reference & 2 pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *, void (TB::*f)(int, TC&, TD*, TE*), TC&, TD*, TE*);

// 2 pointer and one object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *, void (TB::*f)(int, TC*, TD, TE*), TC*, TD, TE*);

// 3 pointer objects arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal3R(int n, TA *, void (TB::*f)(int, TC*, TD*, TE*), TC*, TD*, TE*);

// 2 pointer and one object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void timedParal3(DistTimer &timer, int n, TA *, void (TB::*f)(int, TC*, TD*, TE), TC*, TD*, TE);

// 4 pointer and two object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TF, class TG, class TH>
void timedParal4P2(DistTimer &timer, int n, TA *, void (TB::*f)(int, TC*, TD*,TE*,TF*, TG,TH), TC*, TD*, TE*, TF*, TG,TH);

// 3 pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void timedParal3P(DistTimer &timer, int n, TA *, void (TB::*f)(int, TC*, TD*, TE*), TC*, TD*, TE*);

// four pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal(int n, TA *, void (TA::*f)(int, TB, TC, TD, TE), TB, TC, TD, TE);

// four object arguments timed execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void timedParal(DistTimer &timer, int n, TA *, 
                void (TA::*f)(int, TB, TC, TD, TE), TB, TC, TD, TE);

// four reference arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE&), TB&, TC&, TD&, TE&);

// 3 reference and one object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *, void (TA::*f)(int, TB&, TC, TD&, TE&), TB&, TC, TD&, TE&);

// 3 pointer and one object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void timedParal4(DistTimer &timer, int n, TA *, void (TA::*f)(int, TB*, TC*, TD*, TE), TB*, TC*, TD*, TE);

// 3 reference and one object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE), TB&, TC&, TD&, TE);

// 3 reference and one pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE*), TB&, TC&, TD&, TE*);

// 1 reference and 2 pointer and 1 object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *, void (TA::*f)(int, TB&, TC*, TD, TE*), TB&, TC*, TD, TE*);

// 2 reference and 2 object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD, TE), TB&, TC&, TD, TE);

// 2 reference and 2 object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *, void (TA::*f)(int, TB&, TC, TD&, TE), TB&, TC, TD&, TE);

// 2 reference and 2 object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *, void (TA::*f)(int, TB&, TC, TD, TE&), TB&, TC, TD, TE&);

// 1 reference and 2 object and 1 pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *, void (TA::*f)(int, TB&, TC, TD*, TE), TB&, TC, TD*, TE);

// 1 reference and 3 object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *, void (TA::*f)(int, TB&, TC, TD, TE), TB&, TC, TD, TE);

// 4 pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void execParal4R(int n, TA *, void (TA::*f)(int, TB*, TC*, TD*, TE*), TB*, TC*, TD*, TE*);

// five reference arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE&, TG&), TB&, TC&, TD&, TE&, TG&);

// 1 reference, 3 object and 1 pointer-to-pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *, void (TA::*f)(int, TB&, TC, TD, TE, TG**), TB&, TC, TD, TE, TG**);

// three reference arguments and two pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *, void (TA::*f)(int, TB&, TC*, TD&, TE&, TG*), TB&, TC*, TD&, TE&, TG*);

// four reference arguments and one pointer argument execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE&, TG*), TB&, TC&, TD&, TE&, TG*);

// two reference, two pointer and one object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *, void (TA::*f)(int, TB&, TC*, TD&, TE, TG*), TB&, TC*, TD&, TE, TG*);

// three reference arguments and two object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE, TG), TB&, TC&, TD&, TE, TG);

// one reference, two pointer and two object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal5R(int n, TA *, void (TA::*f)(int, TB&, TC*, TD*, TE, TG), TB&, TC*, TD*, TE, TG);

// three pointer and two object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG>
void execParal(int n, TA *, void (TA::*f)(int, TB*, TC*, TD*, TE, TG), TB&, TC*, TD&, TE, TG*);

// six reference arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH>
void execParal6R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE&, TG&, TH&), TB&, TC&, TD&, TE&, TG&, TH&);

// five reference arguments and one pointer argument execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH>
void execParal6R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD*, TE&, TG&, TH&), TB&, TC&, TD*, TE&, TG&, TH&);

// five reference arguments and one object argument execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH>
void execParal6R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE&, TG&, TH), TB&, TC&, TD&, TE&, TG&, TH);

// four reference arguments and two object argument execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH>
void execParal6R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE&, TG, TH), TB&, TC&, TD&, TE&, TG, TH);

// four reference arguments and two object argument execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH>
void execParal6R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD, TE, TG&, TH&), TB&, TC&, TD, TE, TG&, TH&);

// 1 reference, 3 object and 2 pointer-to-pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH>
void execParal6R(int n, TA *, void (TA::*f)(int, TB&, TC, TD, TE, TG**, TH**), TB&, TC, TD, TE, TG**, TH**);

// 3 reference, 2 object and 1 pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH>
void execParal6R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE, TG*, TH), TB&, TC&, TD&, TE, TG*, TH);

// seven reference arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH, class TI>
void execParal7R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE&, TG&, TH&, TI&), TB&, TC&, TD&, TE&, TG&, TH&, TI&);

// 1 reference, 3 object and 3 pointer-to-pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH, class TI>
void execParal7R(int n, TA *, void (TA::*f)(int, TB&, TC, TD, TE, TG**, TH**, TI**), TB&, TC, TD, TE, TG**, TH**, TI**);

// five reference and two object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH, class TI>
void execParal7R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE&, TG&, TH, TI), TB&, TC&, TD&, TE&, TG&, TH, TI);

// various reference, object and pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH, class TI>
void execParal7R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD, TE*, TG*, TH*, TI*), TB&, TC&, TD, TE*, TG*, TH*, TI*);

// five reference and three object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH, class TI, class TJ>
void execParal8R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE&, TG&, TH, TI, TJ), TB&, TC&, TD&, TE&, TG&, TH, TI, TJ);

// various reference, pointer and object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH, class TI, class TJ>
void execParal8R(int n, TA *, void (TA::*f)(int, TB*, TC&, TD*, TE, TG*, TH*, TI*, TJ*), TB*, TC&, TD*, TE, TG*, TH*, TI*, TJ*);

// four reference and five object arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TG, class TH, class TI, class TJ, class TK>
void execParal9R(int n, TA *, void (TA::*f)(int, TB&, TC&, TD&, TE&, TG, TH, TI, TJ, TK), TB&, TC&, TD&, TE&, TG, TH, TI, TJ, TK);

// four reference arguments timed execute in parallel
template <class TA, class TB, class TC, class TD, class TE>
void timedParal4R(DistTimer &timer, int n, TA *, 
                 void (TA::*f)(int, TB&, TC&, TD&, TE&), TB&, TC&, TD&, TE&);

// five reference arguments timed execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TF>
void timedParal5R(DistTimer &timer, int n, TA *,
                 void (TA::*ff)(int, TB&, TC&, TD&, TE&, TF&), TB&, TC&, TD&, TE&, TF&);

// five reference arguments + 1 pointer argument timed execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void timedParal6R(DistTimer &timer, int n, TA *,
                 void (TA::*ff)(int, TB&, TC&, TD&, TE&, TF&, TG*), TB&, TC&, TD&, TE&, TF&, TG*);

// six reference arguments timed execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void timedParal6R(DistTimer &timer, int n, TA *,
                 void (TA::*ff)(int, TB&, TC&, TD&, TE&, TF&, TG&), TB&, TC&, TD&, TE&, TF&, TG&);

// five pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TF>
void execParal(int n, TA *, void (TA::*f)(int, TB, TC, TD, TE, TF), TB, TC, TD, TE, TF);

// five pointer arguments timed execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TF>
void timedParal(DistTimer&, int n, TA *, void (TA::*f)(int, TB, TC, TD, TE, TF), TB, TC, TD, TE, TF);

// six pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void execParal(int n,  TA *, void (TA::*f)(int, TB, TC, TD, TE, TF, TG), TB, TC, TD, TE, TF, TG);

// six pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TF, class TG>
void timedParal(DistTimer &timer, int n,  TA *, void (TA::*f)(int, TB, TC, TD, TE, TF, TG), TB, TC, TD, TE, TF, TG);

// seven pointer arguments execute in parallel
template <class TA, class TB, class TC, class TD, class TE, class TF, class TG,
          class TH>
void execParal(int n,  TA *, void (TA::*f)(int, TB, TC, TD, TE, TF, TG,TH),
               TB, TC, TD, TE, TF, TG, TH);

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

/*
// PJSA 2-25-03 new helper functions to do parallel sum, min, max etc

template <class Op, class A, class B, class Arg1, class PArg1>
B paralApplyOp<Op>(int n, A **a, B(A::*fct)(T2), PArg1 arg1);
*/


#ifdef _TEMPLATE_FIX_
#include <Threads.d/PHelper.C>
#endif

#endif
