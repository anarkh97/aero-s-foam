// This file contains template versions of commonly 
// used math functions

  template <class T>
  inline T myMax(T a, T b) { return a > b ? a : b; }

  template <class T>
  inline T myMin(T a, T b) { return a < b ? a : b; }

  template <class T>
  inline T abs(T x) { return x > 0 ? x : -x;}

  template <class T>
  inline void
  mySwap(T &a, T &b)
  {
    T  t = b;
       b = a;
       a = t;
  }

