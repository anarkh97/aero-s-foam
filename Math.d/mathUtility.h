#ifndef _MATH_UTILITY_H_
#define _MATH_UTILITY_H_

// Template function declarations for max, min, and abs
template <class T> inline T myMax(T a, T b);
template <class T> inline T myMin(T a, T b);
template <class T> inline T abs(T a);

template <class T> inline void mySwap(T &a, T &b);

template <class T> inline double delta(T a, T b) { return (a==b) ? 1.0 : 0.0; }

#ifdef _TEMPLATE_FIX_
#include <Math.d/mathUtility.C>
#endif

#endif
