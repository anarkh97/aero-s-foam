#ifndef _GLOBALINT_H_
#define _GLOBALINT_H_

#define USE_INT32

#ifdef USE_INT32
typedef int GlobalInt;
#else
typedef long long GlobalInt;
#endif

#endif
