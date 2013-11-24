#ifndef _SENSITIVITYINFO_H_
#define _SENSITIVITYINFO_H_

#include <cstdio>

struct SensitivityInfo {

   enum Type { WeightWRTthickness, StressVMWRTthickness };

   Type type;
   int numParam; // For example, numParam specifies the number of thickness groups
};

#endif
