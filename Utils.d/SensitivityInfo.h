#ifndef _SENSITIVITYINFO_H_
#define _SENSITIVITYINFO_H_

#include <cstdio>

struct SensitivityInfo {

   enum Type { WeightWRTthickness, StressVMWRTthickness, StressVMWRTdisplacement };
   enum Method { Analytic, AutomaticDifferentiation, FiniteDifference };

   Type type;
   Method method;
   int numParam; // For example, numParam specifies the number of thickness groups
                 // when displacement is used as parameter, numParam is equal to Domain::numdof()
};

#endif
