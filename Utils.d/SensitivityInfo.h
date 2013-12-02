#ifndef _SENSITIVITYINFO_H_
#define _SENSITIVITYINFO_H_

#include <cstdio>

struct SensitivityInfo {

   enum Type { WeightWRTthickness, StressVMWRTthickness, StressVMWRTdisplacement };

   Type type;
   int numParam; // For example, numParam specifies the number of thickness groups
                 // when displacement is used as parameter, numParam is equal to Domain::numdof()
};

#endif
