#ifndef _SENSITIVITYINFO_H_
#define _SENSITIVITYINFO_H_

#include <cstdio>

struct SensitivityInfo {

   enum Type { WeightWRTthickness, StressVMWRTthickness, StressVMWRTdisplacement, StressVMWRTcoordinate, StiffnessWRTthickness, LinearStaticWRTthickness };
   enum Method { Analytic, AutomaticDifferentiation, FiniteDifference };

   Type type;
   Method method;
   int surface;  // surface where shell type sensitivity is evaluated
   int numParam; // For example, numParam specifies the number of thickness groups
                 // when displacement is used as parameter, numParam is equal to Domain::numdof()
   void initialize() {
     surface = 2;  // default is midsurface
     numParam = 0;
   }
};

#endif
