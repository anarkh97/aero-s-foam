#ifndef _SENSITIVITYINFO_H_
#define _SENSITIVITYINFO_H_

#include <cstdio>

struct SensitivityInfo {

   enum Type { WeightWRTthickness, WeightWRTshape, StressVMWRTthickness, StressVMWRTdisplacement, 
               StressVMWRTshape, StiffnessWRTthickness, LinearStaticWRTthickness, LinearStaticWRTshape,  
               GravityWRTthickness };
   enum Method { Analytic, AutomaticDifferentiation, FiniteDifference };

   Type type;
   Method method;
   int surface;  // surface where shell type sensitivity is evaluated
   void initialize() {
     surface = 2;  // default is midsurface
   }
};

#endif
