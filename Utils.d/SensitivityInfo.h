#ifndef _SENSITIVITYINFO_H_
#define _SENSITIVITYINFO_H_

#include <cstdio>

struct SensitivityInfo {

   enum Type { WeightWRTthickness, WeightWRTshape,
               StressVMWRTthickness, StressVMWRTshape,
               StressVMWRTmach, StressVMWRTalpha, StressVMWRTbeta,
               DisplacementWRTthickness,
               DisplacementWRTshape,
               DisplacementWRTmach,
               DisplacementWRTalpha,
               DisplacementWRTbeta,
 };

   Type type;
   int surface;  // surface where shell type sensitivity is evaluated
};

#endif
