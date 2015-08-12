#ifndef _SENSITIVITYINFO_H_
#define _SENSITIVITYINFO_H_

#include <cstdio>

struct SensitivityInfo {

   enum Method { Direct, Adjoint };

   enum Type { WeightWRTthickness, WeightWRTshape,
               StressVMWRTthickness, StressVMWRTshape,
               StressVMWRTmach, StressVMWRTalpha, StressVMWRTbeta,
               DisplacementWRTthickness,
               DisplacementWRTshape,
               DisplacementWRTmach,
               DisplacementWRTalpha,
               DisplacementWRTbeta,
 };

   Method method;
   Type type;
   int surface;  // surface where shell type sensitivity is evaluated

   SensitivityInfo() {
     method = SensitivityInfo::Direct;
   }

};

#endif
