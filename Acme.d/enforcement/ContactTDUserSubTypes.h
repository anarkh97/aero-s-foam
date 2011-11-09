
#ifndef ContactTDUserSubTypes_h_
#define ContactTDUserSubTypes_h_

typedef int CONTACT_INIT_MODEL_FN(
  void *enf,     // enforcement object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &istat    // return code
);

typedef int CONTACT_INIT_TIME_STEP_FN(
  void *enf,     // enforcement object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &istat    // return code
);

typedef void CONTACT_INIT_NODE_STATE_DATA_FN(
  void *enf,     // enforcement object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  Real *sdata,   // node state data
  int  &istat    // return code
);

typedef void CONTACT_NUM_NODE_STATE_VARS_FN(
  void *enf,     // enforcement object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &istat    // return code
);

typedef int CONTACT_LIMIT_FORCE_FN(
  void *enf,     // enforcement object
  void *ni,      // interaction object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &type,    // type of interaction, nfi or nsi
  int  &index,   // node index
  Real &gap, 
  Real *rel_disp, 
  Real *slip, 
  Real *normal,
  Real &dt, 
  Real &area, 
  Real *force,
  int  &istat    // return code
);

typedef int CONTACT_INTERACTION_ACTIVE_FN(
  void *enf,     // enforcement object
  void *ni,      // interaction object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &type,    // type of interaction, nfi or nsi
  int  &index,   // node index
  Real &gap,     // gap
  int  &istat    // return code
);

typedef int CONTACT_INTERACTION_TYPE_FN(
  void *enf,     // enforcement object
  void *ni,      // interaction object
  int  &id,      // enforcement model id
  Real *rdata,   // user define real data
  int  *idata,   // user define integer data
  int  &type,    // type of interaction, nfi or nsi
  int  &index,   // node index
  int  &istat    // return code
);

#endif
