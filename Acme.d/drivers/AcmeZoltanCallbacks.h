// $Id: AcmeZoltanCallbacks.h,v 2002.2 2004/06/24 14:51:35 mwglass Exp $

#ifndef CONTACT_NO_MPI

#include "Contact_Defines.h"
#include "lbi_const.h"

extern "C" {
  ZOLTAN_NUM_OBJ_FN            AcmeQueryNumObjects;
  ZOLTAN_OBJ_LIST_FN           AcmeQueryObjectList;
  ZOLTAN_NUM_GEOM_FN           AcmeQueryNumGeomObjects;
  ZOLTAN_GEOM_MULTI_FN         AcmeQueryGeomMultiValues;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     AcmeMigrateEntitySizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     AcmeMigrateEntityMultiPack;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   AcmeMigrateEntityMultiUnpack;
}

#endif
