// $Id$

#include "Contact_Defines.h"
#include "lbi_const.h"

extern "C" {
  ZOLTAN_NUM_OBJ_FN            ContactQueryNumObjects;
  ZOLTAN_OBJ_LIST_FN           ContactQueryObjectList;
  ZOLTAN_NUM_GEOM_FN           ContactQueryNumGeomObjects;
  ZOLTAN_GEOM_MULTI_FN         ContactQueryGeomMultiValues;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactMigrateEntityExportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactMigrateEntityExportPack;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactMigrateEntityUnpack;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactMigrateInteractionSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactMigratePackInteractions;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactMigrateUnpackInteractions;
  
  ZOLTAN_OBJ_SIZE_FN           ContactHostidQuerySize;
  ZOLTAN_PACK_OBJ_FN           ContactHostidQueryPack;
  ZOLTAN_UNPACK_OBJ_FN         ContactHostidQueryUnpack;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactDynamicLoadBalanceSize;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactDynamicLoadBalancePack;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactDynamicLoadBalanceUnpack;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactGhostingExportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactGhostingExportPack;
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactGhostingImportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactGhostingImportPack;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactGhostingUnpack;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactUpdateGhostingExportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactUpdateGhostingExportPack;
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactUpdateGhostingImportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactUpdateGhostingImportPack;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactUpdateGhostingUnpack;
  
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactUpdateTiedImportSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactUpdateTiedImportPack;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactUpdateTiedUnpack;

}

