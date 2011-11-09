// $Id$

#ifndef ContactEnfZoltan_h_
#define ContactEnfZoltan_h_

#ifndef CONTACT_NO_MPI

#include "Contact_Defines.h"
#include "lbi_const.h"

// External call back function to size/pack/unpack objects

extern "C" {

  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactEnfMigrateNodeSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactEnfMigratePackNodes;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactEnfMigrateUnpackNodes;

  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactEnfMigrateSharedNodeSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactEnfMigratePackSharedNodes;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactEnfMigrateUnpackSharedNodes;
 
  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactEnfMigrateFaceSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactEnfMigratePackFaces;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactEnfMigrateUnpackFaces;

  ZOLTAN_OBJ_SIZE_MULTI_FN     ContactEnfMigrateElementSizes;
  ZOLTAN_PACK_OBJ_MULTI_FN     ContactEnfMigratePackElements;
  ZOLTAN_UNPACK_OBJ_MULTI_FN   ContactEnfMigrateUnpackElements;

}

class ContactEnforcement;

class ContactEnfZoltan {

 public:
  
  ContactEnfZoltan( MPI_Comm communicator, Zoltan_Struct*, int& ierr );
  ~ContactEnfZoltan();

  void Set_MultiNodeCallBacks( ContactEnforcement* );
  void Set_MultiSharedNodeCallBacks( ContactEnforcement* );
  void Set_MultiFaceCallBacks( ContactEnforcement* );
  void Set_MultiElementCallBacks( ContactEnforcement* );

  int Compute_Destinations ( int            num_import,
                             ZOLTAN_ID_PTR  import_global_ids,
                             ZOLTAN_ID_PTR  import_local_ids,
                             int*           import_procs,
                             int*           num_export,
                             ZOLTAN_ID_PTR* export_global_ids,
                             ZOLTAN_ID_PTR* export_local_ids,
                             int**          export_procs );
 
  int Help_Migrate ( int           num_import,
                     ZOLTAN_ID_PTR import_global_ids,
                     ZOLTAN_ID_PTR import_local_ids,
                     int*          import_procs,
                     int           num_export,
                     ZOLTAN_ID_PTR export_global_ids,
                     ZOLTAN_ID_PTR export_local_ids,
                     int*          export_procs );
                                                                                
 private:

  Zoltan_Struct* Zoltan_Ptr_;
  int            num_import;
  ZOLTAN_ID_PTR  import_global_ids;
  ZOLTAN_ID_PTR  import_local_ids;
  int*           import_procs;
  int            num_export;
  ZOLTAN_ID_PTR  export_global_ids;
  ZOLTAN_ID_PTR  export_local_ids;
  int*           export_procs;

};

#endif // CONTACT_MPI
#endif // ContactEnfZoltan_h_
