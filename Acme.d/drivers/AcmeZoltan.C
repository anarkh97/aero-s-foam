// $Id: AcmeZoltan.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef CONTACT_NO_MPI

#include <iostream.h>
#include "DriverTopology.h"
#include "AcmeZoltan.h"
#include "AcmeZoltanCallbacks.h"
#include "contact_assert.h"

AcmeZoltan::AcmeZoltan(MPI_Comm communicator, int& error )
{
  int tmpReturn = Zoltan_Initialize( 0, NULL, &version );
  Zoltan_Ptr_ = Zoltan_Create( communicator );
  POSTCONDITION( Zoltan_Ptr_ );
  if( Zoltan_Ptr_ )
    error = ZOLTAN_OK;
  else
    error = ZOLTAN_FATAL;
 num_gid_entries   = 0;
 num_lid_entries   = 0;
 num_import        = 0;
 import_global_ids = NULL;
 import_local_ids  = NULL;
 import_procs      = NULL;
 num_export        = 0;
 export_global_ids = NULL;
 export_local_ids  = NULL;
 export_procs      = NULL;
}

AcmeZoltan::~AcmeZoltan()
{
  Zoltan_Destroy( &Zoltan_Ptr_ );
}

//Load Balance Calls

void AcmeZoltan::Set_GeomCallBacks ( DriverTopology* topology)
{
  Zoltan_Set_Num_Obj_Fn( Zoltan_Ptr_, 
			 AcmeQueryNumObjects, 
			 REINTERPRET_CAST(void*) (topology) );
  
  Zoltan_Set_Obj_List_Fn( Zoltan_Ptr_, 
			  AcmeQueryObjectList, 
			  REINTERPRET_CAST(void*) (topology) );
  
  Zoltan_Set_Num_Geom_Fn( Zoltan_Ptr_, 
			  AcmeQueryNumGeomObjects, 
			  REINTERPRET_CAST(void*) (topology) );
  
  Zoltan_Set_Geom_Multi_Fn( Zoltan_Ptr_, 
                            AcmeQueryGeomMultiValues, 
                            REINTERPRET_CAST(void*) (topology) );
}

void AcmeZoltan::Set_DLBmultiCallBacks( DriverTopology* topology )
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_, 
                                AcmeMigrateEntitySizes, 
                                REINTERPRET_CAST(void*) (topology) );

  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_, 
				AcmeMigrateEntityMultiPack,
				REINTERPRET_CAST(void*) (topology) );

  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_, 
				  AcmeMigrateEntityMultiUnpack,
				  REINTERPRET_CAST(void*) (topology) );
}

int AcmeZoltan::Set_Method( char * string )
{
  char* method = (char*) "LB_METHOD";
  return Zoltan_Set_Param( Zoltan_Ptr_, method, string );
}

int AcmeZoltan::Set_Param( char * param, char * val )
{
  return Zoltan_Set_Param( Zoltan_Ptr_, param, val );
}

int AcmeZoltan::Balance()
{
  //Zoltan_Generate_Files(Zoltan_Ptr_, "heaphy", 1);
  return Zoltan_LB_Balance( Zoltan_Ptr_, &new_decomp, 
			    &num_gid_entries, &num_lid_entries,
			    &num_import, &import_global_ids, 
			    &import_local_ids, &import_procs, 
			    &num_export, &export_global_ids,
			    &export_local_ids, &export_procs );
}

void AcmeZoltan::Evaluate ( int    print_stats,
			       int*   num_objects,
			       float* object_weights,
                               int*   num_cuts,
			       float* cut_weights,
			       int*   num_boundary_objects,
			       int*   num_adj_procs,
			       int*   ierr )
{
  *ierr = Zoltan_LB_Eval( Zoltan_Ptr_, print_stats, num_objects, 
			  object_weights, num_cuts, cut_weights, 
			  num_boundary_objects, num_adj_procs );
}

int AcmeZoltan::Free_Data ()
{
  return Zoltan_LB_Free_Data( &import_global_ids, &import_local_ids, 
			      &import_procs, &export_global_ids, 
			      &export_local_ids, &export_procs );
}

//Decomposition Augmentation

int AcmeZoltan::Point_Assign ( Real* coords, int* proc )
{
  PRECONDITION( sizeof(Real) == sizeof(double) );
  return Zoltan_LB_Point_Assign( Zoltan_Ptr_, coords, proc );
}

int AcmeZoltan::Box_Assign ( Real xmin,
			     Real ymin,
			     Real zmin,
			     Real xmax,
			     Real ymax,
			     Real zmax,
			     int* procs,
			     int* numprocs )
{
  PRECONDITION( sizeof(Real) == sizeof(double) );
  return Zoltan_LB_Box_Assign( Zoltan_Ptr_, xmin, ymin, zmin, xmax, 
			       ymax, zmax, procs, numprocs );
}

int AcmeZoltan::Compute_Destinations ( int            Num_import,
				       ZOLTAN_ID_PTR  Import_global_ids,
				       ZOLTAN_ID_PTR  Import_local_ids,
				       int*           Import_procs,
				       int*           Num_export,
				       ZOLTAN_ID_PTR* Export_global_ids,
				       ZOLTAN_ID_PTR* Export_local_ids,
				       int**          Export_procs )
{
  return Zoltan_Compute_Destinations( Zoltan_Ptr_, 
				      Num_import,       Import_global_ids, 
				      Import_local_ids, Import_procs, 
				      Num_export,       Export_global_ids, 
				      Export_local_ids, Export_procs );
}

int AcmeZoltan::Help_Migrate ( int           Num_import,
			       ZOLTAN_ID_PTR Import_global_ids,
			       ZOLTAN_ID_PTR Import_local_ids,
			       int*          Import_procs,
			       int           Num_export,
			       ZOLTAN_ID_PTR Export_global_ids,
			       ZOLTAN_ID_PTR Export_local_ids,
			       int*          Export_procs )
{
  return Zoltan_Help_Migrate( Zoltan_Ptr_, 
			      Num_import, Import_global_ids, 
			      Import_local_ids, Import_procs,
			      Num_export, Export_global_ids, 
			      Export_local_ids, Export_procs );
}


#endif
