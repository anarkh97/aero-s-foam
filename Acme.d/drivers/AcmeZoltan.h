// $Id: AcmeZoltan.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeZoltan_h_
#define _AcmeZoltan_h_

#ifndef CONTACT_NO_MPI

#include <lbi_const.h>
#include "Contact_Defines.h"

class DriverTopology;

class AcmeZoltan
{

public:

 AcmeZoltan( MPI_Comm communicator, int& ierr );

 ~AcmeZoltan();

 //Load Balance Calls

 void Set_GeomCallBacks  ( DriverTopology* );

 void Set_DLBmultiCallBacks( DriverTopology* );

 int Set_Method( char* string );

 int Set_Param( char* param, char * val );

 int Balance();

 void Evaluate( int    print_stats,
		int*   num_objects,
		float* object_weights,
                int*   num_cuts,
    		float* cut_weights,
		int*   num_boundary_objects,
		int*   num_adj_procs,
		int*   ierr );

 int Free_Data();
  
 //Decomposition Augmentation
 
 int Point_Assign ( Real* coords,
		    int* proc );

 int Box_Assign ( Real xmin,
		  Real ymin,
		  Real zmin,
		  Real xmax,
		  Real ymax,
		  Real zmax,
		  int* procs,
		  int* numprocs );
                  
 //Migration Functions
 
  int Compute_Destinations ( int        num_import,
			     ZOLTAN_ID_PTR  import_global_ids,
			     ZOLTAN_ID_PTR  import_local_ids,
			     int*       import_procs,
			     int*       num_export,
			     ZOLTAN_ID_PTR* export_global_ids,
			     ZOLTAN_ID_PTR* export_local_ids,
			     int**      export_procs );

  int Help_Migrate ( int       num_import,
		     ZOLTAN_ID_PTR import_global_ids,
		     ZOLTAN_ID_PTR import_local_ids,
		     int*      import_procs,
		     int       num_export,
		     ZOLTAN_ID_PTR export_global_ids,
		     ZOLTAN_ID_PTR export_local_ids,
		     int*      export_procs );
  
  inline float     Version() { return version; };
  inline int       New_Decomp() { return new_decomp; };
  inline int       Num_Import() { return num_import; };
  inline ZOLTAN_ID_PTR Import_Global_IDs() { return import_global_ids; };
  inline ZOLTAN_ID_PTR Import_Local_IDs() { return import_local_ids; };
  inline int*      Import_Procs() { return import_procs; };
  inline int       Num_Export() { return num_export; };
  inline ZOLTAN_ID_PTR Export_Global_IDs() { return export_global_ids; };
  inline ZOLTAN_ID_PTR Export_Local_IDs() { return export_local_ids; };
  inline int*      Export_Procs() { return export_procs; };
  inline void      Position(VariableHandle p) { position=p; };
  inline VariableHandle Position() { return position; };
  inline Zoltan_Struct* Get_Data_Structure() { return Zoltan_Ptr_; };

private:

 VariableHandle position;
 
 Zoltan_Struct* Zoltan_Ptr_; 
 float          version;
 int            new_decomp;
 int            num_gid_entries;
 int            num_lid_entries;
 int            num_import;
 ZOLTAN_ID_PTR  import_global_ids;
 ZOLTAN_ID_PTR  import_local_ids;
 int*           import_procs;
 int            num_export;
 ZOLTAN_ID_PTR  export_global_ids;
 ZOLTAN_ID_PTR  export_local_ids;
 int*           export_procs;

};

#endif
#endif
