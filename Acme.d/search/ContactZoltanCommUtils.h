// $Id$

#ifndef ContactZoltanCommUtils_h_
#define ContactZoltanCommUtils_h_

#ifndef CONTACT_NO_MPI

#include <mpi.h>
#include <lbi_const.h>
#include <zz_const.h>
#include "Contact_Defines.h"
#include "contact_assert.h"
#include "ContactZoltanComm.h"
#include "ContactParOStream.h"

class ContactZoltanCommUtils {

 public:
 
    enum CommDirection { IMPORT, EXPORT, BOTH };
 
    ContactZoltanCommUtils(CommDirection dir,
                           Zoltan_Struct* z_ptr,
                           int num, LB_ID_TYPE* gids, 
                           LB_ID_TYPE* lids, int* pids,
                           ContactParOStream& os);
                           
    ContactZoltanCommUtils(ContactZoltanComm* comm1,
                           ContactZoltanComm* comm2,
                           Zoltan_Struct*     z_ptr,
                           ContactParOStream& os);

    ~ContactZoltanCommUtils();
    
    void Resize(ContactParOStream& os);
    
    void Migrate(ContactParOStream& os);
    
    inline int       Num_Import()   {return num_import;};
    inline int       Num_Export()   {return num_export;};
    inline LB_ID_PTR Import_GIDS()  {return import_gids;};
    inline LB_ID_PTR Export_GIDS()  {return export_gids;};
    inline LB_ID_PTR Import_LIDS()  {return import_lids;};
    inline LB_ID_PTR Export_LIDS()  {return export_lids;};
    inline int*      Import_Procs() {return import_pids;};
    inline int*      Export_Procs() {return export_pids;};
  
  private:
    int              type;
    
    CommDirection    direction;
  
    Zoltan_Struct*   zz;
    
    ZOLTAN_COMM_OBJ* import_plan;
    ZOLTAN_COMM_OBJ* export_plan;
    
    int              num_import;
    LB_ID_TYPE*      import_gids;
    LB_ID_TYPE*      import_lids;
    int*             import_pids;
    
    int              num_export;
    LB_ID_TYPE*      export_gids;
    LB_ID_TYPE*      export_lids;
    int*             export_pids;
    
    char*            send_buffer;
    char*            recv_buffer;
    
    int*             send_sizes;
    int*             recv_sizes;
    int*             send_indices;
    int*             recv_indices;

};

#endif

#endif
