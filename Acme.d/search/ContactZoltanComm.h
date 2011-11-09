// $Id$

#ifndef _ContactZoltanComm_h_
#define _ContactZoltanComm_h_

#ifndef CONTACT_NO_MPI

#include "lbi_const.h"
#include "Contact_Defines.h"
#include "ContactZoltanID.h"
#include <vector>
#include <map>

class ContactParOStream;

class ContactZoltanComm {

 public:
 
  enum Zoltan_Dir { ZOLTAN_EXPORT, ZOLTAN_IMPORT, ZOLTAN_UNKNOWN };

  ContactZoltanComm();
  ContactZoltanComm( Zoltan_Dir);
  ~ContactZoltanComm();

  void CleanUp();
  void Initialize(Zoltan_Dir);

  void Add_Import( LB_ID_PTR lid, LB_ID_PTR gid, int proc, int check=1 );
  void Add_Export( LB_ID_PTR lid, LB_ID_PTR gid, int proc, int check=1 );
  
  inline int       Num_Import() {return num_import;};
  inline int       Num_Export() {return num_export;};
  inline LB_ID_PTR Import_GIDS() {return import_gids;};
  inline LB_ID_PTR Export_GIDS() {return export_gids;};
  inline LB_ID_PTR Import_LIDS() {return import_lids;};
  inline LB_ID_PTR Export_LIDS() {return export_lids;};
  inline int*      Import_Procs() {return import_procs;};
  inline int*      Export_Procs() {return export_procs;};
  
  void Set_Export(int num, LB_ID_PTR gids, LB_ID_PTR lids, int* procs);

  
  inline Zoltan_Dir Direction() {return direction;};
  
  void Print(ContactParOStream& os);

 private:
  Zoltan_Dir direction;

  int        num_import;
  int        cur_import_capacity;
  LB_ID_PTR  import_gids;
  LB_ID_PTR  import_lids;
  int*       import_procs;

  int        num_export;
  int        cur_export_capacity;
  LB_ID_PTR  export_gids;
  LB_ID_PTR  export_lids;
  int*       export_procs;


  std::map<std::pair<ContactZoltanGID, int>, int> lookup;
  void add_hash_entry( LB_ID_PTR, int );
  int  find_hash_entry( LB_ID_PTR, int );
};


#endif
#endif
