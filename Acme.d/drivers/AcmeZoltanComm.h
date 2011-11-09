// $Id: AcmeZoltanComm.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeZoltanComm_h_
#define _AcmeZoltanComm_h_

#ifndef CONTACT_NO_MPI

#include "lbi_const.h"
#include "Contact_Defines.h"
#include "AcmeZoltanID.h"

class AcmeZoltanComm {

 public:
 
  enum Zoltan_Dir { ZOLTAN_EXPORT, ZOLTAN_IMPORT };

  AcmeZoltanComm( Zoltan_Dir, int );
  ~AcmeZoltanComm();

  void Add_Import( LB_ID_PTR lid, LB_ID_PTR gid, int proc );
  void Add_Export( LB_ID_PTR lid, LB_ID_PTR gid, int proc, int check=1 );
  
  inline int       Num_Import() {return num_import;};
  inline int       Num_Export() {return num_export;};
  inline LB_ID_PTR Import_GIDS() {return import_gids;};
  inline LB_ID_PTR Export_GIDS() {return export_gids;};
  inline LB_ID_PTR Import_LIDS() {return import_lids;};
  inline LB_ID_PTR Export_LIDS() {return export_lids;};
  inline int*      Import_Procs() {return import_procs;};
  inline int*      Export_Procs() {return export_procs;};
  
  inline void Num_Import(int n) {num_import=n;};
  inline void Num_Export(int n) {num_export=n;};
  inline void Import_GIDS(LB_ID_PTR n) {import_gids=n;};
  inline void Export_GIDS(LB_ID_PTR n) {export_gids=n;};
  inline void Import_LIDS(LB_ID_PTR n) {import_lids=n;};
  inline void Export_LIDS(LB_ID_PTR n) {export_lids=n;};
  inline void Import_Procs(int* n) {import_procs=n;};
  inline void Export_Procs(int* n) {export_procs=n;};
  
  inline Zoltan_Dir Direction(Zoltan_Dir n) {return direction=n;};
  
  struct hash {
    int proc;
    AcmeZoltanLID zoltan_lid;
    AcmeZoltanGID zoltan_gid;
    struct hash *next;
  };

  struct hash_list {
    struct hash* ptr;
    struct hash_list *next;
  };

 private:
  Zoltan_Dir direction;
  int        max_size;
  int        num_import;
  LB_ID_PTR  import_gids;
  LB_ID_PTR  import_lids;
  int*       import_procs;
  int        num_export;
  LB_ID_PTR  export_gids;
  LB_ID_PTR  export_lids;
  int*       export_procs;

  int nbins;
  hash* hash_space;
  hash** bins;
  hash_list* hash_space_list;

  void create_hash( );
  void delete_hash( );
  void add_hash_entry( LB_ID_PTR, LB_ID_PTR, int );
  int  find_hash_entry( LB_ID_PTR, LB_ID_PTR, int );
  int  compare_entry( hash*, AcmeZoltanGID, int );

};
#endif

#endif
