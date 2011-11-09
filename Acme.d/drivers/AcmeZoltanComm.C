// $Id: AcmeZoltanComm.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeZoltanComm.h"
#include "AcmeZoltanID.h"
#include "params.h"
#include "contact_assert.h"
#include <string.h>

#ifndef CONTACT_NO_MPI
	   
AcmeZoltanComm::AcmeZoltanComm( Zoltan_Dir dir, int number )
{
  direction = dir;
  max_size  = number<=0?100:number;
  switch (dir) {
  case ZOLTAN_IMPORT:
    import_gids  = new LB_ID_TYPE [max_size*ACME_ZOLTAN_GID_SIZE];
    import_lids  = new LB_ID_TYPE [max_size*ACME_ZOLTAN_LID_SIZE];
    import_procs = new int [max_size];
    export_gids  = NULL;
    export_lids  = NULL;
    export_procs = NULL;
    num_import   = 0;
    num_export   = 0;
    break;
  case ZOLTAN_EXPORT:
    import_gids  = NULL;
    import_lids  = NULL;
    import_procs = NULL;
    export_gids  = new LB_ID_TYPE [max_size*ACME_ZOLTAN_GID_SIZE];
    export_lids  = new LB_ID_TYPE [max_size*ACME_ZOLTAN_LID_SIZE];
    export_procs = new int [max_size];
    num_import   = 0;
    num_export   = 0;
    break;
  }
  create_hash();
}


AcmeZoltanComm::~AcmeZoltanComm()
{
  switch (direction) {
  case ZOLTAN_IMPORT:
    delete [] import_gids;
    delete [] import_lids;
    delete [] import_procs;
    LB_Free_Data( NULL, NULL, NULL,
                  &export_gids, 
                  &export_lids, 
                  &export_procs );
    break;
  case ZOLTAN_EXPORT:
    delete [] export_gids;
    delete [] export_lids;
    delete [] export_procs;
    LB_Free_Data( &import_gids, 
                  &import_lids, 
                  &import_procs, 
                  NULL, NULL, NULL );
    break;
  }
  delete_hash();
}

void AcmeZoltanComm::Add_Import(LB_ID_PTR lid, LB_ID_PTR gid, int proc)
{
  PRECONDITION(num_import < max_size);
  
  int i;
  // If the node is in the list, just return

  if (find_hash_entry(lid, gid, proc)) return;
  add_hash_entry(lid, gid, proc);
  int nbytes_lid = ACME_ZOLTAN_LID_SIZE*sizeof(LB_ID_TYPE);
  int nbytes_gid = ACME_ZOLTAN_GID_SIZE*sizeof(LB_ID_TYPE);
  memcpy(&import_lids[num_import*ACME_ZOLTAN_LID_SIZE],lid,nbytes_lid);
  memcpy(&import_gids[num_import*ACME_ZOLTAN_GID_SIZE],gid,nbytes_gid);
  import_procs[num_import] = proc;
  num_import++;
  if (num_import == max_size) {
    max_size += 100;
    LB_ID_PTR new_import_gids  = new LB_ID_TYPE [max_size*ACME_ZOLTAN_GID_SIZE];
    LB_ID_PTR new_import_lids  = new LB_ID_TYPE [max_size*ACME_ZOLTAN_LID_SIZE];
    int*      new_import_procs = new int [max_size];
    int size = num_import*sizeof(LB_ID_TYPE);
    memcpy(new_import_gids,import_gids,ACME_ZOLTAN_GID_SIZE*size);
    memcpy(new_import_lids,import_lids,ACME_ZOLTAN_LID_SIZE*size);
    memcpy(new_import_procs,import_procs,num_import*sizeof(int));
    delete [] import_gids;
    delete [] import_lids;
    delete [] import_procs;
    import_gids  = new_import_gids;
    import_lids  = new_import_lids;
    import_procs = new_import_procs;
    delete_hash();
    create_hash();
    LB_ID_PTR lid_ptr = import_lids;
    LB_ID_PTR gid_ptr = import_gids;
    for(i=0 ; i<num_import ; i++ ) {
      add_hash_entry(lid_ptr, gid_ptr, import_procs[i]);
      lid_ptr += ACME_ZOLTAN_LID_SIZE;
      gid_ptr += ACME_ZOLTAN_GID_SIZE;
    }
  }
}

void AcmeZoltanComm::Add_Export(LB_ID_PTR lid, LB_ID_PTR gid, 
                                   int proc, int check)
{
  PRECONDITION(num_export < max_size);
  
  int i;
  
  if (check) {
    if (find_hash_entry(lid, gid, proc)) return;
    add_hash_entry(lid, gid, proc);
  }
  int nbytes_lid = ACME_ZOLTAN_LID_SIZE*sizeof(LB_ID_TYPE);
  int nbytes_gid = ACME_ZOLTAN_GID_SIZE*sizeof(LB_ID_TYPE);
  memcpy(&export_lids[num_export*ACME_ZOLTAN_LID_SIZE],lid,nbytes_lid);
  memcpy(&export_gids[num_export*ACME_ZOLTAN_GID_SIZE],gid,nbytes_gid);
  export_procs[num_export] = proc;
  num_export++;
  if (num_export == max_size) {
    max_size += 100;
    LB_ID_PTR new_export_gids  = new LB_ID_TYPE [max_size*ACME_ZOLTAN_GID_SIZE];
    LB_ID_PTR new_export_lids  = new LB_ID_TYPE [max_size*ACME_ZOLTAN_LID_SIZE];
    int*      new_export_procs = new int [max_size];
    int size = num_export*sizeof(LB_ID_TYPE);
    memcpy(new_export_gids,export_gids,ACME_ZOLTAN_GID_SIZE*size);
    memcpy(new_export_lids,export_lids,ACME_ZOLTAN_LID_SIZE*size);
    memcpy(new_export_procs,export_procs,num_export*sizeof(int));
    delete [] export_gids;
    delete [] export_lids;
    delete [] export_procs;
    export_gids  = new_export_gids;
    export_lids  = new_export_lids;
    export_procs = new_export_procs;
    if (check) {
      delete_hash();
      create_hash();
      LB_ID_PTR lid_ptr = export_lids;
      LB_ID_PTR gid_ptr = export_gids;
      for(i=0 ; i<num_export ; i++ ) {
        add_hash_entry(lid_ptr, gid_ptr, export_procs[i]);
        lid_ptr += ACME_ZOLTAN_LID_SIZE;
        gid_ptr += ACME_ZOLTAN_GID_SIZE;
      }
    }
  }
}

void AcmeZoltanComm::create_hash()
{
  int i;

  hash_space = NULL;
  hash_space_list = NULL;

#define BIN_FRACTION 0.25
#define BIN_MINIMUM 100

  // Create space for bins
  nbins = (int) (BIN_FRACTION*max_size);
  nbins = MAX( nbins, BIN_MINIMUM );

#undef BIN_FRACTION
#undef BIN_MINIMUM

  // Round up nbins to make it more prime-like
  if ( !(nbins % 2) ) nbins += 1;
  if ( !(nbins % 3) ) nbins += 2;
  if ( !(nbins % 6) ) nbins += 6;

  // Allocate the bins
  bins = new hash*[nbins];

  // Initialize the bins to empty
  for( i=0 ; i<nbins ; i++ )
    bins[i] = NULL;

  // Create hash space
  int n = max_size;
  if( n<=0 ) n=1;

  hash_space = new hash[n];

  hash_list* new_list = new hash_list;
  new_list->next = hash_space_list;
  new_list->ptr = hash_space;
  hash_space_list = new_list;
}

void AcmeZoltanComm::delete_hash()
{
  hash_list* next_list;
  while( hash_space_list != NULL ){
    delete [] hash_space_list->ptr;
    next_list = hash_space_list->next;
    delete hash_space_list;
    hash_space_list = next_list;
  }
  delete [] bins;
}

void AcmeZoltanComm::add_hash_entry( LB_ID_PTR lid, LB_ID_PTR gid, int proc )
{
  hash* newpt;
  hash* ptr;
  hash** previous;
  int ibin;
  AcmeZoltanLID zoltan_lid(lid);
  AcmeZoltanGID zoltan_gid(gid);

  ibin = zoltan_gid % nbins;
  previous = bins + ibin;
  ptr = bins[ibin];
  
  while( ptr != NULL && compare_entry(ptr, zoltan_gid, proc)<0){
    previous = &(ptr->next);
    ptr = ptr->next;
  }
  
  if( ptr == NULL || compare_entry(ptr, zoltan_gid, proc)>0) {
    newpt = hash_space++;
    *previous = newpt;
    newpt->proc       = proc;
    newpt->zoltan_lid = zoltan_lid;
    newpt->zoltan_gid = zoltan_gid;
    newpt->next = ptr;
  }
}

int AcmeZoltanComm::find_hash_entry( LB_ID_PTR lid, LB_ID_PTR gid, int proc )
{
  hash* ptr;
  int ibin;
  AcmeZoltanGID zoltan_gid(gid);

  ibin = zoltan_gid % nbins;
  ptr = bins[ibin];
    
  while( ptr != NULL && compare_entry(ptr, zoltan_gid, proc)<0){
    ptr = ptr->next;
  }
  
  if( ptr == NULL || compare_entry(ptr, zoltan_gid, proc)>0) return 0;
  
  return 1;
}

int  AcmeZoltanComm::compare_entry(hash* ptr,              
                                   AcmeZoltanGID zoltan_gid, 
                                   int proc)
{
  if (direction==ZOLTAN_EXPORT) {
    if (proc<ptr->proc) return -1;
    else if (proc>ptr->proc) return 1;
  }
  if (zoltan_gid<ptr->zoltan_gid) return -1;
  else if (zoltan_gid>ptr->zoltan_gid) return 1;
  return 0;
}
#endif

