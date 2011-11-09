// $Id: AcmeZoltanCallbacks.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef CONTACT_NO_MPI

#include "lbi_const.h"
#include "AcmeZoltanCallbacks.h"
#include "DriverTopology.h"
#include "AcmeZoltan.h"
#include "AcmeZoltanID.h"
#include "AcmeZoltanComm.h"
#include "AcmeEntityHash.h"
#include "AcmeNode.h"
#include "AcmeEntityBlock.h"
#include "AcmeHexElemL8.h"
#include "AcmeQuadFaceL4.h"
#include "AcmeQuadFaceQ8.h"
#include "AcmeTriFaceL3.h"
#include "AcmeTriFaceQ6.h"
#include "params.h"
#include "Contact_Communication.h"
#include "contact_assert.h"
#include <iostream.h>
#include <stdio.h>
#include <string.h>

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      The following three routines are for the RCB decomposition.
//
//-------------------------------------------------------------------------

int AcmeQueryNumObjects(void *data, int *ierr)
{
  DriverTopology* topology = (DriverTopology *)data;
  int number_of_nodes = topology->NumNodes();
  AcmeNode** Nodes = 
    REINTERPRET_CAST(AcmeNode**)(topology->NodeList()->EntityList()); 
   
  int NumberOfLocalObjects  = 0;
  for (int i=0; i<number_of_nodes; i++) {
    if (Nodes[i]->Ownership() == AcmeEntity::OWNED &&
        Nodes[i]->State() == AcmeEntity::ACTIVE) {
      NumberOfLocalObjects++;
    }
  }
  *ierr = ZOLTAN_OK;
  return NumberOfLocalObjects;
}

void AcmeQueryObjectList(void *data, 
                         int num_gid_entries, int num_lid_entries, 
                         ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                         int wdim, float *objwgts, int *ierr)
{
  PRECONDITION(num_gid_entries==ACME_ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ACME_ZOLTAN_LID_SIZE);
  DriverTopology* topology = (DriverTopology *)data;
  int number_of_nodes = topology->NumNodes();
  AcmeNode** Nodes = 
    REINTERPRET_CAST(AcmeNode**)(topology->NodeList()->EntityList());  
  
  int k=0;
  PRECONDITION(num_gid_entries==ACME_ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ACME_ZOLTAN_LID_SIZE);
  for (int i=0; i<number_of_nodes; i++) {
    if (Nodes[i]->Ownership() == AcmeEntity::OWNED &&
        Nodes[i]->State() == AcmeEntity::ACTIVE) {
      Nodes[i]->ZoltanLID(&lids[k*ACME_ZOLTAN_LID_SIZE]);
      Nodes[i]->ZoltanGID(&gids[k*ACME_ZOLTAN_GID_SIZE]);
      if (wdim!=0) {
        for (int j=0; j<wdim; j++) {
          objwgts[k*wdim+j] = 1.0;
        }
      }
      k++;
    }
  }
  *ierr = ZOLTAN_OK;
}     

int AcmeQueryNumGeomObjects(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  DriverTopology* topology = (DriverTopology *)data;
  return topology->Dimensionality();
}

void AcmeQueryGeomMultiValues(void *data, 
                              int num_gid_entries, 
                              int num_lid_entries,
                              int num_obj, 
                              ZOLTAN_ID_PTR gids, 
                              ZOLTAN_ID_PTR lids,
                              int num_dim, 
                              double *geom_vec, 
                              int *ierr)
{
  PRECONDITION(num_gid_entries==ACME_ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ACME_ZOLTAN_LID_SIZE);
  DriverTopology* topology = (DriverTopology *)data;
  AcmeZoltan*   zoltan   = topology->Get_Zoltan(); 
  AcmeNode**    Nodes    = 
    REINTERPRET_CAST(AcmeNode**)(topology->NodeList()->EntityList()); 
  
  PRECONDITION(num_dim==topology->Dimensionality());

  for (int i=0; i<num_obj; i++) {
    ZOLTAN_ID_PTR lid = lids + i*num_lid_entries;
    int entity_index = AcmeZoltanLID::Index(lid);
    PRECONDITION(entity_index>=0);
    PRECONDITION(entity_index<topology->NumNodes());
    POSTCONDITION(Nodes[entity_index]);
    for (int j=0; j<num_dim; j++) {
      geom_vec[i*num_dim+j] = Nodes[entity_index]->Coordinates()[j];
    }
  }
  *ierr = ZOLTAN_OK;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      The following three routines are for migrating
//      entities for the dynamic load balance.
//
//-------------------------------------------------------------------------

void AcmeMigrateEntitySizes(void *data,  
			    int num_gid_entries, int num_lid_entries, 
			    int num_ids,
			    ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
			    int* sizes,
			    int *ierr )
{
  DriverTopology* topology = (DriverTopology *)data;
  PRECONDITION(num_gid_entries==ACME_ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ACME_ZOLTAN_LID_SIZE);
  AcmeNode** nodes = 
    REINTERPRET_CAST(AcmeNode**)(topology->NodeList()->EntityList());
  AcmeFace** faces = 
    REINTERPRET_CAST(AcmeFace**)(topology->FaceList()->EntityList());
  AcmeElem** elems = 
    REINTERPRET_CAST(AcmeElem**)(topology->ElemList()->EntityList());

  for( int i=0 ; i<num_ids ; i++ ){
    ZOLTAN_ID_PTR lid = lids+num_lid_entries*i;
    int entity_type   = AcmeZoltanLID::Type(lid);
    int entity_index  = AcmeZoltanLID::Index(lid);
    int& size         = sizes[i];
    switch (entity_type) {
    case AcmeEntity::CT_NODE:
      {
	AcmeNode* node = nodes[entity_index];
        POSTCONDITION(node);
	size = node->Size();
	break;
      }
    case AcmeEntity::CT_FACE:
      {
	AcmeFace* face = faces[entity_index];
        POSTCONDITION(face);
	size = face->Size();
	break;
      }
    case AcmeEntity::CT_ELEM:
      {
	AcmeElem* elem = elems[entity_index];
        POSTCONDITION(elem);
	size = elem->Size();
	break;
      }
    default:
      POSTCONDITION(false);
      break;
    }
  }
  *ierr = ZOLTAN_OK;
}

void AcmeMigrateEntityMultiPack(void *data, 
			        int num_gid_entries, int num_lid_entries,
			        int num_ids,
			        ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
			        int* procs, int* sizes, int* indices,
			        char *Buf, int *ierr)
{
  PRECONDITION(num_gid_entries==ACME_ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ACME_ZOLTAN_LID_SIZE);
  DriverTopology* topology = (DriverTopology *)data;
  AcmeNode** nodes = 
    REINTERPRET_CAST(AcmeNode**)(topology->NodeList()->EntityList());
  AcmeFace** faces = 
    REINTERPRET_CAST(AcmeFace**)(topology->FaceList()->EntityList());
  AcmeElem** elems = 
    REINTERPRET_CAST(AcmeElem**)(topology->ElemList()->EntityList());
  
  for( int i=0 ; i<num_ids ; i++ ){
    ZOLTAN_ID_PTR lid = lids + num_lid_entries*i;
    int entity_type   = AcmeZoltanLID::Type(lid);
    int entity_index  = AcmeZoltanLID::Index(lid);
    char* buf         = Buf+indices[i];
    switch (entity_type) {
    case AcmeEntity::CT_NODE:
      {
        AcmeNode* node = nodes[entity_index];
	POSTCONDITION(node);
        node->Pack(buf);
      }
      break;
    case AcmeEntity::CT_FACE:
      {
        AcmeFace* face = faces[entity_index];
	POSTCONDITION(face);
	face->Pack(buf);
      }
      break;
    case AcmeEntity::CT_ELEM:
      {
        AcmeElem* elem = elems[entity_index];
	POSTCONDITION(elem);
	elem->Pack(buf);
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
  }
  *ierr = ZOLTAN_OK;
}

void AcmeMigrateEntityMultiUnpack(void *data, int num_gid_entries,
				  int num_ids, ZOLTAN_ID_PTR gids, 
				  int* sizes,  int* indices,
				  char *Buf, int *ierr)
{
  PRECONDITION(num_gid_entries==ACME_ZOLTAN_GID_SIZE);
  DriverTopology* topology = (DriverTopology *)data;

  for( int i=0 ; i<num_ids ; i++ ){
    ZOLTAN_ID_PTR gid = gids + i*num_gid_entries;
    int entity_type   = AcmeZoltanGID::Type(gid);
    char* buf         = Buf + indices[i];
    int*  ibuf        = REINTERPRET_CAST(int*) (buf);
    int block         = ibuf[AcmeEntity::BLOCK_ID];
    switch (entity_type) {
    case AcmeEntity::CT_NODE:
      PRECONDITION( ibuf[0] == AcmeEntity::CT_NODE);
      PRECONDITION( block>=0 && block<topology->NumNodeBlocks() );
      topology->NodeBlock(block)->Insert_Entity(buf);
    break;
    case AcmeEntity::CT_FACE:
      PRECONDITION( ibuf[0] == AcmeEntity::CT_FACE );
      PRECONDITION( block>=0 && block<topology->NumFaceBlocks() );
      topology->FaceBlock(block)->Insert_Entity(buf);
    break;
    case AcmeEntity::CT_ELEM:
      PRECONDITION( ibuf[0] == AcmeEntity::CT_ELEM );
      PRECONDITION( block>=0 && block<topology->NumElemBlocks() );
      topology->ElemBlock(block)->Insert_Entity(buf);
    break;
    default:
      POSTCONDITION(false);
      break;
    }
  }
  *ierr = ZOLTAN_OK;
}

#endif



