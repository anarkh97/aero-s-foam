// $Id: AcmeEntityBlock.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeEntityBlock.h"
#include "AcmeNode.h"
#include "AcmeFace.h"
#include "AcmeElem.h"
#include "AcmeQuadFaceL4.h"
#include "AcmeQuadFaceQ8.h"
#include "AcmeTriFaceL3.h"
#include "AcmeTriFaceQ6.h"
#include "AcmeHexElemL8.h"

#include "contact_assert.h"

#include <iostream.h>

AcmeEntityBlock::AcmeEntityBlock( AcmeEntity::BaseEntityType Base_Type,
                                  AcmeEntity::DerivedEntityType Derived_Type, 
				  int  block_id,
				  int  num_entities,
                                  int* host_ids)
  : base_type(Base_Type), derived_type(Derived_Type), 
    id(block_id), number_of_entities(num_entities), 
    num_attributes(0), shell_layer_offset(0.0),is_a_shell(false)
{
  int j;
  
  entity_list        = new AcmeBlockEntityList(number_of_entities);
  num_entities_added = number_of_entities;
  PRECONDITION(Derived_Type!=AcmeEntity::ACME_UNKNOWN);
  
  switch( Derived_Type ){
  case AcmeEntity::ACME_NODE:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeNode* node = new AcmeNode(AcmeEntity::NODE,block_id,j,host_ids[j]);
      entity_list->Insert(node);
    }
    break;
  case AcmeEntity::ACME_POINT:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeNode* node = new AcmeNode(AcmeEntity::POINT,block_id,j,host_ids[j]);
      entity_list->Insert(node);
    }
    break;
#ifdef CONTACT_2D_SUPPORT
  case AcmeEntity::ACME_LINEFACEL2:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeLineFaceL2(block_id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
  case AcmeEntity::ACME_LINEFACEL3:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeLineFaceL3(block_id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
#endif
  case AcmeEntity::ACME_QUADFACEL4:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeQuadFaceL4(block_id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
  case AcmeEntity::ACME_QUADFACEQ8:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeQuadFaceQ8(block_id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
  case AcmeEntity::ACME_TRIFACEL3:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeTriFaceL3(block_id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
  case AcmeEntity::ACME_TRIFACEQ6:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeTriFaceQ6(block_id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
  case AcmeEntity::ACME_HEXELEML8:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeElem* elem = new AcmeHexElemL8(block_id,j,host_ids[j]);
      entity_list->Insert(elem);
    }
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    cerr << "Unknown entity type in AcmeEntityBlock: " << Derived_Type << endl;
#endif
    POSTCONDITION( 0 );
  }
}

AcmeEntityBlock::~AcmeEntityBlock()
{
  Delete_Entities();
  delete entity_list;
}

void
AcmeEntityBlock::Add_Entities( int num_entities, int* host_ids)
{
  PRECONDITION( num_entities >= 0 );

  int j;
  
  number_of_entities += num_entities;
  
  switch( derived_type ){
  case AcmeEntity::ACME_NODE:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeNode* node = new AcmeNode(AcmeEntity::NODE,id,j,host_ids[j]);
      entity_list->Insert(node);
    }
    break;
  case AcmeEntity::ACME_POINT:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeNode* node = new AcmeNode(AcmeEntity::POINT,id,j,host_ids[j]);
      entity_list->Insert(node);
    }
    break;
#ifdef CONTACT_2D_SUPPORT
  case AcmeEntity::ACME_LINEFACEL2:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeLineFaceL2(id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
  case AcmeEntity::ACME_LINEFACEL3:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeLineFaceL3(id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
#endif
  case AcmeEntity::ACME_QUADFACEL4:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeQuadFaceL4(id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
  case AcmeEntity::ACME_QUADFACEQ8:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeQuadFaceQ8(id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
  case AcmeEntity::ACME_TRIFACEL3:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeTriFaceL3(id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
  case AcmeEntity::ACME_TRIFACEQ6:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeFace* face = new AcmeTriFaceQ6(id,j,host_ids[j]);
      entity_list->Insert(face);
    }
    break;
  case AcmeEntity::ACME_HEXELEML8:
    for( j=0 ; j<num_entities ; j++ ) {
      AcmeElem* elem = new AcmeHexElemL8(id,j,host_ids[j]);
      entity_list->Insert(elem);
    }
    break;
  }
}

void 
AcmeEntityBlock::Delete_Entity_List()
{
  Delete_Entities();
  entity_list->CleanUp();
  number_of_entities = 0;
  num_entities_added = 0;
}

void AcmeEntityBlock::Delete_Entities()
{
  if( entity_list ){
    switch( derived_type ){
    case AcmeEntity::ACME_NODE:
      entity_list->IteratorStart();
      while (AcmeEntity* entity = entity_list->IteratorForward()) {
        AcmeNode* node = static_cast<AcmeNode*>(entity);
	node->~AcmeNode();
      }
      break;
    case AcmeEntity::ACME_POINT:
      entity_list->IteratorStart();
      while (AcmeEntity* entity = entity_list->IteratorForward()) {
        AcmeNode* node = static_cast<AcmeNode*>(entity);
        node->~AcmeNode();
      }
      break;
    }
  }
}

void AcmeEntityBlock::Insert_Entity( AcmeEntity* entity )
{
  PRECONDITION( entity );
  entity_list->Insert(entity);
  number_of_entities = entity_list->NumEntities();
}

#ifndef CONTACT_NO_MPI

void AcmeEntityBlock::Insert_Entity( char* buffer )
{
  PRECONDITION( buffer );
  entity_list->Insert(buffer);
  number_of_entities = entity_list->NumEntities();
}

#endif

void AcmeEntityBlock::Delete_Entity( AcmeEntity* entity )
{
  PRECONDITION( entity );
  entity_list->Delete(entity);
  number_of_entities = entity_list->NumEntities();
}
