// $Id: AcmeEntity.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeEntity.h"
#include "AcmeNode.h"
#include "AcmeFace.h"
#include "AcmeElem.h"
#include "contact_assert.h"
#include <stddef.h>
#include <string.h>

AcmeEntity::AcmeEntity( AcmeEntity::BaseEntityType Type, 
                        int Block_ID, 
			int index_in_block,
			int exo_id )
  : entity_type(Type),
    num_data(0), data(NULL), num_attributes(0), attributes(NULL),
    block_id(Block_ID),
    block_array_index(index_in_block), exodus_id(exo_id),
    owner(-1), temp_tag(0), proc_array_index(-1),
    state(ACTIVE), action(NONE)
{ }

AcmeEntity::~AcmeEntity()
{}

int AcmeEntity::Size()
{
  return( NUMBER_PACKED_VARS*sizeof(int)+
          (num_data+num_attributes)*sizeof(Real) );
}

void AcmeEntity::Pack( char* buffer )
{
  // Note that the first location (0) is reserved for the entity type.
  // and the second location (1) is reserved for derived type.
  // The later will be filled in by the derived classes.
  int* i_buffer = REINTERPRET_CAST(int*) (buffer);
  i_buffer[BASE_TYPE]   = entity_type;
  i_buffer[EXO_ID]      = exodus_id;
  i_buffer[OWNERSHIP]   = ownership;
  i_buffer[OWNER]       = owner;
  i_buffer[IS_ACTIVE]   = state;
  i_buffer[TAG]         = temp_tag;
  i_buffer[BLOCK_INDEX] = block_array_index;
  i_buffer[PROC_INDEX]  = proc_array_index;
  i_buffer[HOST_INDEX]  = host_array_index;
  i_buffer[BLOCK_ID]    = block_id;
  i_buffer[NUM_DATA]    = num_data;
  i_buffer[NUM_ATTR]    = num_attributes;
  Real* data_loc = REINTERPRET_CAST(Real*) 
    (buffer+NUMBER_PACKED_VARS*sizeof(int));
  memcpy(data_loc,data,num_data*sizeof(Real));
  Real* attr_loc = REINTERPRET_CAST(Real*) 
    (buffer+NUMBER_PACKED_VARS*sizeof(int)+num_data*sizeof(Real));
  memcpy(attr_loc,attributes,num_attributes*sizeof(Real));
  
}

void AcmeEntity::Unpack( char* buffer )
{
  int* i_buffer = REINTERPRET_CAST(int*) (buffer);
  PRECONDITION( i_buffer[BASE_TYPE] == entity_type );
  exodus_id         = i_buffer[EXO_ID];
  ownership         = (ProcessorOwnership)(i_buffer[OWNERSHIP]);
  owner             = i_buffer[OWNER];
  state             = (EntityState)(i_buffer[IS_ACTIVE]);
  temp_tag          = i_buffer[TAG];
  block_array_index = i_buffer[BLOCK_INDEX];
  proc_array_index  = i_buffer[PROC_INDEX];
  host_array_index  = i_buffer[HOST_INDEX];
  block_id          = i_buffer[BLOCK_ID];
  num_data          = i_buffer[NUM_DATA];
  num_attributes    = i_buffer[NUM_ATTR];
  data       = new Real[num_data];
  attributes = new Real[num_attributes];
  Real* data_loc = REINTERPRET_CAST(Real*) 
    (buffer+NUMBER_PACKED_VARS*sizeof(int));
  memcpy(data,data_loc,num_data*sizeof(Real));
  Real* attr_loc = REINTERPRET_CAST(Real*) 
    (buffer+NUMBER_PACKED_VARS*sizeof(int)+num_data*sizeof(Real));
  memcpy(attributes,attr_loc,num_attributes*sizeof(Real));
}
