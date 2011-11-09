// $Id: AcmeFace.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeFace.h"
#include "AcmeNode.h"
#include "AcmeTopologyEntityList.h"
#include "contact_assert.h"
#include <stddef.h>
#include <iostream.h>
#include <string.h>
#include <math.h>

AcmeFace::AcmeFace(AcmeEntity::AcmeFaceType Type, 
		   int Block_ID, 
                   int Host_Index_in_Block, 
                   int Exo_ID, 
                   AcmeNode** NodeList, 
                   int* NodeIds)
  : AcmeEntity( AcmeEntity::CT_FACE, Block_ID, 
		Host_Index_in_Block, Exo_ID )
{
  type      = Type;
  node_list = NodeList;
  node_ids  = NodeIds;
}

AcmeFace::~AcmeFace()
{
}

void AcmeFace::ConnectNodes(AcmeTopologyEntityList* NodeList)
{
  for (int i=0; i<number_node_connections; i++) {
    AcmeNode* node = static_cast<AcmeNode*>(NodeList->Find(node_ids[i]));
    node_list[i]   = node;
    POSTCONDITION(node);
  }
}

#ifndef CONTACT_NO_MPI
int AcmeFace::Size()
{
  return( AcmeEntity::Size() + 
	  number_node_connections*sizeof(int) );
}

void AcmeFace::Pack( char* buffer )
{
  int* i_buf = REINTERPRET_CAST(int*) (buffer);
  // ContactTopologyEntity packs in location 0 as AcmeFace and here we pack
  // in the derived type in location 1.
  i_buf[1] = type;
  AcmeEntity::Pack( buffer );
  // Add the node ids
  i_buf = REINTERPRET_CAST(int*) (buffer+AcmeEntity::Size());
  for( int i=0 ; i<number_node_connections ; i++ ){
    *i_buf++ = node_ids[i];
  }
}

void AcmeFace::Unpack( char* buffer )
{
  AcmeEntity::Unpack( buffer );
  PRECONDITION( ((int*)buffer)[1] == type );
  // Get the node ids
  int* i_buf = REINTERPRET_CAST(int*) (buffer+AcmeEntity::Size());
  for( int i=0 ; i<number_node_connections ; i++ ){
    node_ids[i] = *i_buf++;
  }
}
#endif
