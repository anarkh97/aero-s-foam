// $Id: AcmeNode.C,v 2002.2 2004/06/23 17:08:55 mwglass Exp $

#include "AcmeNode.h"
#include "AcmeFace.h"
#include "AcmeElem.h"
#include "contact_assert.h"
#include <stddef.h>
#include <new.h>

AcmeNode::AcmeNode(AcmeEntity::AcmeNodeType Type,
                   int Block_Index, 
		   int Host_Index_in_Block,
		   int Exo_ID ) 
  : AcmeEntity( AcmeEntity::CT_NODE, Block_Index, Host_Index_in_Block, Exo_ID)
{
  type = Type;
  number_face_connections = 0;
  faces = NULL;
  number_elem_connections = 0;
  elems = NULL;
}

AcmeNode::~AcmeNode()
{
  if( elems ) delete [] elems;
  if( faces ) delete [] faces;
}

void AcmeNode::Delete_Elem_Connections( )
{
  if( elems ) delete [] elems;
  elems = NULL;
}

void AcmeNode::Allocate_Elem_Connections( int number )
{
  PRECONDITION( number_elem_connections >= 0 );
  number_elem_connections = number;
  if( number_elem_connections ){
    elems = new AcmeElem* [number_elem_connections];
    for( int i=0 ; i<number_elem_connections ; i++ )
      elems[i] = NULL;
  }
}

void AcmeNode::Connect_Elem( int number, AcmeElem* Elem )
{
  PRECONDITION( number>=0 && number<number_elem_connections );
  PRECONDITION( Elem );
  elems[number] = Elem;
}

void AcmeNode::Delete_Face_Connections( )
{
  if( faces ) delete [] faces;
  faces = NULL;
}

void AcmeNode::Allocate_Face_Connections( int number )
{
  PRECONDITION( number_face_connections >= 0 );
  number_face_connections = number;
  if( number_face_connections ){
    faces = new AcmeFace* [number];
    for( int i=0 ; i<number_face_connections ; i++ ) {
      faces[i] = NULL;
    }
  }
}

void AcmeNode::Connect_Face( int number, AcmeFace* Face )
{
  PRECONDITION( number>=0 && number<number_face_connections );
  PRECONDITION( Face );
  faces[number] = Face;
}
