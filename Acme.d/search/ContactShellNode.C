// $Id$

#include "ContactShellNode.h"
#include "ContactFixedSizeAllocator.h"
#include <cstddef>
#include <new>

ContactShellNode::ContactShellNode( ContactFixedSizeAllocator* alloc,
                                    ContactSearch::ContactNode_Type Type, 
				    int Block_Index, 
				    int Index_in_Block )
  : ContactNode( alloc, Type, Block_Index, Index_in_Block, 
		 CT_SHELL_NODE )
{
  previous_lofting[0] = 0.0;
  previous_lofting[1] = 0.0;
  previous_lofting[2] = 0.0;
}

ContactShellNode::~ContactShellNode()
{
}

bool ContactShellNode::ConnectedToAllShellFaces() {
  int num_faces = Number_Face_Connections();
  for(int iface = 0; iface < num_faces; ++iface) {
    if(!ContactSearch::Is_a_Shell_Face(GetFace(iface)->FaceType())) return false;
  }
  return true;
}


ContactShellNode* 
ContactShellNode::new_ContactShellNode( ContactFixedSizeAllocator* alloc,
					ContactSearch::ContactNode_Type Type,
					int Block_Index, 
					int Index_in_Block )
{
  return new (alloc[ContactSearch::ALLOC_ContactShellNode].New_Frag()) 
    ContactShellNode( alloc, Type, Block_Index, Index_in_Block );
}

void ContactShellNode_SizeAllocator( ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactShellNode),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactShellNode allocator" );
}
