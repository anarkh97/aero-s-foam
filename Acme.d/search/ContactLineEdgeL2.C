// $Id$

#include "ContactLineEdgeL2.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactNode.h"
#include <new>

ContactLineEdgeL2::ContactLineEdgeL2( int Blk_Index, int Host_Index_in_Blk ) 
  : ContactEdge( ContactSearch::LINEEDGEL2, Blk_Index, Host_Index_in_Blk, nodes )
{}


ContactLineEdgeL2* ContactLineEdgeL2::new_ContactLineEdgeL2(
		    ContactFixedSizeAllocator& alloc,
		    int Block_Index, int Host_Index_in_Block )
{
  return new (alloc.New_Frag())
    ContactLineEdgeL2( Block_Index, Host_Index_in_Block );
}

void ContactLineEdgeL2_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactLineEdgeL2),
		100,    // block size
		0 );  // initial block size
  alloc.Set_Name( "ContactLineEdgeL2 allocator" );
}

ContactLineEdgeL2::~ContactLineEdgeL2() {}


