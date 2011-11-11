// $Id$

#include "ContactLineEdgeQ3.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactNode.h"
#include <new>

ContactLineEdgeQ3::ContactLineEdgeQ3( int Blk_Index, int Index_in_Blk ) 
  : ContactEdge<Real>( ContactSearch::LINEEDGEQ3, Blk_Index, Index_in_Blk, nodes )
{}

ContactLineEdgeQ3* ContactLineEdgeQ3::new_ContactLineEdgeQ3(
		    ContactFixedSizeAllocator& alloc,
		    int Block_Index, int Index_in_Block )
{
  return new (alloc.New_Frag())
    ContactLineEdgeQ3( Block_Index, Index_in_Block );
}

void ContactLineEdgeQ3_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactLineEdgeQ3),
		100,    // block size
		0 );  // initial block size
  alloc.Set_Name( "ContactLineEdgeQ3 allocator" );
}

ContactLineEdgeQ3::~ContactLineEdgeQ3() {}
