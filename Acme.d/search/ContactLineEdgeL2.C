// $Id$

#ifndef ContactLineEdgeL2_C_
#define ContactLineEdgeL2_C_

#include "ContactLineEdgeL2.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactNode.h"
#include <new>

template<typename DataType>
ContactLineEdgeL2<DataType>::ContactLineEdgeL2( int Blk_Index, int Host_Index_in_Blk ) 
  : ContactEdge<DataType>( ContactSearch::LINEEDGEL2, Blk_Index, Host_Index_in_Blk, nodes )
{}


template<typename DataType>
ContactLineEdgeL2<DataType>* ContactLineEdgeL2<DataType>::new_ContactLineEdgeL2(
		    ContactFixedSizeAllocator& alloc,
		    int Block_Index, int Host_Index_in_Block )
{
  return new (alloc.New_Frag())
    ContactLineEdgeL2<DataType>( Block_Index, Host_Index_in_Block );
}

template<typename DataType>
void ContactLineEdgeL2_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactLineEdgeL2<DataType>),
		100,    // block size
		0 );  // initial block size
  alloc.Set_Name( "ContactLineEdgeL2<DataType> allocator" );
}

template<typename DataType>
ContactLineEdgeL2<DataType>::~ContactLineEdgeL2() {}

#endif  // #define ContactLineEdgeL2_C_
