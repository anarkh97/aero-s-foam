// $Id$

#include "ContactShellQuadFaceL4.h"
#include "ContactFixedSizeAllocator.h"
#include <new>

template<typename DataType>
ContactShellQuadFaceL4<DataType>::ContactShellQuadFaceL4(ContactFixedSizeAllocator* alloc,
                                               int Block_Index, 
					       int Index_in_Block, int key ) 
  : ContactQuadFaceL4<DataType>( alloc, Block_Index, Index_in_Block, key )
{
  // Reset the type (because QuadFaceL4 set it to its type)
  this->face_type = ContactSearch::SHELLQUADFACEL4;
}

template<typename DataType>
ContactShellQuadFaceL4<DataType>* ContactShellQuadFaceL4<DataType>::new_ContactShellQuadFaceL4(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactShellQuadFaceL4].New_Frag())
             ContactShellQuadFaceL4<DataType>(alloc, Block_Index, 
                                    Index_in_Block, key);
}

template<typename DataType>
void ContactShellQuadFaceL4_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactShellQuadFaceL4<DataType>),
                100,  // block size
                0,   // initial block size
                sizeof(DataType) );
  alloc.Set_Name( "ContactShellQuadFaceL4 allocator" );
}

template<typename DataType>
ContactShellQuadFaceL4<DataType>::~ContactShellQuadFaceL4() 
{
}

