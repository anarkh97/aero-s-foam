// $Id$

#include "ContactShellTriFaceL3.h"
#include "ContactFixedSizeAllocator.h"
#include <new>

template<typename DataType>
ContactShellTriFaceL3<DataType>::ContactShellTriFaceL3(ContactFixedSizeAllocator* alloc,
                                             int Block_Index, 
					     int Index_in_Block, int key ) 
  : ContactTriFaceL3<DataType>( alloc, Block_Index, Index_in_Block, key )
{
  // Reset the type (because TriFaceL3 set it to its type)
  this->face_type = ContactSearch::SHELLTRIFACEL3;
}

template<typename DataType>
ContactShellTriFaceL3<DataType>* ContactShellTriFaceL3<DataType>::new_ContactShellTriFaceL3(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactShellTriFaceL3].New_Frag())
             ContactShellTriFaceL3<DataType>(alloc, Block_Index, Index_in_Block, key);
}

template<typename DataType>
void ContactShellTriFaceL3_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactShellTriFaceL3<DataType>),
                100,  // block size
                0,   // initial block size
                sizeof(DataType) );
  alloc.Set_Name( "ContactShellTriFaceL3 allocator" );
}

template<typename DataType>
ContactShellTriFaceL3<DataType>::~ContactShellTriFaceL3() 
{
}

