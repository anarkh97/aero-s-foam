// $Id$

#include "ContactShellTriFaceL3.h"
#include "ContactFixedSizeAllocator.h"
#include <new>

ContactShellTriFaceL3::ContactShellTriFaceL3(ContactFixedSizeAllocator* alloc,
                                             int Block_Index, 
					     int Index_in_Block, int key ) 
  : ContactTriFaceL3<Real>( alloc, Block_Index, Index_in_Block, key )
{
  // Reset the type (because TriFaceL3 set it to its type)
  face_type = ContactSearch::SHELLTRIFACEL3;
}

ContactShellTriFaceL3* ContactShellTriFaceL3::new_ContactShellTriFaceL3(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactShellTriFaceL3].New_Frag())
             ContactShellTriFaceL3(alloc, Block_Index, Index_in_Block, key);
}

void ContactShellTriFaceL3_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactShellTriFaceL3),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactShellTriFaceL3 allocator" );
}

ContactShellTriFaceL3::~ContactShellTriFaceL3() 
{
}

