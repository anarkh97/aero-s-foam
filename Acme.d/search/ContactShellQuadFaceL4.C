// $Id$

#include "ContactShellQuadFaceL4.h"
#include "ContactFixedSizeAllocator.h"
#include <new>

ContactShellQuadFaceL4::ContactShellQuadFaceL4(ContactFixedSizeAllocator* alloc,
                                               int Block_Index, 
					       int Index_in_Block, int key ) 
  : ContactQuadFaceL4( alloc, Block_Index, Index_in_Block, key )
{
  // Reset the type (because QuadFaceL4 set it to its type)
  face_type = ContactSearch::SHELLQUADFACEL4;
}

ContactShellQuadFaceL4* ContactShellQuadFaceL4::new_ContactShellQuadFaceL4(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactShellQuadFaceL4].New_Frag())
             ContactShellQuadFaceL4(alloc, Block_Index, 
                                    Index_in_Block, key);
}

void ContactShellQuadFaceL4_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactShellQuadFaceL4),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactShellQuadFaceL4 allocator" );
}

ContactShellQuadFaceL4::~ContactShellQuadFaceL4() 
{
}

