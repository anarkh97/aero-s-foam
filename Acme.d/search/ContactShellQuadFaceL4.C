// This file is part of a modified version of ACME: Algorithms for Contact in
// a Multiphysics Environment, derived from version 2.7f
//
// Copyright (C) 2007 Sandia Corporation
// Copyright (C) 2011 Stanford University
//
// ACME is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ACME.  If not, see <http://www.gnu.org/licenses/>.


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

