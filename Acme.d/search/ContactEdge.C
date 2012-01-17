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


#include "ContactEdge.h"
#include "ContactNode.h"
#include "ContactFace.h"
#include <cstddef>
#include <iostream>
#include <cstring>
#include <cmath>

bool ContactEdge::array_init = false;
int ContactEdge::NODES_PER_EDGE[ContactSearch::NEDGE_TYPES] = {0};

ContactEdge::ContactEdge( ContactSearch::ContactEdge_Type Type, 
                          int Block_Index, int Host_Index_in_Block,
                          ContactNode **node_list_) 
  : ContactTopologyEntity( Block_Index,Host_Index_in_Block, DataArray, CT_EDGE),
    node_list(node_list_)
{
  edge_type      = Type;
  number_face_connections = 0;
  face_info.owner = -1; 
  face_info.owner_proc_array_index = -1;
  face_info.host_gid[0] = -1; 
  face_info.host_gid[1] = -1;
}

ContactEdge::~ContactEdge() {}

void ContactEdge::ConnectFace( ContactFace* face )
{
  PRECONDITION( number_face_connections>=0 && number_face_connections<2 );
  PRECONDITION( face );
  faces[number_face_connections++] = face;
}

void ContactEdge::Initialize_Lookup_Arrays() {
  NODES_PER_EDGE[ContactSearch::NO_EDGES  ] = 0;
  NODES_PER_EDGE[ContactSearch::LINEEDGEL2] = 2;
  NODES_PER_EDGE[ContactSearch::LINEEDGEQ3] = 3;
  array_init = true;
}

void ContactEdge::Smooth_Normal( VariableHandle FACE_NORMAL, 
				 Real* coords, Real* smooth_normal )
{
  smooth_normal[0]    = 0.0;
  smooth_normal[1]    = 0.0;
  smooth_normal[2]    = 0.0;
  for( int i=0 ; i<number_face_connections ; ++i ){
    Real* face_normal = faces[i]->Variable(FACE_NORMAL);
    smooth_normal[0] += face_normal[0];
    smooth_normal[1] += face_normal[1];
    smooth_normal[2] += face_normal[2];
  }
  if (!shared && number_face_connections>1) {
    Real mag = std::sqrt( smooth_normal[0]*smooth_normal[0] +
                          smooth_normal[1]*smooth_normal[1] +
                          smooth_normal[2]*smooth_normal[2] );
    if( mag > 0.0 ){
      mag = 1.0/mag;
      smooth_normal[0] *= mag;
      smooth_normal[1] *= mag;
      smooth_normal[2] *= mag;
    }
  }
}

void ContactEdge::Compute_Smoothed_Normal( VariableHandle FACE_NORMAL )
{
  Real* smooth_normal = &DataArray[NUMBER_SCALAR_VARS+3*Smoothed_Normal];
  smooth_normal[0] = 0.0;
  smooth_normal[1] = 0.0;
  smooth_normal[2] = 0.0;
  for( int i=0 ; i<number_face_connections ; ++i ){
    Real* face_normal = faces[i]->Variable(FACE_NORMAL);
    smooth_normal[0] += face_normal[0];
    smooth_normal[1] += face_normal[1];
    smooth_normal[2] += face_normal[2];
  }
  Real mag = std::sqrt( smooth_normal[0]*smooth_normal[0] +
		        smooth_normal[1]*smooth_normal[1] +
		        smooth_normal[2]*smooth_normal[2] );
  if( mag > 0.0 ){
    mag = 1.0/mag;
    smooth_normal[0] *= mag;
    smooth_normal[1] *= mag;
    smooth_normal[2] *= mag;
  }
}

void ContactEdge::Pack( char* buffer )
{
  int* i_buf = reinterpret_cast<int*>(buffer);
  *i_buf = edge_type;
  // ContactTopologyEntity packs in location 0 as ContactFace and here we pack
  // in the derived type in location 1.
  i_buf[1] = edge_type;
  ContactTopologyEntity::Pack( buffer, DataArray_Length() );
  // Add the global ids of the nodes
  i_buf = reinterpret_cast<int*>(buffer+ContactTopologyEntity::Size(DataArray_Length()));
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Edge() ; ++i ){
    ContactTopologyEntity* entity = 
      static_cast<ContactTopologyEntity*>(Node(i));
    cnt += PackConnection(entity, &i_buf[cnt]);
  }
}

void ContactEdge::Copy( ContactEdge* src )
{
  ContactTopologyEntity::Copy( src, DataArray_Length() );
  face_info = src->face_info;
}
