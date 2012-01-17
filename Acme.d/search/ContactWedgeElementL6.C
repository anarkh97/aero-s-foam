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


#include <algorithm>

#include "allocators.h"
#include "ContactNode.h"
#include "ContactLineEdgeL2.h"
#include "ContactTriFaceL3.h"
#include "ContactQuadFaceL4.h"
#include "ContactWedgeElementL6.h"
#include "ContactFixedSizeAllocator.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <new>

ContactWedgeElemL6::ContactWedgeElemL6( int Block_Index, 
				        int Host_Index_in_Block, int key ) 
  : ContactElem( ContactSearch::WEDGEELEML6,Block_Index,Host_Index_in_Block,key) 
{
  for (int i=0; i<Nodes_Per_Element(); ++i) {
    nodes[i]    = NULL;
    node_ids[i] = -1;
  }
  for (int i=0; i<Edges_Per_Element(); ++i) {
    edges[i]    = NULL;
    edge_ids[i] = -1;
  }
  for (int i=0; i<Faces_Per_Element(); ++i) {
    faces[i]    = NULL;
    face_ids[i] = -1;
  }
}

ContactWedgeElemL6* ContactWedgeElemL6::new_ContactWedgeElemL6(
                        ContactFixedSizeAllocator& alloc,
                        int Block_Index, int Host_Index_in_Block, int key)
{
  return new (alloc.New_Frag())
             ContactWedgeElemL6(Block_Index, Host_Index_in_Block, key);
}

void ContactWedgeElemL6_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactWedgeElemL6),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactWedgeElemL6 allocator" );
}

ContactWedgeElemL6::~ContactWedgeElemL6() {}

void ContactWedgeElemL6::BuildTopology(int nID, int eID, int fID,
				       ContactFixedSizeAllocator* allocators)
{
  int i;
  ContactNode* node;
  ContactEdge* edge;
  ContactFace* face;
  
  int NextID = nID;
  for( i=0 ; i<Nodes_Per_Element() ; ++i ) {
    node = ContactNode::new_ContactNode(allocators,
                                        ContactSearch::NODE, 
                                        ++NextID );
    ConnectNode(i, node);
  }
  NextID = eID;
  for( i=0 ; i<Edges_Per_Element() ; ++i ) {
    edge = ContactLineEdgeL2::new_ContactLineEdgeL2( 
                        allocators[ContactSearch::ALLOC_ContactLineEdgeL2],
                        ContactSearch::LINEEDGEL2, ++NextID);
    ConnectEdge(i, edge);
  }
  edge = Edge(0);
  edge->ConnectNode(0, Node(0));
  edge->ConnectNode(1, Node(1));
  edge = Edge(1);
  edge->ConnectNode(0, Node(1));
  edge->ConnectNode(1, Node(4));
  edge = Edge(2);
  edge->ConnectNode(0, Node(4));
  edge->ConnectNode(1, Node(3));
  edge = Edge(3);
  edge->ConnectNode(0, Node(3));
  edge->ConnectNode(1, Node(0));
  edge = Edge(4);
  edge->ConnectNode(0, Node(1));
  edge->ConnectNode(1, Node(2));
  edge = Edge(5);
  edge->ConnectNode(0, Node(2));
  edge->ConnectNode(1, Node(5));
  edge = Edge(6);
  edge->ConnectNode(0, Node(5));
  edge->ConnectNode(1, Node(4));
  edge = Edge(7);
  edge->ConnectNode(0, Node(2));
  edge->ConnectNode(1, Node(0));
  edge = Edge(8);
  edge->ConnectNode(0, Node(3));
  edge->ConnectNode(1, Node(5));
  
  NextID = fID;
  for( i=0 ; i<3 ; ++i ) {
    face = ContactQuadFaceL4::new_ContactQuadFaceL4(allocators, ++NextID );
    ConnectFace(i, face);
  }
  for( i=3 ; i<5 ; ++i ) {
    face = ContactTriFaceL3::new_ContactTriFaceL3(allocators, ++NextID );
    ConnectFace(i, face);
  }
  face = Face(0);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(1));
  face->ConnectNode(2, Node(4));
  face->ConnectNode(3, Node(3));
  face = Face(1);
  face->ConnectNode(0, Node(1));
  face->ConnectNode(1, Node(2));
  face->ConnectNode(2, Node(5));
  face->ConnectNode(3, Node(4));
  face = Face(2);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(3));
  face->ConnectNode(2, Node(5));
  face->ConnectNode(3, Node(2));
  face = Face(3);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(2));
  face->ConnectNode(2, Node(1));
  face = Face(4);
  face->ConnectNode(0, Node(3));
  face->ConnectNode(1, Node(4));
  face->ConnectNode(2, Node(5));
}

void ContactWedgeElemL6::DeleteTopology(ContactFixedSizeAllocator* allocators)
{
  int i;
  for( i=0 ; i<Nodes_Per_Element() ; ++i ) {
    ContactNode* node = Node(i);
    node->~ContactNode();
    allocators[ContactSearch::ALLOC_ContactNode].Delete_Frag(node);
  }
  for( i=0 ; i<Edges_Per_Element() ; ++i ) {
    ContactEdge* edge = Edge(i);
    edge->~ContactEdge();
    allocators[ContactSearch::ALLOC_ContactLineEdgeL2].Delete_Frag(edge);
  }
  for( i=0 ; i<3 ; ++i ) {
    ContactFace* face = Face(i);
    face->~ContactFace();
    allocators[ContactSearch::ALLOC_ContactQuadFaceL4].Delete_Frag(face);
  }
  for( i=3 ; i<Faces_Per_Element() ; ++i ) {
    ContactFace* face = Face(i);
    face->~ContactFace();
    allocators[ContactSearch::ALLOC_ContactTriFaceL3].Delete_Frag(face);
  }
}

void ContactWedgeElemL6::UpdateTopology(ContactFace* face, 
					VariableHandle POSITION,
					VariableHandle FACE_NORMAL,
					VariableHandle NODE_NORMAL,
					Real tol, bool use_node_normals)
{
  int i;
  int num_nodes = face->Nodes_Per_Face();
  for( i=0 ; i<num_nodes ; ++i ){
    Real* projection;
    ContactNode* face_node    = face->Node(i);
    ContactNode* elem_node1   = Node(i);
    ContactNode* elem_node2   = Node(i+num_nodes);
    Real* face_node_position  = face_node->Variable(POSITION);
    Real* elem_node_position1 = elem_node1->Variable(POSITION);
    Real* elem_node_position2 = elem_node2->Variable(POSITION);
    if (use_node_normals) {
      projection = face_node->Variable(NODE_NORMAL);
    } else {
      projection = face->Variable(FACE_NORMAL);
    }
    for( int k=0 ; k<3 ; ++k ){
      elem_node_position1[k] = face_node_position[k]-projection[k]*tol;
      elem_node_position2[k] = face_node_position[k]+projection[k]*tol;
    }
  }
  for( i=0 ; i<Faces_Per_Element() ; ++i ) {
    Face(i)->Compute_Normal(POSITION,FACE_NORMAL);
  }
}

bool ContactWedgeElemL6::Is_Local_Coordinates_Inside_Element( Real* local_coords )
{
  if( local_coords[0] >=  0.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >=  0.0 && local_coords[1] <= 1.0 &&
      local_coords[2] >=  0.0 && local_coords[2] <= 1.0 &&
      local_coords[3] >= -1.0 && local_coords[3] <= 1.0 )
    return true;
  return false;
}

bool ContactWedgeElemL6::Is_Local_Coordinates_Near_Element( Real* local_coords, Real tolerance )
{
  Real low_coord  = -(1.+tolerance);
  Real high_coord = 1.+tolerance;
  if( local_coords[0] >= low_coord && local_coords[0] <= high_coord &&
      local_coords[1] >= low_coord && local_coords[1] <= high_coord &&
      local_coords[2] >= low_coord && local_coords[2] <= high_coord )
    return true;
  return false;
}

void ContactWedgeElemL6::Evaluate_Shape_Functions( Real* local_coords,
						   Real* shape_functions )
{
  Compute_Shape_Functions(local_coords, shape_functions);
}

void ContactWedgeElemL6::Compute_Global_Coordinates( VariableHandle POSITION,
						     Real* local_coords,
						     Real* global_coords )
{
  Real node_positions[6][3];
  for(int i=0; i<Nodes_Per_Element(); ++i ){
    Real* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void ContactWedgeElemL6::Compute_Local_Coordinates( Real Config_Param,
						    VariableHandle POSITION0, 
						    VariableHandle POSITION1, 
						    VariableHandle FACE_NORMAL,
						    Real* global_coords,
						    Real* local_coords )
{
  int i, j;
  Real node_positions[6][3];
  if (Config_Param == 0.0) {
    for (i=0; i<Nodes_Per_Element(); ++i) {
      Real* node_position = Node(i)->Variable(POSITION0);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else if (Config_Param == 1.0) {
    for (i=0; i<Nodes_Per_Element(); ++i) {
      Real* node_position = Node(i)->Variable(POSITION1);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else {
    Real alpha = 1.0 - Config_Param, beta = Config_Param;
    for (i=0; i<Nodes_Per_Element(); ++i) {
      Real* node_position0 = Node(i)->Variable(POSITION0);
      Real* node_position1 = Node(i)->Variable(POSITION1);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = alpha*node_position0[j]+beta*node_position1[j];
      }
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

void ContactWedgeElemL6::Compute_Local_Coordinates( VariableHandle POSITION,
						    Real* global_coords,
						    Real* local_coords )
{
  int i, j;
  Real node_positions[6][3];
  for (i=0; i<Nodes_Per_Element(); ++i) {
    Real* node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

/*************************************************************************/
/*************************************************************************/
/*                                                                       */
/* The following functions are supplied for doing generic computations   */
/* on a linear W6 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary W6 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

void ContactWedgeElemL6::Compute_Shape_Functions( Real* local_coords,
						  Real* shape_functions )
{
  shape_functions[0] = 0.50*local_coords[0]*(1.0-local_coords[3]);
  shape_functions[1] = 0.50*local_coords[1]*(1.0-local_coords[3]);
  shape_functions[2] = 0.50*local_coords[2]*(1.0-local_coords[3]);
  shape_functions[3] = 0.50*local_coords[0]*(1.0+local_coords[3]);
  shape_functions[4] = 0.50*local_coords[1]*(1.0+local_coords[3]);
  shape_functions[5] = 0.50*local_coords[2]*(1.0+local_coords[3]);
}

void ContactWedgeElemL6::Compute_Shape_Derivatives( Real* local_coords,
						    Real  shape_derivs[3][6] )
{
  shape_derivs[0][0] =  0.50*(1.0-local_coords[3]);
  shape_derivs[0][1] =  0.0;
  shape_derivs[0][2] = -0.50*(1.0-local_coords[3]);
  shape_derivs[0][3] =  0.50*(1.0+local_coords[3]);
  shape_derivs[0][4] =  0.0;
  shape_derivs[0][5] = -0.50*(1.0+local_coords[3]);

  shape_derivs[1][0] =  0.0;
  shape_derivs[1][1] =  0.50*(1.0-local_coords[3]);
  shape_derivs[1][2] = -0.50*(1.0-local_coords[3]);
  shape_derivs[1][3] =  0.0;
  shape_derivs[1][4] =  0.50*(1.0+local_coords[3]);
  shape_derivs[1][5] = -0.50*(1.0+local_coords[3]);

  shape_derivs[2][0] = -0.50*local_coords[0];
  shape_derivs[2][1] = -0.50*local_coords[1];
  shape_derivs[2][2] = -0.50*local_coords[2];
  shape_derivs[2][3] =  0.50*local_coords[0];
  shape_derivs[2][4] =  0.50*local_coords[1];
  shape_derivs[2][5] =  0.50*local_coords[2];
}

void ContactWedgeElemL6::Compute_Local_Coords( Real node_positions[8][3], 
					       Real* global_coords,
					       Real* local_coords )
{
  //
  // 1st check for coincidence with one of the face nodes
  //
  int  i, j;
  int  nnodes=6;
  Real spatial_tolerance = 1.0e-10;
  for (i=0; i<nnodes; ++i) {
    Real dx = node_positions[i][0]-global_coords[0];
    Real dy = node_positions[i][1]-global_coords[1];
    Real dz = node_positions[i][2]-global_coords[2];
    Real d  = std::sqrt(dx*dx+dy*dy+dz*dz);
    if (d<spatial_tolerance) break;
  }
  switch (i) {
  case 0:
    local_coords[0] =  1.0;
    local_coords[1] =  0.0;
    local_coords[2] =  0.0;
    local_coords[3] = -1.0;
    break;
  case 1:
    local_coords[0] =  0.0;
    local_coords[1] =  1.0;
    local_coords[2] =  0.0;
    local_coords[3] = -1.0;
    break;
  case 2:
    local_coords[0] =  0.0;
    local_coords[1] =  0.0;
    local_coords[2] =  1.0;
    local_coords[3] = -1.0;
    break;
  case 3:
    local_coords[0] =  1.0;
    local_coords[1] =  0.0;
    local_coords[2] =  0.0;
    local_coords[3] =  1.0;
    break;
  case 4:
    local_coords[0] =  0.0;
    local_coords[1] =  1.0;
    local_coords[2] =  0.0;
    local_coords[3] =  1.0;
    break;
  case 5:
    local_coords[0] =  0.0;
    local_coords[1] =  0.0;
    local_coords[2] =  1.0;
    local_coords[3] =  1.0;
    break;
  }
  if (i<nnodes) return;
  //
  // else use newton's method to iterate
  //
  int  iterations=0;
  int  max_iterations=200;
  bool converged = false;
  Real tolerance = 1.0e-12;
  Real u, u0=0.0, u1, du;
  Real v, v0=0.0, v1, dv;
  Real w, w0=0.0, w1, dw;
  Real f[3], J[3][3], invJ[3][3];
  Real shape_derivatives[3][6], coords[4];
  while (!converged && iterations<max_iterations) {
    coords[0] = u0;
    coords[1] = v0;
    coords[2] = 1.0-u0-v0;
    coords[3] = w0;
    Compute_Global_Coords(node_positions , coords, f );
    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives( coords, shape_derivatives );
    for (i=0; i<3; ++i) {
      J[0][i] = 0.0;
      J[1][i] = 0.0;
      J[2][i] = 0.0;
      for (j=0; j<6; ++j) {
        J[0][i] += shape_derivatives[i][j]*node_positions[j][0];
        J[1][i] += shape_derivatives[i][j]*node_positions[j][1];
        J[2][i] += shape_derivatives[i][j]*node_positions[j][2];
      }
    }
    
    Real detJ  =  1.0/(J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-
                       J[0][1]*(J[1][0]*J[2][2]-J[2][0]*J[1][2])+
                       J[0][2]*(J[1][0]*J[2][1]-J[2][0]*J[1][1]));

    invJ[0][0] =  (J[1][1]*J[2][2]-J[1][2]*J[2][1])*detJ;
    invJ[0][1] = -(J[0][1]*J[2][2]-J[2][1]*J[0][2])*detJ;
    invJ[0][2] =  (J[1][2]*J[0][1]-J[0][2]*J[1][1])*detJ;
    invJ[1][0] = -(J[1][0]*J[2][2]-J[2][0]*J[1][2])*detJ;
    invJ[1][1] =  (J[0][0]*J[2][2]-J[0][2]*J[2][0])*detJ;
    invJ[1][2] = -(J[0][0]*J[1][2]-J[1][0]*J[0][2])*detJ;
    invJ[2][0] =  (J[1][0]*J[2][1]-J[2][0]*J[1][1])*detJ;
    invJ[2][1] = -(J[0][0]*J[2][1]-J[2][0]*J[0][1])*detJ;
    invJ[2][2] =  (J[0][0]*J[1][1]-J[0][1]*J[1][0])*detJ;

    // APPLY NEWTON ALGORITHM

    u  = f[0]-global_coords[0];
    v  = f[1]-global_coords[1];
    w  = f[2]-global_coords[2];
    u1 = u0-(invJ[0][0]*u+invJ[0][1]*v+invJ[0][2]*w);
    v1 = v0-(invJ[1][0]*u+invJ[1][1]*v+invJ[1][2]*w);
    w1 = w0-(invJ[2][0]*u+invJ[2][1]*v+invJ[2][2]*w);
    du = std::fabs(u1-u0);
    dv = std::fabs(v1-v0);
    dw = std::fabs(w1-w0);
    u0 = u1;
    v0 = v1;
    w0 = w1;
    if (du<tolerance && dv<tolerance && dw<tolerance) converged = true;
    ++iterations;
  }
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  if (!converged) {
    std::cerr << "ContactWedgeElemL6::Compute_Local_Coordinates() did not converge" 
	 << std::endl;
  }
#endif
  POSTCONDITION(converged);
  if (u0<1.0+spatial_tolerance) {
    u0 = std::min(u0, 1.0);
  }
  if (u0>-spatial_tolerance) {
    u0 = std::max(u0, 0.0);
  }
  if (v0<1.0+spatial_tolerance) {
    v0 = std::min(v0, 1.0);
  }
  if (v0>-spatial_tolerance) {
    v0 = std::max(v0, 0.0);
  }
  if (std::fabs(w0)<1.0+spatial_tolerance) {
    w0 = std::min(w0, 1.0);
    w0 = std::max(w0,-1.0);
  }
  local_coords[0] = u0;
  local_coords[1] = v0;
  local_coords[2] = 1.0-u0-v0;
  local_coords[3] = w0;
}

void ContactWedgeElemL6::Compute_Global_Coords( Real node_positions[6][3],
						Real local_coords[4],
						Real global_coords[3] )
{
  Real N[6];
  int  nnodes=6;
  global_coords[0] = 0.0;
  global_coords[1] = 0.0;
  global_coords[2] = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    for (int j=0; j<3; ++j) {
      global_coords[j] += N[i]*node_positions[i][j];
    }
  }
}

void ContactWedgeElemL6::Interpolate_Scalar( Real  local_coords[4],
					     Real  node_scalars[6],
					     Real& interpolated_scalar )
{
  Real N[6];
  int  nnodes=6;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

void ContactWedgeElemL6::Interpolate_Vector( Real local_coords[4],
					     Real node_vectors[6][3],
					     Real interpolated_vector[3] )
{
  Real N[6];
  int  nnodes=6;
  interpolated_vector[0] = 0.0;
  interpolated_vector[1] = 0.0;
  interpolated_vector[2] = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    for (int j=0; j<3; ++j) {
      interpolated_vector[j] += N[i]*node_vectors[i][j];
    }
  }
}



