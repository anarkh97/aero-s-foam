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
#include "ContactUtilities.h"
#include "ContactTriFaceL3.h"
#include "ContactQuadFaceL4.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactNode.h"
#include <iostream>
#include <cmath>
#include <new>

using acme::Dot;
using acme::Cross;
using acme::Normalize;
using acme::Magnitude;
using acme::ScalarTripleProduct;

ContactTriFaceL3::ContactTriFaceL3( ContactFixedSizeAllocator* alloc,
                                    int Block_Index, 
				    int Index_in_Block,int key ) 
  : ContactFace( alloc, ContactSearch::TRIFACEL3,
                 Block_Index, Index_in_Block, key, 
                 nodes, edges, Node_Info, Edge_Info) 
{}

ContactTriFaceL3* ContactTriFaceL3::new_ContactTriFaceL3(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactTriFaceL3].New_Frag())
             ContactTriFaceL3(alloc, Block_Index, Index_in_Block, key);
}

void ContactTriFaceL3_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactTriFaceL3),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactTriFaceL3 allocator" );
}

ContactTriFaceL3::~ContactTriFaceL3() {}

void ContactTriFaceL3::Compute_Normal( VariableHandle CURRENT_POSITION,
				       VariableHandle FACE_NORMAL )
{
  ContactNode* node0 = Node(0);
  ContactNode* node1 = Node(1);
  ContactNode* node2 = Node(2);

  Real* Position0 = node0->Variable(CURRENT_POSITION);
  Real* Position1 = node1->Variable(CURRENT_POSITION);
  Real* Position2 = node2->Variable(CURRENT_POSITION);

  // Compute the vector from node 0 to node 1 & from node 0 to node 2
  Real Vec01[3],Vec02[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];

  // Compute the face normal as the cross product of the two vectors
  Real* face_normal = Variable(FACE_NORMAL);

  Cross(Vec01, Vec02, face_normal);
  Normalize(face_normal);
}

void ContactTriFaceL3::Compute_Normal(VariableHandle POSITION,
				       Real* normal, Real* local_coords )
{
  ContactNode* node0 = Node(0);
  ContactNode* node1 = Node(1);
  ContactNode* node2 = Node(2);

  Real* Position0    = node0->Variable(POSITION);
  Real* Position1    = node1->Variable(POSITION);
  Real* Position2    = node2->Variable(POSITION);

  // Compute the vector from node 0 to node 1 & from node 0 to node 2
  Real Vec01[3],Vec02[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];

  // Compute the face normal as the cross product of the two vectors
  Cross(Vec01,Vec02,normal);
  Normalize(normal);
}

//
// KHP:Q2: The normal vector variable is the third argument in this
// KHP:Q2: method while it is the second argument in the previous Compute_Normal
// KHP:Q2: method. It would be nice to have the arguments be in a consistent order.
//
void ContactTriFaceL3::Compute_Normal(Real** nodal_positions,
				      Real* local_coords, 
				      Real* normal )
{
  Real* Position0    = nodal_positions[0];
  Real* Position1    = nodal_positions[1];
  Real* Position2    = nodal_positions[2];

  // Compute the vector from node 0 to node 1 & from node 0 to node 2
  Real Vec01[3],Vec02[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];

  // Compute the face normal as the cross product of the two vectors
  Cross(Vec01,Vec02,normal);
  Normalize(normal);
}

void ContactTriFaceL3::Compute_CharacteristicLength( VariableHandle CURRENT_POSITION,
						     VariableHandle CHARACTERISTIC_LENGTH )
{
  ContactNode* node0 = Node(0);
  ContactNode* node1 = Node(1);
  ContactNode* node2 = Node(2);

  Real* Position0 = node0->Variable(CURRENT_POSITION);
  Real* Position1 = node1->Variable(CURRENT_POSITION);
  Real* Position2 = node2->Variable(CURRENT_POSITION);

  // Compute the vector from node 0 to node 1 & from node 0 to node 2
  Real Vec01[3],Vec02[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];

  // Compute the face normal as the cross product of the two vectors
  Real face_normal[3];
  Cross(Vec01,Vec02,face_normal);

  Real* characteristiclength = Variable(CHARACTERISTIC_LENGTH);
  //
  // KHP:Q3: Couldn't we compute the square of the characteristic length
  // KHP:Q3: and avoid the sqrt operation?
  //
  *characteristiclength = Magnitude(face_normal) * 0.5;
}

void ContactTriFaceL3::Compute_Centroid( VariableHandle CURRENT_POSITION,
					 VariableHandle CENTROID )
{
  Real* centroid = Variable(CENTROID);
  Real* node_position = Node(0)->Variable(CURRENT_POSITION);
  centroid[0] = node_position[0];
  centroid[1] = node_position[1];
  centroid[2] = node_position[2];
  for( int i=1 ; i<3 ; ++i ){
    node_position = Node(i)->Variable(CURRENT_POSITION);
    centroid[0] += node_position[0];
    centroid[1] += node_position[1];
    centroid[2] += node_position[2];
  }
  Real scale   = 1.0/3.0;
  centroid[0] *= scale;
  centroid[1] *= scale;
  centroid[2] *= scale;
}

void ContactTriFaceL3::Get_Edge_Nodes( int i, ContactNode** node )
{
  PRECONDITION( i>=0 && i<3 );
  switch( i ){
  case 0:
    node[0] = Node(0);
    node[1] = Node(1);
    break;
  case 1:
    node[0] = Node(1);
    node[1] = Node(2);
    break;
  case 2:
    node[0] = Node(2);
    node[1] = Node(0);
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Invalid Edge request in ContactTriFaceL3::Get_Face_Nodes" << std::endl;
#endif
    POSTCONDITION( 0 );
  }
}

int ContactTriFaceL3::Get_Edge_Number( ContactNode** edge_nodes )
{
  PRECONDITION( edge_nodes[0] && edge_nodes[1] );
  PRECONDITION( edge_nodes[0] != edge_nodes[1] );
  int i1=-1,i2=-1;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    ContactNode* node = Node(i);
    if( edge_nodes[0] == node ) i1=i;
    if( edge_nodes[1] == node ) i2=i;
  }
  PRECONDITION( 0<=i1 && i1<3 );
  PRECONDITION( 0<=i2 && i2<3 );
  switch( i1 ){
  case 0:
    if( i2 == 1 ) return(0);
    if( i2 == 2 ) return(2);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  case 1:
    if( i2 == 0 ) return(0);
    if( i2 == 2 ) return(1);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  case 2:
    if( i2 == 0 ) return(2);
    if( i2 == 1 ) return(1);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Error: You probably have a bad connectivity" << std::endl;
#endif
    POSTCONDITION( 0 );
    break;
  }
  POSTCONDITION( 0 );
  return(-1);
}

int ContactTriFaceL3::Get_Edge_Number( Real* local_coords )
{
  Real tol(1.0e-8);
  if (std::fabs(local_coords[1])<tol && std::fabs(local_coords[3])<tol) return(0);
  if (std::fabs(local_coords[0])<tol && std::fabs(local_coords[2])<tol) return(1);
  if (std::fabs((1.0-local_coords[0]-local_coords[1]))<tol &&
      std::fabs((1.0-local_coords[2]-local_coords[3]))<tol) return(2);
  return( -1 );
}

void
ContactTriFaceL3::Compute_Edge_Normal( VariableHandle POSITION,
				       VariableHandle FACE_NORMAL,
				       int Edge,
				       Real* edge_normal)
{
  ContactNode *edge_node[2];
  Get_Edge_Nodes( Edge, edge_node );

  // Compute the tangent direction as a vector from node 1 to node 2
  Real* position1 = edge_node[0]->Variable(POSITION);
  Real* position2 = edge_node[1]->Variable(POSITION);
  Real edge_tangent[3];
  edge_tangent[0] = position2[0] - position1[0];
  edge_tangent[1] = position2[1] - position1[1];
  edge_tangent[2] = position2[2] - position1[2];

  // Compute the normal direction as the cross product of the edge tangent
  // and the face normal.
  Real* face_normal = Variable(FACE_NORMAL);
  Cross(edge_tangent, face_normal, edge_normal);
  Normalize(edge_normal);
}

void ContactTriFaceL3::Get_Close_Edges( Real* local_coords, int& number,
					int& edge_1_id, int& edge_2_id ){
  Real a0 = local_coords[0];
  Real a1 = local_coords[1];
  Real a2 = local_coords[2];

  int edge_case = 0;

  if( a0 <= 0.1 ) edge_case += 1;  // Near Edge 1
  if( a1 <= 0.1 ) edge_case += 2;  // Near Edge 2
  if( a2 <= 0.1 ) edge_case += 4;  // Near Edge 0

  switch( edge_case ){
  case(0):
    // not near any edge
    number = 0;
    break;
  case(1):
    // near edge 1 only
    number = 1;
    edge_1_id = 1;
    break;
  case(2):
    // near edge 2 only
    number = 1;
    edge_1_id = 2;
    break;
  case(4):
    // near edge 0 only
    number = 1;
    edge_1_id = 0;
    break;
  case(3):
    // near edges 1 & 2
    number = 2;
    edge_1_id = 1;
    edge_2_id = 2;
    break;
  case(5): 
    // near edge 0 & 1
    number = 2;
    edge_1_id = 0;
    edge_2_id = 1;
    break;
  case(6):
    // near edges 0 & 2
    number = 2;
    edge_1_id = 0;
    edge_2_id = 2;
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "error can only be near at most 2 edges" << std::endl;
#endif
    POSTCONDITION( 0 );
  }
}

bool ContactTriFaceL3::Is_Inside_Face( Real* local_coords )
{
  if( local_coords[0] >= 0.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >= 0.0 && local_coords[1] <= 1.0 &&
      local_coords[2] >= 0.0 && local_coords[2] <= 1.0 )
    return true;
  return false;
}

ContactFace* ContactTriFaceL3::Neighbor( Real* local_coords )
{
  PRECONDITION( 0 );
  return (ContactFace*) NULL;
}

void ContactTriFaceL3::FacetDecomposition(int& nfacets,
          Real* coordinates0, Real* normals0, VariableHandle POSITION0,
          Real* coordinates1, Real* normals1, VariableHandle POSITION1,
          Real* coordinates2, Real* normals2, VariableHandle POSITION2) 
{
  nfacets = 1;
  for(int j=0 ; j<3 ; ++j ){
	Real* node_position  = Node(j)->Variable(POSITION0);
	coordinates0[0+3*j]  = node_position[0];
	coordinates0[1+3*j]  = node_position[1];
	coordinates0[2+3*j]  = node_position[2];
  }

  ContactNode* node0 = Node(0);
  ContactNode* node1 = Node(1);
  ContactNode* node2 = Node(2);

  Real* Position0 = node0->Variable(POSITION0);
  Real* Position1 = node1->Variable(POSITION0);
  Real* Position2 = node2->Variable(POSITION0);

  // Compute the vector from node 0 to node 1 & from node 0 to node 2
  Real Vec01[3],Vec02[3];
  Vec01[0] = Position1[0] - Position0[0];
  Vec01[1] = Position1[1] - Position0[1];
  Vec01[2] = Position1[2] - Position0[2];
  Vec02[0] = Position2[0] - Position0[0];
  Vec02[1] = Position2[1] - Position0[1];
  Vec02[2] = Position2[2] - Position0[2];

  // Compute the face normal as the cross product of the two vectors
  Cross(Vec01,Vec02,normals0);
  Normalize(normals0);
 
  if (coordinates1!=NULL) {
    for(int j=0 ; j<3 ; ++j ){
      Real* node_position = Node(j)->Variable(POSITION1);
      coordinates1[0+3*j]  = node_position[0];
      coordinates1[1+3*j]  = node_position[1];
      coordinates1[2+3*j]  = node_position[2];
    }

    node0     = Node(0);
    node1     = Node(1);
    node2     = Node(2);

    Position0 = node0->Variable(POSITION1);
    Position1 = node1->Variable(POSITION1);
    Position2 = node2->Variable(POSITION1);

    // Compute the vector from node 0 to node 1 & from node 0 to node 2
    Vec01[0] = Position1[0] - Position0[0];
    Vec01[1] = Position1[1] - Position0[1];
    Vec01[2] = Position1[2] - Position0[2];
    Vec02[0] = Position2[0] - Position0[0];
    Vec02[1] = Position2[1] - Position0[1];
    Vec02[2] = Position2[2] - Position0[2];

    // Compute the face normal as the cross product of the two vectors
    Cross(Vec01,Vec02,normals1);
    Normalize(normals1);
  }

  if (coordinates2!=NULL) {
    for(int j=0 ; j<3 ; ++j ){
      Real* node_position = Node(j)->Variable(POSITION2);
      coordinates2[0+3*j]  = node_position[0];
      coordinates2[1+3*j]  = node_position[1];
      coordinates2[2+3*j]  = node_position[2];
    }

    node0     = Node(0);
    node1     = Node(1);
    node2     = Node(2);

    Position0 = node0->Variable(POSITION2);
    Position1 = node1->Variable(POSITION2);
    Position2 = node2->Variable(POSITION2);

    // Compute the vector from node 0 to node 1 & from node 0 to node 2
    Vec01[0] = Position1[0] - Position0[0];
    Vec01[1] = Position1[1] - Position0[1];
    Vec01[2] = Position1[2] - Position0[2];
    Vec02[0] = Position2[0] - Position0[0];
    Vec02[1] = Position2[1] - Position0[1];
    Vec02[2] = Position2[2] - Position0[2];

    // Compute the face normal as the cross product of the two vectors
    Cross(Vec01,Vec02,normals2);
    Normalize(normals2);
  }
}

void ContactTriFaceL3::FacetStaticRestriction(int nfacets, Real* coordinates, 
					      Real* normals, Real* ctrcl_facets,
					      Real* ctrcl)
{
  for (int i=0; i<LENGTH; ++i) {
	ctrcl[i] = ctrcl_facets[i];
  }
}
void ContactTriFaceL3::FacetDynamicRestriction(int nfacets, 
                                               Real* ctrcl_facets, 
                                               Real* ctrcl)
{
  for (int i=0; i<LENGTH; ++i) {
	ctrcl[i] = ctrcl_facets[i];
  }
}



void ContactTriFaceL3::Smooth_Normal( VariableHandle CURRENT_POSITION,
				      VariableHandle NODE_NORMAL,
				      VariableHandle FACE_NORMAL,
				      VariableHandle CURVATURE,
    		       ContactSearch::Smoothing_Resolution resolution,
				      Real percentage,
				      Real* coordinates,
				      Real* smooth_normal,
                                      Real critical_curvature )
{
  /*                   2
                     h/ \g
                     /\G/\
                    /  c  \
                   /  / \  \
                  /  /   \  \
                 /  /     \  \
                / D/   A   \C \
               /  /         \  \
              /  /           \  \
             i--a-------------b--f
            / E/       B       \F \
           0--d-----------------e--1
              
      Region    A: Use surface normal
              B-D: Smooth along one edge
              E-G: Smooth along two edges
  */
  Real upper_bound = 1.0 - percentage;
  Real lower_bound = percentage/2.0;
  Real* face_normal = Variable(FACE_NORMAL);
  Real node_position[4][3];
  Real node_normal[4][3];
  Real contact_point[3];
  Real coords[3];
  bool smooth_edge0,smooth_edge1;
  Real curvature0,curvature1;

  if( coordinates[0] >= lower_bound && coordinates[0] <= upper_bound &&
      coordinates[1] >= lower_bound && coordinates[1] <= upper_bound &&
      coordinates[2] >= lower_bound && coordinates[2] <= upper_bound ){
    // Region A
    smooth_normal[0] = face_normal[0];
    smooth_normal[1] = face_normal[1];
    smooth_normal[2] = face_normal[2];
    return;
  } 
  else if( coordinates[0] >= lower_bound && coordinates[1] >= lower_bound &&
	   coordinates[2] <  lower_bound){
    /* Region B
                  1-------------0
                 /       B       \
                2-----------------3
    */     
    if( fabs(GetEdgeCurvature(0)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    } 
    // Compute the data for node 0
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute the data for node 2
    coords[0] = 1.0-lower_bound;
    coords[1] = lower_bound;
    coords[2] = 0.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[2][0]) );
    GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute the data for node 3
    coords[0] = lower_bound;
    coords[1] = 1.0-lower_bound;
    coords[2] = 0.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    GetEdgeSmoothedNormal(0, &(node_normal[3][0]));
    //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );
  }
  else if( coordinates[0] <  lower_bound && coordinates[1] >= lower_bound &&
	   coordinates[2] >= lower_bound ){
    /* Region C
                        3            
                       /\
                      0  \
                       \  \
                        \  \
                         \  \
                          \C \
                           \  \
                            \  \
                             1--2
    */
    if( fabs(GetEdgeCurvature(1)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    } 
    // Compute the data for node 0
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    coords[2] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute the data for node 2
    coords[0] = 0.0;
    coords[1] = 1.0-lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[2][0]) );
    GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute the data for node 3
    coords[0] = 0.0;
    coords[1] = lower_bound;
    coords[2] = 1.0-lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    GetEdgeSmoothedNormal(1, &(node_normal[3][0]));
    //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );

  }
  else if( coordinates[0] >= lower_bound && coordinates[1] <  lower_bound &&
	   coordinates[2] >= lower_bound ){
    /* Region D
                       2
                       /\
                      /  1
                     /  /
                    /  /
                   /  /
                  / D/
                 /  /
                /  /
               3--0
    */
    if( fabs(GetEdgeCurvature(2)) > critical_curvature) {
      // No smoothing is needed
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    } 
    // Compute the data for node 0
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    coords[2] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    node_normal[1][0] = face_normal[0];
    node_normal[1][1] = face_normal[1];
    node_normal[1][2] = face_normal[2];
    // Compute the data for node 2
    coords[0] = lower_bound;
    coords[1] = 0.0;
    coords[2] = 1.0-lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[2][0]) );
    GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[2][0]), &(node_normal[2][0]) );
    // Compute the data for node 3
    coords[0] = 1.0-lower_bound;
    coords[1] = 0.0;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    GetEdgeSmoothedNormal(2, &(node_normal[3][0]));
    //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		    &(node_position[3][0]), &(node_normal[3][0]) );

  }
  else if( coordinates[0] > upper_bound && coordinates[1] < lower_bound &&
	   coordinates[2] < lower_bound ){
    /* Region E
      
               1--0
              / E/
             2--3
    */
    curvature0 = GetEdgeCurvature(2);
    curvature1 = GetEdgeCurvature(0);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( fabs(curvature0) > critical_curvature) smooth_edge0 = false;
    if( fabs(curvature1) > critical_curvature) smooth_edge1 = false;
    if( !smooth_edge0 && !smooth_edge1 ){
      // No smoothing is needed along either edge, so return the face normal
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute the data for node 0
    coords[0] = upper_bound;
    coords[1] = lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = 1.0-lower_bound;
    coords[1] = 0.0;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    if( smooth_edge0 ) {
      GetEdgeSmoothedNormal(2, &(node_normal[1][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }
    // Compute the data for node 2
    node_position[2][0] = Node(0)->Variable(CURRENT_POSITION)[0];
    node_position[2][1] = Node(0)->Variable(CURRENT_POSITION)[1];
    node_position[2][2] = Node(0)->Variable(CURRENT_POSITION)[2];
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0] = Node(0)->Variable(NODE_NORMAL)[0];
      node_normal[2][1] = Node(0)->Variable(NODE_NORMAL)[1];
      node_normal[2][2] = Node(0)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } 
    // Compute the data for node 3
    coords[0] = 1.0-lower_bound;
    coords[1] = lower_bound;
    coords[2] = 0.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    if( smooth_edge1 ) {
      GetEdgeSmoothedNormal(0, &(node_normal[3][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
  }
  else if( coordinates[0] < lower_bound && coordinates[1] > upper_bound &&
	   coordinates[2] < lower_bound ){
    /* Region F
       
             0--3
              \F \
               1--2
    */
    curvature0 = GetEdgeCurvature(0);
    curvature1 = GetEdgeCurvature(1);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( fabs(curvature0) > critical_curvature) smooth_edge0 = false;
    if( fabs(curvature1) > critical_curvature) smooth_edge1 = false;
    if( !smooth_edge0 && !smooth_edge1 ){
      // No smoothing is needed along either edge, so return the face normal
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Compute the data for node 0
    coords[0] = lower_bound;
    coords[1] = upper_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = lower_bound;
    coords[1] = 1.0-lower_bound;
    coords[2] = 0.0;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    if( smooth_edge0 ) {
      GetEdgeSmoothedNormal(0, &(node_normal[1][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }
    // Compute the data for node 2
    node_position[2][0] = Node(1)->Variable(CURRENT_POSITION)[0];
    node_position[2][1] = Node(1)->Variable(CURRENT_POSITION)[1];
    node_position[2][2] = Node(1)->Variable(CURRENT_POSITION)[2];
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0] = Node(1)->Variable(NODE_NORMAL)[0];
      node_normal[2][1] = Node(1)->Variable(NODE_NORMAL)[1];
      node_normal[2][2] = Node(1)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      GetEdgeSmoothedNormal(0, &(node_normal[2][0]));
      //Edge(0)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    }
    // Compute the data for node 3
    coords[0] = 0.0;
    coords[1] = 1.0-lower_bound;
    coords[2] = lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    if( smooth_edge1 ) {
      GetEdgeSmoothedNormal(1, &(node_normal[3][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
  }
  else if( coordinates[0] < lower_bound && coordinates[1] < lower_bound &&
	   coordinates[2] > upper_bound ){
    /* Region G
       
            2
          3/ \1
           \G/
            0
    */
    curvature0 = GetEdgeCurvature(1);
    curvature1 = GetEdgeCurvature(2);
    smooth_edge0 = true;
    smooth_edge1 = true;
    if( fabs(curvature0) > critical_curvature) smooth_edge0 = false;
    if( fabs(curvature1) > critical_curvature) smooth_edge1 = false;
    if( !smooth_edge0 && !smooth_edge1 ){
      // No smoothing is needed along either edge, so return the face normal
      smooth_normal[0] = face_normal[0];
      smooth_normal[1] = face_normal[1];
      smooth_normal[2] = face_normal[2];
      return;
    }
    // Get the global coordinates of the contact point
    Compute_Global_Coordinates( CURRENT_POSITION, coordinates, contact_point );
    // Compute the data for node 0
    coords[0] = lower_bound;
    coords[1] = lower_bound;
    coords[2] = upper_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[0][0]) );
    node_normal[0][0] = face_normal[0];
    node_normal[0][1] = face_normal[1];
    node_normal[0][2] = face_normal[2];
    // Compute the data for node 1
    coords[0] = 0.0;
    coords[1] = lower_bound;
    coords[2] = 1.0-lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[1][0]) );
    if( smooth_edge0 ) {
      GetEdgeSmoothedNormal(1, &(node_normal[1][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[1][0]), &(node_normal[1][0]) );
    } else {
      node_normal[1][0] = face_normal[0];
      node_normal[1][1] = face_normal[1];
      node_normal[1][2] = face_normal[2];
    }
    // Compute the data for node 2
    node_position[2][0] = Node(2)->Variable(CURRENT_POSITION)[0];
    node_position[2][1] = Node(2)->Variable(CURRENT_POSITION)[1];
    node_position[2][2] = Node(2)->Variable(CURRENT_POSITION)[2];
    if( (smooth_edge0 && smooth_edge1) ||
	resolution == ContactSearch::USE_NODE_NORMAL){
      node_normal[2][0] = Node(2)->Variable(NODE_NORMAL)[0];
      node_normal[2][1] = Node(2)->Variable(NODE_NORMAL)[1];
      node_normal[2][2] = Node(2)->Variable(NODE_NORMAL)[2];
    } else if( smooth_edge0 ){
      GetEdgeSmoothedNormal(1, &(node_normal[2][0]));
      //Edge(1)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    } else if( smooth_edge1 ){
      GetEdgeSmoothedNormal(2, &(node_normal[2][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[2][0]), &(node_normal[2][0]) );
    }
    // Compute the data for node 3
    coords[0] = lower_bound;
    coords[1] = 0.0;
    coords[2] = 1.0-lower_bound;
    Compute_Global_Coordinates( CURRENT_POSITION, coords, &(node_position[3][0]) );
    if( smooth_edge1 ) {
      GetEdgeSmoothedNormal(2, &(node_normal[3][0]));
      //Edge(2)->Smooth_Normal( FACE_NORMAL,
	//		      &(node_position[3][0]), &(node_normal[3][0]) );
    } else {
      node_normal[3][0] = face_normal[0];
      node_normal[3][1] = face_normal[1];
      node_normal[3][2] = face_normal[2];
    }
  }
  // Get the global coordinates of the contact point
  Compute_Global_Coordinates( CURRENT_POSITION, coordinates, contact_point );

  // Compute Contact Point in sub area
  Real local_coords[3];
  ContactQuadFaceL4::Compute_Quad_Local_Coords( node_position,
					       contact_point,
					       local_coords);

  // Now interpolate to get the normal
  ContactQuadFaceL4::Interpolate_Vector( local_coords, 
					 node_normal, 
					 smooth_normal );

  // Now make smooth_normal a unit vector
  Normalize(smooth_normal);
}

void ContactTriFaceL3::Compute_Node_Areas( VariableHandle POSITION,
	                                   VariableHandle /*FACE_NORMAL*/,
                                           Real* node_areas )
{
  // Get the nodal positions
  Real* pos_0 = Node(0)->Variable(POSITION);
  Real* pos_1 = Node(1)->Variable(POSITION);
  Real* pos_2 = Node(2)->Variable(POSITION);

  // Construct vectors between the nodes
  Real vec_0_1[3], vec_0_2[3];
  vec_0_1[0] = pos_1[0] - pos_0[0];
  vec_0_1[1] = pos_1[1] - pos_0[1];
  vec_0_1[2] = pos_1[2] - pos_0[2];
  vec_0_2[0] = pos_2[0] - pos_0[0];
  vec_0_2[1] = pos_2[1] - pos_0[1];
  vec_0_2[2] = pos_2[2] - pos_0[2];

  // Take the cross product
  Real cross[3];
  Cross(vec_0_1, vec_0_2, cross);

  // The total area is then the 1/2 the magnitude of the vector.
  // The node area is then 1/3 of the total area.
  Real node_area = Magnitude(cross) / 6.0;

  node_areas[0] = node_area;
  node_areas[1] = node_area;
  node_areas[2] = node_area;
}


int ContactTriFaceL3::FaceEdge_Intersection(VariableHandle POSITION,
                                            ContactEdge* edge,Real* coords)
{
  int intersection=0;
#if 0
//#if CONTACT_DEBUG_PRINT_LEVEL>=1
  std::cerr << "ContactTriFaceL3::FaceEdge_Intersection not yet implemented\n";
//#endif
  POSTCONDITION( 0 );
  return intersection;
#else
  // Do bounding box check
  int i, j;
  Real edge_min[3], edge_max[3];
  Real face_min[3], face_max[3];
  Real* node_position = Node(0)->Variable(POSITION);
  for (j=0; j<3; ++j) {
    face_min[j] = node_position[j];
    face_max[j] = node_position[j];
  }
  for (i=1; i<Nodes_Per_Face(); ++i) {
    node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      face_min[j] = std::min(face_min[j],node_position[j]);
      face_max[j] = std::max(face_max[j],node_position[j]);
    }
  }
  node_position = edge->Node(0)->Variable(POSITION);
  for (j=0; j<3; ++j) {
    edge_min[j] = node_position[j];
    edge_max[j] = node_position[j];
  }
  node_position = edge->Node(1)->Variable(POSITION);
  for (j=0; j<3; ++j) {
    edge_min[j] = std::min(edge_min[j],node_position[j]);
    edge_max[j] = std::max(edge_max[j],node_position[j]);
  }
  if ( edge_min[0]>face_max[0] || 
       edge_min[1]>face_max[1] ||  
       edge_min[2]>face_max[2] ||
       edge_max[0]<face_min[0] || 
       edge_max[1]<face_min[1] || 
       edge_max[2]<face_min[2] ) return 0;
       
  // Compute edge ray
  Real  edge_pnt[3];
  Real  edge_dir[3];
  Real* edge_node_position0 = edge->Node(0)->Variable(POSITION);
  Real* edge_node_position1 = edge->Node(1)->Variable(POSITION);
  for (j=0; j<3; ++j) {
    edge_pnt[j] = edge_node_position0[j];
    edge_dir[j] = edge_node_position1[j]-edge_node_position0[j];
  }
  Normalize(edge_dir);
  Real tmax;
  if (std::fabs(edge_dir[0])>=std::fabs(edge_dir[1]) && std::fabs(edge_dir[0])>=std::fabs(edge_dir[2])) {
    tmax = (edge_node_position1[0]-edge_node_position0[0])/edge_dir[0];
  } else
  if (std::fabs(edge_dir[1])>=std::fabs(edge_dir[0]) && std::fabs(edge_dir[1])>=std::fabs(edge_dir[2])) {
    tmax = (edge_node_position1[1]-edge_node_position0[1])/edge_dir[1];
  } else
  if (std::fabs(edge_dir[2])>=std::fabs(edge_dir[0]) && std::fabs(edge_dir[2])>=std::fabs(edge_dir[1])) {
    tmax = (edge_node_position1[2]-edge_node_position0[2])/edge_dir[2];
  }
  tmax *= 1.05;
  
  Real* node_position0 = Node(0)->Variable(POSITION);
  Real* node_position1 = Node(1)->Variable(POSITION);
  Real* node_position2 = Node(2)->Variable(POSITION);

  // Compute the face normal
  Real dx1  = node_position1[0] - node_position0[0];
  Real dy1  = node_position1[1] - node_position0[1];
  Real dz1  = node_position1[2] - node_position0[2];
  Real dx2  = node_position2[0] - node_position0[0];
  Real dy2  = node_position2[1] - node_position0[1];
  Real dz2  = node_position2[2] - node_position0[2];

  Real normal[3];
  Cross(dx1,dy1,dz1, // vector a
              dx2,dy2,dz2, // vector b
	      normal);     // vector c = a X b;

  Normalize(normal);
  
  // determine (edge ray)/(face plane) intersection
  Real n_dot_d = Dot(normal, edge_dir);

  // if (edge ray) and (tri3 plane) are parallel => no intersection
  if (std::fabs(n_dot_d)<1.0e-10) {
    return 0;
  }

  Real q[3] = {node_position0[0]-edge_pnt[0],
               node_position0[1]-edge_pnt[1],
               node_position0[2]-edge_pnt[2]};

  Real t = Dot(normal,q)/n_dot_d;

  if (t<0.0 || t>tmax) {
    return 0;
  }

  Real P[3] = { edge_pnt[0] + edge_dir[0]*t,
                edge_pnt[1] + edge_dir[1]*t,
                edge_pnt[2] + edge_dir[2]*t};

  // now determine if this point is inside the face
  Real u0, u1, u2, v0, v1, v2;
  Real nx = normal[0]*normal[0];
  Real ny = normal[1]*normal[1];
  Real nz = normal[2]*normal[2];
  if (nx>=ny && nx>=nz) {
    u0 = P[1] - node_position0[1];
    v0 = P[2] - node_position0[2];
    u1 = node_position1[1] - node_position0[1];
    u2 = node_position2[1] - node_position0[1];
    v1 = node_position1[2] - node_position0[2];
    v2 = node_position2[2] - node_position0[2];
  }
  else if (ny>=nx && ny>=nz) { 
    u0 = P[2] - node_position0[2];
    v0 = P[0] - node_position0[0];
    u1 = node_position1[2] - node_position0[2];
    u2 = node_position2[2] - node_position0[2];
    v1 = node_position1[0] - node_position0[0];
    v2 = node_position2[0] - node_position0[0];
  }
  else if (nz>=nx && nz>=ny) {
    u0 = P[0] - node_position0[0];
    v0 = P[1] - node_position0[1];
    u1 = node_position1[0] - node_position0[0];
    u2 = node_position2[0] - node_position0[0];
    v1 = node_position1[1] - node_position0[1];
    v2 = node_position2[1] - node_position0[1];
  }
  /*==========================================*/
  /* BEGIN BARYCENTRIC INTERSECTION ALGORITHM */
  /*==========================================*/
  Real alpha, beta;
  if (u1!=0.0)    {  /* common case */
      beta = (v0*u1 - u0*v1)/(v2*u1 - u2*v1);
      if ((beta>=0.0) && (beta<=1.0)) {
          alpha        = (u0 - beta*u2)/u1;
          intersection = ((alpha>=0.0) && ((alpha+beta)<=1.0));
      }
  } else {           /* uncommon case */
      beta = u0/u2;
      if ((beta>=0.0) && (beta<=1.0)) {
          alpha        = (v0 - beta*v2)/v1;
          intersection = ((alpha>=0.0) && ((alpha+beta)<=1.0));
      }
  }
  if (intersection) {
    coords[0] = P[0];
    coords[1] = P[1];
    coords[2] = P[2];
  }
  return (intersection);
#endif
}

void ContactTriFaceL3::Evaluate_Shape_Functions( Real* local_coords,
						 Real* shape_functions )
{
  
  PRECONDITION( std::fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  Compute_Shape_Functions( local_coords, shape_functions );
}

void ContactTriFaceL3::Compute_Global_Coordinates( VariableHandle POSITION,
						   Real* local_coords,
						   Real* global_coords )
{
  Real node_positions[3][3];
  for(int i=0; i<3; ++i ){
    Real* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void ContactTriFaceL3::Compute_Local_Coordinates( Real Config_Param,
						  VariableHandle POSITION0,
						  VariableHandle POSITION1,
						  VariableHandle FACE_NORMAL,
						  Real* global_coords,
						  Real* local_coords )
{
  int  i, j;
  Real node_positions[3][3];
  if (Config_Param == 0.0) {
    for (i=0; i<Nodes_Per_Face(); ++i) {
      Real* node_position = Node(i)->Variable(POSITION0);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else if (Config_Param == 1.0) {
    for (i=0; i<Nodes_Per_Face(); ++i) {
      Real* node_position = Node(i)->Variable(POSITION1);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else {
    Real alpha = 1.0 - Config_Param, beta = Config_Param;
    for (i=0; i<Nodes_Per_Face(); ++i) {
      Real* node_position0 = Node(i)->Variable(POSITION0);
      Real* node_position1 = Node(i)->Variable(POSITION1);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = alpha*node_position0[j]+beta*node_position1[j];
      }
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

void ContactTriFaceL3::Compute_Local_Coordinates( VariableHandle POSITION,
						  Real* global_coords,
						  Real* local_coords )
{
  int  i, j;
  Real node_positions[3][3];
  for (i=0; i<Nodes_Per_Face(); ++i) {
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
/* on a linear T3 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary T3 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

void ContactTriFaceL3::Compute_Shape_Functions( Real local_coords[3],
						Real shape_functions[3] )
{
  
  PRECONDITION( std::fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  shape_functions[0] = local_coords[0];
  shape_functions[1] = local_coords[1];
  shape_functions[2] = local_coords[2];
}

void ContactTriFaceL3::Compute_Shape_Derivatives( Real local_coords[3],
					          Real shape_derivs[2][3] )
{
  
  PRECONDITION( std::fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  shape_derivs[0][0] =  1.0;
  shape_derivs[0][1] =  0.0;
  shape_derivs[0][2] = -1.0;
  shape_derivs[1][0] =  0.0;
  shape_derivs[1][1] =  1.0;
  shape_derivs[1][2] = -1.0;
}

void ContactTriFaceL3::Compute_Local_Coords( Real node_positions[MAX_NODES_PER_FACE][3],
					     Real global_coords[3],
					     Real local_coords[3] )
{
  int  i;
  Real spatial_tolerance = 1.0e-10;
  //
  // check for coincidence with one of the face nodes
  //
  // KHP: Couldn't we avoid the sqrt here by squaring both
  // sides of the spatial_tolerance comparison?
  //
  for (i=0; i<3; ++i) {
    Real dx = node_positions[i][0]-global_coords[0];
    Real dy = node_positions[i][1]-global_coords[1];
    Real dz = node_positions[i][2]-global_coords[2];
    Real d  = std::sqrt(dx*dx+dy*dy+dz*dz);
    if (d<spatial_tolerance) break;
  }
  switch (i) {
  case 0:
    local_coords[0] = 1.0;
    local_coords[1] = 0.0;
    local_coords[2] = 0.0;
    break;
  case 1:
    local_coords[0] = 0.0;
    local_coords[1] = 1.0;
    local_coords[2] = 0.0;
    break;
  case 2:
    local_coords[0] = 0.0;
    local_coords[1] = 0.0;
    local_coords[2] = 1.0;
    break;
  }
  if (i<3) return;

  /*
                      2
                     /|\
                    / | \
                   /a1|a0\
                  /  /c\  \
                 / /     \ \
                //    a2   \\ 
               0-------------1
  */

  // Compute vectors from the point to each node
  Real v_c_0[3] = { node_positions[0][0] - global_coords[0],
                    node_positions[0][1] - global_coords[1],
                    node_positions[0][2] - global_coords[2] };

  Real v_c_1[3] = { node_positions[1][0] - global_coords[0],
                    node_positions[1][1] - global_coords[1],
                    node_positions[1][2] - global_coords[2] };

  Real v_c_2[3] = { node_positions[2][0] - global_coords[0],
                    node_positions[2][1] - global_coords[1],
                    node_positions[2][2] - global_coords[2] };

  // Compute the normal (as v_0_1 x v_0_2)
  Real v_0_1[3] = { node_positions[1][0] - node_positions[0][0],
                    node_positions[1][1] - node_positions[0][1],
                    node_positions[1][2] - node_positions[0][2] };
  Real v_0_2[3] = { node_positions[2][0] - node_positions[0][0],
                    node_positions[2][1] - node_positions[0][1],
                    node_positions[2][2] - node_positions[0][2] };

  Real Normal[3];
  Cross(v_0_1, v_0_2, Normal);

  // Compute the areas
  Real a0 = ScalarTripleProduct(v_c_1, v_c_2, Normal);
  Real a1 = ScalarTripleProduct(v_c_2, v_c_0, Normal);
  Real a2 = ScalarTripleProduct(v_c_0, v_c_1, Normal);

  // Compute the local coordinates
  Real area = a0 + a1 + a2;
  Real L1   = a0/area;
  Real L2   = a1/area;
  // If it's close to any of the edges, snap to it
  if (L1<1.0+spatial_tolerance) {
    L1 = std::min(L1, 1.0);
  }
  if (L1>-spatial_tolerance) {
    L1 = std::max(L1, 0.0);
  }
  if (L2<1.0+spatial_tolerance) {
    L2 = std::min(L2, 1.0);
  }
  if (L2>-spatial_tolerance) {
    L2 = std::max(L2, 0.0);
  }
  local_coords[0] = L1;
  local_coords[1] = L2;
  local_coords[2] = 1.0-L1-L2;
}

void ContactTriFaceL3::Compute_Global_Coords( Real node_positions[3][3],
					      Real local_coords[3],
					      Real global_coords[3] )
{
  Real N[3];
  int  nnodes=3;
  global_coords[0] = 0.0;
  global_coords[1] = 0.0;
  global_coords[2] = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; ++i ){
    for (int j=0; j<3; ++j) {
      global_coords[j] += N[i]*node_positions[i][j];
    }
  }
}

void ContactTriFaceL3::Interpolate_Scalar( Real  local_coords[3],
					   Real  node_scalars[3],
					   Real& interpolated_scalar )
{
  Real N[3];
  int  nnodes=3;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

void ContactTriFaceL3::Interpolate_Vector( Real local_coords[3],
					   Real node_vectors[3][3],
					   Real interpolated_vector[3] )
{
  Real N[3];
  int  nnodes=3;
  interpolated_vector[0] = 0.0;
  interpolated_vector[1] = 0.0;
  interpolated_vector[2] = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; ++i ){
    for (int j=0; j<3; ++j) {
      interpolated_vector[j] += N[i]*node_vectors[i][j];
    }
  }
}
