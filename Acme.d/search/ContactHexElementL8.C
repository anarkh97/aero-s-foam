// $Id$

#ifndef ContactHexElementL8_C_
#define ContactHexElementL8_C_

#include "allocators.h"
#include "ContactNode.h"
#include "ContactLineEdgeL2.h"
#include "ContactQuadFaceL4.h"
#include "ContactHexElementL8.h"
#include "ContactFixedSizeAllocator.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <new>

template<typename DataType>
ContactHexElemL8<DataType>::ContactHexElemL8( int Block_Index, 
				    int Host_Index_in_Block, int key ) 
  : ContactElem<DataType>( ContactSearch::HEXELEML8,Block_Index,Host_Index_in_Block,key) 
{}

template<typename DataType>
ContactHexElemL8<DataType>* ContactHexElemL8<DataType>::new_ContactHexElemL8(
                        ContactFixedSizeAllocator& alloc,
                        int Block_Index, int Host_Index_in_Block, int key)
{
  return new (alloc.New_Frag())
             ContactHexElemL8<DataType>(Block_Index, Host_Index_in_Block, key);
}

template<typename DataType>
void ContactHexElemL8_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactHexElemL8<DataType>),
                100,  // block size
                0,   // initial block size
                sizeof(DataType) );
  alloc.Set_Name( "ContactHexElemL8<DataType> allocator" );
}

template<typename DataType>
ContactHexElemL8<DataType>::~ContactHexElemL8() {}

template<typename DataType>
void ContactHexElemL8<DataType>::BuildTopology(int nID, int eID, int fID,
				     ContactFixedSizeAllocator* allocators)
{
  int i;
  ContactNode<DataType>* node;
  ContactEdge<DataType>* edge;
  ContactFace<DataType>* face;
  
  int NextID = nID;
  for( i=0 ; i<Nodes_Per_Element() ; ++i ) {
    node = ContactNode<DataType>::new_ContactNode(allocators,
                                        ContactSearch::NODE, 
                                        ++NextID );
    ConnectNode(i, node);
  }
  NextID = eID;
  for( i=0 ; i<Edges_Per_Element() ; ++i ) {
    edge = ContactLineEdgeL2<DataType>::new_ContactLineEdgeL2( 
                        allocators[ContactSearch::ALLOC_ContactLineEdgeL2],
                        ContactSearch::LINEEDGEL2, ++NextID);
    ConnectEdge(i, edge);
  }
  edge = Edge(0);
  edge->ConnectNode(0, Node(0));
  edge->ConnectNode(1, Node(1));
  edge = Edge(1);
  edge->ConnectNode(0, Node(1));
  edge->ConnectNode(1, Node(5));
  edge = Edge(2);
  edge->ConnectNode(0, Node(5));
  edge->ConnectNode(1, Node(4));
  edge = Edge(3);
  edge->ConnectNode(0, Node(4));
  edge->ConnectNode(1, Node(0));
  edge = Edge(4);
  edge->ConnectNode(0, Node(1));
  edge->ConnectNode(1, Node(2));
  edge = Edge(5);
  edge->ConnectNode(0, Node(2));
  edge->ConnectNode(1, Node(6));
  edge = Edge(6);
  edge->ConnectNode(0, Node(6));
  edge->ConnectNode(1, Node(5));
  edge = Edge(7);
  edge->ConnectNode(0, Node(2));
  edge->ConnectNode(1, Node(3));
  edge = Edge(8);
  edge->ConnectNode(0, Node(3));
  edge->ConnectNode(1, Node(7));
  edge = Edge(9);
  edge->ConnectNode(0, Node(7));
  edge->ConnectNode(1, Node(6));
  edge = Edge(10);
  edge->ConnectNode(0, Node(4));
  edge->ConnectNode(1, Node(7));
  edge = Edge(11);
  edge->ConnectNode(0, Node(3));
  edge->ConnectNode(1, Node(0));
  
  NextID = fID;
  for( i=0 ; i<Faces_Per_Element() ; ++i ) {
    face = ContactQuadFaceL4<DataType>::new_ContactQuadFaceL4(allocators, ++NextID );
    ConnectFace(i, face);
  }
  face = Face(0);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(1));
  face->ConnectNode(2, Node(5));
  face->ConnectNode(3, Node(4));
  face = Face(1);
  face->ConnectNode(0, Node(1));
  face->ConnectNode(1, Node(2));
  face->ConnectNode(2, Node(6));
  face->ConnectNode(3, Node(5));
  face = Face(2);
  face->ConnectNode(0, Node(2));
  face->ConnectNode(1, Node(3));
  face->ConnectNode(2, Node(7));
  face->ConnectNode(3, Node(6));
  face = Face(3);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(4));
  face->ConnectNode(2, Node(7));
  face->ConnectNode(3, Node(3));
  face = Face(4);
  face->ConnectNode(0, Node(0));
  face->ConnectNode(1, Node(3));
  face->ConnectNode(2, Node(2));
  face->ConnectNode(3, Node(1));
  face = Face(5);
  face->ConnectNode(0, Node(4));
  face->ConnectNode(1, Node(5));
  face->ConnectNode(2, Node(6));
  face->ConnectNode(3, Node(7));
}

template<typename DataType>
void ContactHexElemL8<DataType>::DeleteTopology(ContactFixedSizeAllocator* allocators)
{
  int i;
  for( i=0 ; i<Nodes_Per_Element() ; ++i ) {
    ContactNode<DataType>* node = Node(i);
    node->~ContactNode<DataType>();
    allocators[ContactSearch::ALLOC_ContactNode].Delete_Frag(node);
  }
  for( i=0 ; i<Edges_Per_Element() ; ++i ) {
    ContactEdge<DataType>* edge = Edge(i);
    edge->~ContactEdge<DataType>();
    allocators[ContactSearch::ALLOC_ContactLineEdgeL2].Delete_Frag(edge);
  }
  for( i=0 ; i<Faces_Per_Element() ; ++i ) {
    ContactFace<DataType>* face = Face(i);
    face->~ContactFace<DataType>();
    allocators[ContactSearch::ALLOC_ContactQuadFaceL4].Delete_Frag(face);
  }
}

template<typename DataType>
void ContactHexElemL8<DataType>::UpdateTopology(ContactFace<DataType>* face, 
				      VariableHandle POSITION,
				      VariableHandle FACE_NORMAL,
				      VariableHandle NODE_NORMAL,
				      Real tol, bool use_node_normals)
{
  int i;
  int num_nodes = face->Nodes_Per_Face();
  for( i=0 ; i<num_nodes ; ++i ){
    DataType* projection;
    ContactNode<DataType>* face_node    = face->Node(i);
    ContactNode<DataType>* elem_node1   = Node(i);
    ContactNode<DataType>* elem_node2   = Node(i+num_nodes);
    DataType* face_node_position  = face_node->Variable(POSITION);
    DataType* elem_node_position1 = elem_node1->Variable(POSITION);
    DataType* elem_node_position2 = elem_node2->Variable(POSITION);
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

template<typename DataType>
bool
ContactHexElemL8<DataType>::Is_Local_Coordinates_Inside_Element( DataType* local_coords )
{
  if( local_coords[0] >= -1.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >= -1.0 && local_coords[1] <= 1.0 &&
      local_coords[2] >= -1.0 && local_coords[2] <= 1.0 )
    return true;
  return false;
}

template<typename DataType>
bool
ContactHexElemL8<DataType>::Is_Local_Coordinates_Near_Element( DataType* local_coords, DataType tolerance )
{
  DataType low_coord  = -(1.+tolerance);
  DataType high_coord = 1.+tolerance;
  if( local_coords[0] >= low_coord && local_coords[0] <= high_coord &&
      local_coords[1] >= low_coord && local_coords[1] <= high_coord &&
      local_coords[2] >= low_coord && local_coords[2] <= high_coord )
    return true;
  return false;
}

template<typename DataType>
void ContactHexElemL8<DataType>::Evaluate_Shape_Functions( DataType* local_coords,
						  DataType* shape_functions )
{
  Compute_Shape_Functions(local_coords, shape_functions);
}

template<typename DataType>
void ContactHexElemL8<DataType>::Compute_Global_Coordinates( VariableHandle POSITION,
						   DataType* local_coords,
						   DataType* global_coords )
{
  DataType node_positions[8][3];
  for(int i=0; i<Nodes_Per_Element(); ++i ){
    DataType* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

template<typename DataType>
void ContactHexElemL8<DataType>::Compute_Local_Coordinates( DataType Config_Param,
						  VariableHandle POSITION0, 
						  VariableHandle POSITION1, 
						  VariableHandle FACE_NORMAL,
						  DataType* global_coords,
						  DataType* local_coords )
{
  int i, j;
  DataType node_positions[8][3];
  if (Config_Param == 0.0) {
    for (i=0; i<Nodes_Per_Element(); ++i) {
      DataType* node_position = Node(i)->Variable(POSITION0);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else if (Config_Param == 1.0) {
    for (i=0; i<Nodes_Per_Element(); ++i) {
      DataType* node_position = Node(i)->Variable(POSITION1);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = node_position[j];
      }
    }
  } else {
    DataType alpha = 1.0 - Config_Param, beta = Config_Param;
    for (i=0; i<Nodes_Per_Element(); ++i) {
      DataType* node_position0 = Node(i)->Variable(POSITION0);
      DataType* node_position1 = Node(i)->Variable(POSITION1);
      for (j=0; j<3; ++j) {
        node_positions[i][j] = alpha*node_position0[j]+beta*node_position1[j];
      }
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

template<typename DataType>
void ContactHexElemL8<DataType>::Compute_Local_Coordinates( VariableHandle POSITION,
						  DataType* global_coords,
						  DataType* local_coords )
{
  int i, j;
  DataType node_positions[8][3];
  for (i=0; i<Nodes_Per_Element(); ++i) {
    DataType* node_position = Node(i)->Variable(POSITION);
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
/* on a linear H8 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary H8 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

template<typename DataType>
void ContactHexElemL8<DataType>::Compute_Shape_Functions( DataType* local_coords,
						DataType* shape_functions )
{
  shape_functions[0] = 0.125*(1.0-local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[1] = 0.125*(1.0+local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[2] = 0.125*(1.0+local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[3] = 0.125*(1.0-local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[4] = 0.125*(1.0-local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[5] = 0.125*(1.0+local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[6] = 0.125*(1.0+local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[7] = 0.125*(1.0-local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0+local_coords[2]);
}

template<typename DataType>
void ContactHexElemL8<DataType>::Compute_Shape_Derivatives( DataType* local_coords,
						  DataType shape_derivs[3][8] )
{
  shape_derivs[0][0] = -0.125*(1.0-local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][1] =  0.125*(1.0-local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][2] =  0.125*(1.0+local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][3] = -0.125*(1.0+local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][4] = -0.125*(1.0-local_coords[1])*(1.0+local_coords[2]);
  shape_derivs[0][5] =  0.125*(1.0-local_coords[1])*(1.0+local_coords[2]);
  shape_derivs[0][6] =  0.125*(1.0+local_coords[1])*(1.0+local_coords[2]);
  shape_derivs[0][7] = -0.125*(1.0+local_coords[1])*(1.0+local_coords[2]);

  shape_derivs[1][0] = -0.125*(1.0-local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][1] = -0.125*(1.0+local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][2] =  0.125*(1.0+local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][3] =  0.125*(1.0-local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][4] = -0.125*(1.0-local_coords[0])*(1.0+local_coords[2]);
  shape_derivs[1][5] = -0.125*(1.0+local_coords[0])*(1.0+local_coords[2]);
  shape_derivs[1][6] =  0.125*(1.0+local_coords[0])*(1.0+local_coords[2]);
  shape_derivs[1][7] =  0.125*(1.0-local_coords[0])*(1.0+local_coords[2]);

  shape_derivs[2][0] = -0.125*(1.0-local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][1] = -0.125*(1.0+local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][2] = -0.125*(1.0+local_coords[0])*(1.0+local_coords[1]);
  shape_derivs[2][3] = -0.125*(1.0-local_coords[0])*(1.0+local_coords[1]);
  shape_derivs[2][4] =  0.125*(1.0-local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][5] =  0.125*(1.0+local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][6] =  0.125*(1.0+local_coords[0])*(1.0+local_coords[1]);
  shape_derivs[2][7] =  0.125*(1.0-local_coords[0])*(1.0+local_coords[1]);
}

template<typename DataType>
void 
ContactHexElemL8<DataType>::Compute_Local_Coords( DataType node_positions[8][3], 
					DataType global_coords[3],
					DataType local_coords[3] )
{
  using std::sqrt;
  using std::abs;
  using std::min;
  using std::max;
  //
  // 1st check for coincidence with one of the face nodes
  //
  int  i, j;
  int  nnodes=8;
  DataType spatial_tolerance = 1.0e-10;

  // are we on a node?
  for (i=0; i<nnodes; ++i) {
    DataType dx = node_positions[i][0]-global_coords[0];
    DataType dy = node_positions[i][1]-global_coords[1];
    DataType dz = node_positions[i][2]-global_coords[2];
    DataType d  = sqrt(dx*dx+dy*dy+dz*dz);
    if (d<spatial_tolerance) break;
  }
  switch (i) {
  case 0:
    local_coords[0] = -1.0;
    local_coords[1] = -1.0;
    local_coords[2] = -1.0;
    break;
  case 1:
    local_coords[0] =  1.0;
    local_coords[1] = -1.0;
    local_coords[2] = -1.0;
    break;
  case 2:
    local_coords[0] =  1.0;
    local_coords[1] =  1.0;
    local_coords[2] = -1.0;
    break;
  case 3:
    local_coords[0] = -1.0;
    local_coords[1] =  1.0;
    local_coords[2] = -1.0;
    break;
  case 4:
    local_coords[0] = -1.0;
    local_coords[1] = -1.0;
    local_coords[2] =  1.0;
    break;
  case 5:
    local_coords[0] =  1.0;
    local_coords[1] = -1.0;
    local_coords[2] =  1.0;
    break;
  case 6:
    local_coords[0] =  1.0;
    local_coords[1] =  1.0;
    local_coords[2] =  1.0;
    break;
  case 7:
    local_coords[0] = -1.0;
    local_coords[1] =  1.0;
    local_coords[2] =  1.0;
    break;
  }
  if (i<nnodes) return;
  //
  // else use newton's method to iterate
  //
  int  iterations=0;
  int  max_iterations=400;
  bool converged = false;
  DataType tolerance = 1.0e-12;
  DataType u, u0=0.0, u1, du;
  DataType v, v0=0.0, v1, dv;
  DataType w, w0=0.0, w1, dw;
  DataType f[3], J[3][3], invJ[3][3];
  DataType shape_derivatives[3][8];
  while (!converged && iterations<max_iterations) {
    local_coords[0] = u0;
    local_coords[1] = v0;
    local_coords[2] = w0;
    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives( local_coords, shape_derivatives );
    for (i=0; i<3; ++i) {
      J[0][i] = 0.0;
      J[1][i] = 0.0;
      J[2][i] = 0.0;
      for (j=0; j<8; ++j) {
        J[0][i] += shape_derivatives[i][j]*node_positions[j][0];
        J[1][i] += shape_derivatives[i][j]*node_positions[j][1];
        J[2][i] += shape_derivatives[i][j]*node_positions[j][2];
      }
    }
    
    DataType detJ  =  1.0/(J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-
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

    Compute_Global_Coords( node_positions, local_coords, f );
    u  = f[0]-global_coords[0];
    v  = f[1]-global_coords[1];
    w  = f[2]-global_coords[2];
    u1 = u0-(invJ[0][0]*u+invJ[0][1]*v+invJ[0][2]*w);
    v1 = v0-(invJ[1][0]*u+invJ[1][1]*v+invJ[1][2]*w);
    w1 = w0-(invJ[2][0]*u+invJ[2][1]*v+invJ[2][2]*w);
    du = abs(u1-u0);
    dv = abs(v1-v0);
    dw = abs(w1-w0);
    u0 = u1;
    v0 = v1;
    w0 = w1;
    if (du<tolerance && dv<tolerance && dw<tolerance) converged = true;
    ++iterations;
  }
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  if (!converged) {
    std::cerr << "ContactHexElemL8<DataType>::Compute_Local_Coordinates() did not converge" 
	 << std::endl;
    std::cerr << "                     after "<<max_iterations
         << " iterations:  du = "<<du
         <<";  dv = "<<dv
         <<";  dw = "<<dw<<std::endl;
  }
#endif
  POSTCONDITION(converged);
  // If it's close to any of the edges, snap to it
  if (abs(u0)<1.0+spatial_tolerance) {
    u0 = min(u0, 1.0);
    u0 = max(u0,-1.0);
  }
  if (abs(v0)<1.0+spatial_tolerance) {
    v0 = min(v0, 1.0);
    v0 = max(v0,-1.0);
  }
  if (abs(w0)<1.0+spatial_tolerance) {
    w0 = min(w0, 1.0);
    w0 = max(w0,-1.0);
  }
  local_coords[0] = u0;
  local_coords[1] = v0;
  local_coords[2] = w0;
}

template<>
inline void
ContactHexElemL8<ActiveScalar>::Compute_Local_Coords( ActiveScalar active_node_positions[8][3],
                                                      ActiveScalar active_global_coords[3],
                                                      ActiveScalar active_local_coords[3] )
{
  int  i, j;
  const int  nnodes=8;

  double node_positions[nnodes][3], global_coords[3], local_coords[3];
  for(i=0; i<3; ++i) {
    global_coords[i] = active_global_coords[i].value();
    for(j=0; j<nnodes; ++j)
      node_positions[j][i] = active_node_positions[j][i].value();
  }

  using std::sqrt;
  using std::abs;
  using std::min;
  using std::max;

  double spatial_tolerance = 1.0e-10;

  // are we on a node?
  for (i=0; i<nnodes; ++i) {
    double dx = node_positions[i][0]-global_coords[0];
    double dy = node_positions[i][1]-global_coords[1];
    double dz = node_positions[i][2]-global_coords[2];
    double d  = sqrt(dx*dx+dy*dy+dz*dz);
    if (d<spatial_tolerance) break;
  }
  switch (i) {
  case 0:
    local_coords[0] = -1.0;
    local_coords[1] = -1.0;
    local_coords[2] = -1.0;
    break;
  case 1:
    local_coords[0] =  1.0;
    local_coords[1] = -1.0;
    local_coords[2] = -1.0;
    break;
  case 2:
    local_coords[0] =  1.0;
    local_coords[1] =  1.0;
    local_coords[2] = -1.0;
    break;
  case 3:
    local_coords[0] = -1.0;
    local_coords[1] =  1.0;
    local_coords[2] = -1.0;
    break;
  case 4:
    local_coords[0] = -1.0;
    local_coords[1] = -1.0;
    local_coords[2] =  1.0;
    break;
  case 5:
    local_coords[0] =  1.0;
    local_coords[1] = -1.0;
    local_coords[2] =  1.0;
    break;
  case 6:
    local_coords[0] =  1.0;
    local_coords[1] =  1.0;
    local_coords[2] =  1.0;
    break;
  case 7:
    local_coords[0] = -1.0;
    local_coords[1] =  1.0;
    local_coords[2] =  1.0;
    break;
  }
  if (i>=nnodes) {
    //
    // use newton's method to iterate (values only)
    //
    int  iterations=0;
    int  max_iterations=400;
    bool converged = false;
    double tolerance = 1.0e-12;
    double u, u0=0.0, u1, du;
    double v, v0=0.0, v1, dv;
    double w, w0=0.0, w1, dw;
    double f[3], J[3][3], invJ[3][3];
    double shape_derivatives[3][8], shape_functions[8];;
    while (!converged && iterations<max_iterations) {
      local_coords[0] = u0;
      local_coords[1] = v0;
      local_coords[2] = w0;

      // BUILD JACOBIAN AND INVERT
      shape_derivatives[0][0] = -0.125*(1.0-local_coords[1])*(1.0-local_coords[2]);
      shape_derivatives[0][1] =  0.125*(1.0-local_coords[1])*(1.0-local_coords[2]);
      shape_derivatives[0][2] =  0.125*(1.0+local_coords[1])*(1.0-local_coords[2]);
      shape_derivatives[0][3] = -0.125*(1.0+local_coords[1])*(1.0-local_coords[2]);
      shape_derivatives[0][4] = -0.125*(1.0-local_coords[1])*(1.0+local_coords[2]);
      shape_derivatives[0][5] =  0.125*(1.0-local_coords[1])*(1.0+local_coords[2]);
      shape_derivatives[0][6] =  0.125*(1.0+local_coords[1])*(1.0+local_coords[2]);
      shape_derivatives[0][7] = -0.125*(1.0+local_coords[1])*(1.0+local_coords[2]);

      shape_derivatives[1][0] = -0.125*(1.0-local_coords[0])*(1.0-local_coords[2]);
      shape_derivatives[1][1] = -0.125*(1.0+local_coords[0])*(1.0-local_coords[2]);
      shape_derivatives[1][2] =  0.125*(1.0+local_coords[0])*(1.0-local_coords[2]);
      shape_derivatives[1][3] =  0.125*(1.0-local_coords[0])*(1.0-local_coords[2]);
      shape_derivatives[1][4] = -0.125*(1.0-local_coords[0])*(1.0+local_coords[2]);
      shape_derivatives[1][5] = -0.125*(1.0+local_coords[0])*(1.0+local_coords[2]);
      shape_derivatives[1][6] =  0.125*(1.0+local_coords[0])*(1.0+local_coords[2]);
      shape_derivatives[1][7] =  0.125*(1.0-local_coords[0])*(1.0+local_coords[2]);

      shape_derivatives[2][0] = -0.125*(1.0-local_coords[0])*(1.0-local_coords[1]);
      shape_derivatives[2][1] = -0.125*(1.0+local_coords[0])*(1.0-local_coords[1]);
      shape_derivatives[2][2] = -0.125*(1.0+local_coords[0])*(1.0+local_coords[1]);
      shape_derivatives[2][3] = -0.125*(1.0-local_coords[0])*(1.0+local_coords[1]);
      shape_derivatives[2][4] =  0.125*(1.0-local_coords[0])*(1.0-local_coords[1]);
      shape_derivatives[2][5] =  0.125*(1.0+local_coords[0])*(1.0-local_coords[1]);
      shape_derivatives[2][6] =  0.125*(1.0+local_coords[0])*(1.0+local_coords[1]);
      shape_derivatives[2][7] =  0.125*(1.0-local_coords[0])*(1.0+local_coords[1]);

      for (i=0; i<3; ++i) {
        J[0][i] = 0.0;
        J[1][i] = 0.0;
        J[2][i] = 0.0;
        for (j=0; j<8; ++j) {
          J[0][i] += shape_derivatives[i][j]*node_positions[j][0];
          J[1][i] += shape_derivatives[i][j]*node_positions[j][1];
          J[2][i] += shape_derivatives[i][j]*node_positions[j][2];
        }
      }
    
      double detJ  =  1.0/(J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-
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
      shape_functions[0] = 0.125*(1.0-local_coords[0])*(1.0-local_coords[1])*(1.0-local_coords[2]);
      shape_functions[1] = 0.125*(1.0+local_coords[0])*(1.0-local_coords[1])*(1.0-local_coords[2]);
      shape_functions[2] = 0.125*(1.0+local_coords[0])*(1.0+local_coords[1])*(1.0-local_coords[2]);
      shape_functions[3] = 0.125*(1.0-local_coords[0])*(1.0+local_coords[1])*(1.0-local_coords[2]);
      shape_functions[4] = 0.125*(1.0-local_coords[0])*(1.0-local_coords[1])*(1.0+local_coords[2]);
      shape_functions[5] = 0.125*(1.0+local_coords[0])*(1.0-local_coords[1])*(1.0+local_coords[2]);
      shape_functions[6] = 0.125*(1.0+local_coords[0])*(1.0+local_coords[1])*(1.0+local_coords[2]);
      shape_functions[7] = 0.125*(1.0-local_coords[0])*(1.0+local_coords[1])*(1.0+local_coords[2]);
      f[0] = 0.0; 
      f[1] = 0.0; 
      f[2] = 0.0;
      for( int i=0 ; i<nnodes ; ++i ){
        for (int j=0; j<3; ++j) {
          f[j] += shape_functions[i]*node_positions[i][j];
        }
      }

      u  = f[0]-global_coords[0];
      v  = f[1]-global_coords[1];
      w  = f[2]-global_coords[2];
      u1 = u0-(invJ[0][0]*u+invJ[0][1]*v+invJ[0][2]*w);
      v1 = v0-(invJ[1][0]*u+invJ[1][1]*v+invJ[1][2]*w);
      w1 = w0-(invJ[2][0]*u+invJ[2][1]*v+invJ[2][2]*w);
      du = abs(u1-u0);
      dv = abs(v1-v0);
      dw = abs(w1-w0);
      u0 = u1;
      v0 = v1;
      w0 = w1;
      if (du<tolerance && dv<tolerance && dw<tolerance) converged = true;
      ++iterations;
    }
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    if (!converged) {
      std::cerr << "ContactHexElemL8<ActiveScalar>::Compute_Local_Coordinates() did not converge"
                << std::endl;
      std::cerr << "                     after "<<max_iterations
                << " iterations:  du = "<<du
                <<";  dv = "<<dv
                <<";  dw = "<<dw<<std::endl;
    }
#endif
    POSTCONDITION(converged);
    // If it's close to any of the edges, snap to it
    if (abs(u0)<1.0+spatial_tolerance) {
      u0 = min(u0, 1.0);
      u0 = max(u0,-1.0);
    }
    if (abs(v0)<1.0+spatial_tolerance) {
      v0 = min(v0, 1.0);
      v0 = max(v0,-1.0);
    }
    if (abs(w0)<1.0+spatial_tolerance) {
      w0 = min(w0, 1.0);
      w0 = max(w0,-1.0);
    }

    local_coords[0] = u0;
    local_coords[1] = v0;
    local_coords[2] = w0;
  }

  //
  // use newton's method to iterate starting from solution of previous solve so only 1 iteration will be done
  //
  int  iterations=0;
  int  max_iterations=400;
  bool converged = false;
  ActiveScalar tolerance = 1.0e-12;
  ActiveScalar u, u0=local_coords[0], u1, du;
  ActiveScalar v, v0=local_coords[1], v1, dv;
  ActiveScalar w, w0=local_coords[2], w1, dw;
  ActiveScalar f[3], J[3][3], invJ[3][3];
  ActiveScalar shape_derivatives[3][8];
  while (!converged && iterations<max_iterations) {
    active_local_coords[0] = u0;
    active_local_coords[1] = v0;
    active_local_coords[2] = w0;

    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives( active_local_coords, shape_derivatives );
    for (i=0; i<3; ++i) {
      J[0][i] = 0.0;
      J[1][i] = 0.0;
      J[2][i] = 0.0;
      for (j=0; j<8; ++j) {
        J[0][i] += shape_derivatives[i][j]*active_node_positions[j][0];
        J[1][i] += shape_derivatives[i][j]*active_node_positions[j][1];
        J[2][i] += shape_derivatives[i][j]*active_node_positions[j][2];
      }
    }
    
    ActiveScalar detJ  =  1.0/(J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-
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
    Compute_Global_Coords( active_node_positions, active_local_coords, f );
    u  = f[0]-active_global_coords[0];
    v  = f[1]-active_global_coords[1];
    w  = f[2]-active_global_coords[2];
    u1 = u0-(invJ[0][0]*u+invJ[0][1]*v+invJ[0][2]*w);
    v1 = v0-(invJ[1][0]*u+invJ[1][1]*v+invJ[1][2]*w);
    w1 = w0-(invJ[2][0]*u+invJ[2][1]*v+invJ[2][2]*w);
    du = abs(u1-u0);
    dv = abs(v1-v0);
    dw = abs(w1-w0);
    u0 = u1;
    v0 = v1;
    w0 = w1;
    if (du<tolerance && dv<tolerance && dw<tolerance) converged = true;
    ++iterations;
  }
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  if (!converged) {
    std::cerr << "ContactHexElemL8<ActiveScalar>::Compute_Local_Coordinates() did not converge"
              << std::endl;
    std::cerr << "                     after "<<max_iterations
              << " iterations:  du = "<<du
              <<";  dv = "<<dv
              <<";  dw = "<<dw<<std::endl;
  }
#endif
  POSTCONDITION(converged);
  active_local_coords[0] = u0;
  active_local_coords[1] = v0;
  active_local_coords[2] = w0;
}

template<typename DataType>
void ContactHexElemL8<DataType>::Compute_Global_Coords( DataType node_positions[8][3],
					      DataType local_coords[3],
					      DataType global_coords[3] )
{
  DataType N[8];
  int  nnodes=8;
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

template<typename DataType>
void  ContactHexElemL8<DataType>::Interpolate_Scalar( DataType  local_coords[3],
					    DataType  node_scalars[8],
					    DataType& interpolated_scalar )
{
  DataType N[8];
  int  nnodes=8;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

template<typename DataType>
void  ContactHexElemL8<DataType>::Interpolate_Vector( DataType local_coords[3],
					    DataType node_vectors[8][3],
					    DataType interpolated_vector[3] )
{
  DataType N[8];
  int  nnodes=8;
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


/*************************************************************************
 
         C A R T E S I A N    H E X   E L E M E N T

*************************************************************************/

ContactCartesianHexElementL8::ContactCartesianHexElementL8( 
                                     ContactFixedSizeAllocator* alloc,
                                     int Block_Index, 
				     int Index_in_Block, int key ) 
  : ContactElement( alloc, ContactSearch::CARTESIANHEXELEMENTL8,
		    Block_Index, Index_in_Block, key) 
{}

ContactCartesianHexElementL8* 
ContactCartesianHexElementL8::new_ContactCartesianHexElementL8(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactCartesianHexElementL8].New_Frag())
             ContactCartesianHexElementL8(alloc, Block_Index, 
                                          Index_in_Block, key);
}

void ContactCartesianHexElementL8_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactCartesianHexElementL8),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactCartesianHexElementL8 allocator" );
}

ContactCartesianHexElementL8::~ContactCartesianHexElementL8() {}

void
ContactCartesianHexElementL8::TetDice( int &ntets, Real thex[][4][3], 
                                       VariableHandle POSITION ) 
{
  switch(ntets) {
  case 5:
    {
      // works fine for planar faces
      ntets = 5;
      int node_index[5][4] = {{0,1,2,5},
			      {2,3,0,7},
			      {7,6,5,2},
			      {5,4,7,0},
			      {0,5,2,7}};
      
      for (int i=0; i<5; ++i) {
	for (int j=0; j<4; ++j) {
	  Real* position = nodes[node_index[i][j]]->Variable(POSITION);
	  for (int k=0; k<3; ++k) {
	    thex[i][j][k] = position[k];
	  }
	}
      }        
      break;
    }
  case 24:
    {
      // Use 24 tets to dice hex
      Real*     node_position;
      Real      centroid3[3];
      Real      centroid4[3];
      int       tetcount = 0;
      int node_index[6][4] = {{0,1,5,4},
			      {1,2,6,5},
			      {2,3,7,6},
			      {0,4,7,3},
			      {0,3,2,1},
			      {4,5,6,7}};
      
      // construct hex-centroid -- this is used for every face decomposition
      centroid4[0] = 0.;
      centroid4[1] = 0.;
      centroid4[2] = 0.;
      for( int i=0; i<8; ++i ) { // 8 nodes per hex
	node_position = nodes[i]->Variable(POSITION);
	centroid4[0] += node_position[0];
	centroid4[1] += node_position[1];
	centroid4[2] += node_position[2];
      }
      centroid4[0] = 0.125*centroid4[0];
      centroid4[1] = 0.125*centroid4[1];
      centroid4[2] = 0.125*centroid4[2];
      
      // now compute decomposition
      tetcount = 0;
      for( int j=0; j<6; ++j ) {// 6 faces on a hex
	centroid3[0] = 0.;
	centroid3[1] = 0.;
	centroid3[2] = 0.;
	for( int i=0; i<4; ++i ) { // 4 nodes per quad face
	  node_position = nodes[node_index[j][i]]->Variable(POSITION);
	  centroid3[0] += node_position[0];
	  centroid3[1] += node_position[1];
	  centroid3[2] += node_position[2];
	}
	centroid3[0] = 0.25*centroid3[0];
	centroid3[1] = 0.25*centroid3[1];
	centroid3[2] = 0.25*centroid3[2];
	
	for(int i = 0; i < 4; i ++){
	  int iplus = (i+1)%4;
	  node_position = nodes[node_index[j][i]]->Variable(POSITION);
	  thex[tetcount][0][0] = node_position[0];
	  thex[tetcount][0][1] = node_position[1];
	  thex[tetcount][0][2] = node_position[2];
	  
	  node_position = nodes[node_index[j][iplus]]->Variable(POSITION);
	  thex[tetcount][1][0] = node_position[0];
	  thex[tetcount][1][1] = node_position[1];
	  thex[tetcount][1][2] = node_position[2];
	  
	  thex[tetcount][2][0] = centroid4[0];
	  thex[tetcount][2][1] = centroid4[1];
	  thex[tetcount][2][2] = centroid4[2];
	  
	  thex[tetcount][3][0] = centroid3[0];
	  thex[tetcount][3][1] = centroid3[1];
	  thex[tetcount][3][2] = centroid3[2];
	  
	  ++tetcount;
	}
      } // end of loop on hex faces
      break;
    }
  default:
    {
      //should never get here and don't want to test for every interaction
      //so just do hard abort
      std::exit(1);
    }
  }
}

void
ContactCartesianHexElementL8::Compute_Volume( VariableHandle NODE_POSITION,
					      VariableHandle ELEMENT_VOLUME )
{
  using std::abs;
  Real* volume = Variable(ELEMENT_VOLUME);
  Real* n0_pos = Node(0)->Variable(NODE_POSITION);
  Real* n6_pos = Node(6)->Variable(NODE_POSITION);
  *volume = abs( (n0_pos[0]-n6_pos[0]) * (n0_pos[1]-n6_pos[1]) * 
		 (n0_pos[2]-n6_pos[2]) );
}

void
ContactCartesianHexElementL8::Compute_Centroid( VariableHandle NODE_POSITION,
					        VariableHandle CENTROID )
{
  Real* centroid = Variable(CENTROID);
  centroid[0] = 0.0;
  centroid[1] = 0.0;
  centroid[2] = 0.0;
  for (int i=0; i<8; ++i) {
    Real* position = Node(i)->Variable(NODE_POSITION);
    centroid[0] += position[0];
    centroid[1] += position[1];
    centroid[2] += position[2];
  }
  centroid[0] /= 8.0;
  centroid[1] /= 8.0;
  centroid[2] /= 8.0;
}


bool
ContactCartesianHexElementL8::Is_Local_Coordinates_Inside_Element( Real* local_coords )
{
  if( local_coords[0] >= -1.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >= -1.0 && local_coords[1] <= 1.0 &&
      local_coords[2] >= -1.0 && local_coords[2] <= 1.0 )
    return true;
  return false;
}

bool
ContactCartesianHexElementL8::Is_Local_Coordinates_Near_Element( Real* local_coords, Real tolerance )
{
  Real low_coord  = -(1.+tolerance);
  Real high_coord = 1.+tolerance;
  if( local_coords[0] >= low_coord && local_coords[0] <= high_coord &&
      local_coords[1] >= low_coord && local_coords[1] <= high_coord &&
      local_coords[2] >= low_coord && local_coords[2] <= high_coord )
    return true;
  return false;
}

void ContactCartesianHexElementL8::Evaluate_Shape_Functions( Real* local_coords,
						    Real* shape_functions )
{
  Compute_Shape_Functions(local_coords, shape_functions);
}

void ContactCartesianHexElementL8::Compute_Global_Coordinates( VariableHandle POSITION,
							       Real* local_coords,
							       Real* global_coords )
{
  Real node_positions[8][3];
  for(int i=0; i<Nodes_Per_Element(); ++i ){
    Real* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void ContactCartesianHexElementL8::Compute_Local_Coordinates( VariableHandle POSITION,
						     Real* global_coords,
						     Real* local_coords )
{
  int i, j;
  Real node_positions[8][3];
  for (i=0; i<Nodes_Per_Element(); ++i) {
    Real* node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}

void ContactCartesianHexElementL8::Interpolate_Scalar_Value( Real* local_coords,
							     Real* node_scalars,
							     Real& interpolated_scalar )
{
  Interpolate_Scalar( local_coords, node_scalars, interpolated_scalar );
  return;
}

/*************************************************************************/
/*************************************************************************/
/*                                                                       */
/* The following functions are supplied for doing generic computations   */
/* on a linear H8 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary H8 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

void ContactCartesianHexElementL8::Compute_Shape_Functions( Real* local_coords,
							    Real* shape_functions )
{
  shape_functions[0] = 0.125*(1.0-local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[1] = 0.125*(1.0+local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[2] = 0.125*(1.0+local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[3] = 0.125*(1.0-local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[4] = 0.125*(1.0-local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[5] = 0.125*(1.0+local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[6] = 0.125*(1.0+local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[7] = 0.125*(1.0-local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0+local_coords[2]);
}

void 
ContactCartesianHexElementL8::Compute_Local_Coords( Real node_positions[8][3], 
						    Real global_coords[3],
						    Real local_coords[3] )
{
  //
  // 1st check for coincidence with one of the face nodes
  //
  int  i;
  int  nnodes=8;


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
    local_coords[0] = -1.0;
    local_coords[1] = -1.0;
    local_coords[2] = -1.0;
    break;
  case 1:
    local_coords[0] =  1.0;
    local_coords[1] = -1.0;
    local_coords[2] = -1.0;
    break;
  case 2:
    local_coords[0] =  1.0;
    local_coords[1] =  1.0;
    local_coords[2] = -1.0;
    break;
  case 3:
    local_coords[0] = -1.0;
    local_coords[1] =  1.0;
    local_coords[2] = -1.0;
    break;
  case 4:
    local_coords[0] = -1.0;
    local_coords[1] = -1.0;
    local_coords[2] =  1.0;
    break;
  case 5:
    local_coords[0] =  1.0;
    local_coords[1] = -1.0;
    local_coords[2] =  1.0;
    break;
  case 6:
    local_coords[0] =  1.0;
    local_coords[1] =  1.0;
    local_coords[2] =  1.0;
    break;
  case 7:
    local_coords[0] = -1.0;
    local_coords[1] =  1.0;
    local_coords[2] =  1.0;
    break;
  }
  if (i<nnodes) return;

  //
  local_coords[0] = 2.*(global_coords[0] - node_positions[0][0])/
    (node_positions[1][0] - node_positions[0][0]) - 1.;
  local_coords[1] = 2.*(global_coords[1] - node_positions[0][1])/
    (node_positions[3][1] - node_positions[0][1]) - 1.;
  local_coords[2] = 2.*(global_coords[2] - node_positions[0][2])/
    (node_positions[4][2] - node_positions[0][2]) - 1.;
}

void ContactCartesianHexElementL8::Compute_Global_Coords( Real node_positions[8][3],
							  Real local_coords[3],
							  Real global_coords[3] )
{
  Real N[8];
  int  nnodes=8;
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

void  ContactCartesianHexElementL8::Interpolate_Scalar( Real  local_coords[3],
							Real  node_scalars[8],
							Real& interpolated_scalar )
{
  Real N[8];
  int  nnodes=8;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}


/*************************************************************************
 
         T R I L I N E A R    H E X   E L E M E N T

*************************************************************************/

ContactHexElementL8::ContactHexElementL8( ContactFixedSizeAllocator* alloc, 
                                          int Block_Index, 
				          int Index_in_Block, int key ) 
  : ContactElement( alloc, ContactSearch::HEXELEMENTL8, 
                    Block_Index, Index_in_Block,
		    key) 
{}

ContactHexElementL8* 
ContactHexElementL8::new_ContactHexElementL8(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key)
{
  return new (alloc[ContactSearch::ALLOC_ContactHexElementL8].New_Frag())
             ContactHexElementL8(alloc, Block_Index, Index_in_Block, key);
}

void ContactHexElementL8_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactHexElementL8),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactHexElementL8 allocator" );
}

ContactHexElementL8::~ContactHexElementL8() {}

void
ContactHexElementL8::TetDice( int &ntets, Real thex[5][4][3], 
                              VariableHandle POSITION) 
{
  switch(ntets) {
  case 5:
    {
      // works fine for planar faces
      int node_index[5][4] = {{0,1,2,5},
			      {2,3,0,7},
			      {7,6,5,2},
			      {5,4,7,0},
			      {0,5,2,7}};
      
      for (int i=0; i<5; ++i) {
	for (int j=0; j<4; ++j) {
	  Real* position = nodes[node_index[i][j]]->Variable(POSITION);
	  for (int k=0; k<3; ++k) {
	    thex[i][j][k] = position[k];
	  }
	}
      }        
      break;
    }
  case 24:
    {
      // Use 24 tets to dice hex
      Real*     node_position;
      Real      centroid3[3];
      Real      centroid4[3];
      int       tetcount = 0;
      int node_index[6][4] = {{0,1,5,4},
			      {1,2,6,5},
			      {2,3,7,6},
			      {0,4,7,3},
			      {0,3,2,1},
			      {4,5,6,7}};
      
      // construct hex-centroid -- this is used for every face decomposition
      centroid4[0] = 0.;
      centroid4[1] = 0.;
      centroid4[2] = 0.;
      for( int i=0; i<8; ++i ) { // 8 nodes per hex
	node_position = nodes[i]->Variable(POSITION);
	centroid4[0] += node_position[0];
	centroid4[1] += node_position[1];
	centroid4[2] += node_position[2];
      }
      centroid4[0] = 0.125*centroid4[0];
      centroid4[1] = 0.125*centroid4[1];
      centroid4[2] = 0.125*centroid4[2];
      
      // now compute decomposition
      tetcount = 0;
      for( int j=0; j<6; ++j ) {// 6 faces on a hex
	centroid3[0] = 0.;
	centroid3[1] = 0.;
	centroid3[2] = 0.;
	for( int i=0; i<4; ++i ) { // 4 nodes per quad face
	  node_position = nodes[node_index[j][i]]->Variable(POSITION);
	  centroid3[0] += node_position[0];
	  centroid3[1] += node_position[1];
	  centroid3[2] += node_position[2];
	}
	centroid3[0] = 0.25*centroid3[0];
	centroid3[1] = 0.25*centroid3[1];
	centroid3[2] = 0.25*centroid3[2];
	
	for(int i = 0; i < 4; i ++){
	  int iplus = (i+1)%4;
	  node_position = nodes[node_index[j][i]]->Variable(POSITION);
	  thex[tetcount][0][0] = node_position[0];
	  thex[tetcount][0][1] = node_position[1];
	  thex[tetcount][0][2] = node_position[2];
	  
	  node_position = nodes[node_index[j][iplus]]->Variable(POSITION);
	  thex[tetcount][1][0] = node_position[0];
	  thex[tetcount][1][1] = node_position[1];
	  thex[tetcount][1][2] = node_position[2];
	  
	  thex[tetcount][2][0] = centroid4[0];
	  thex[tetcount][2][1] = centroid4[1];
	  thex[tetcount][2][2] = centroid4[2];
	  
	  thex[tetcount][3][0] = centroid3[0];
	  thex[tetcount][3][1] = centroid3[1];
	  thex[tetcount][3][2] = centroid3[2];
	  
	  ++tetcount;
	}
      } // end of loop on hex faces
      break;
    }
  default:
    {
      //should never get here and don't want to test for every interaction
      //so just do hard abort
      std::exit(1);
    }

  }
}

void ContactHexElementL8::Compute_Volume( VariableHandle NODE_POSITION,
					  VariableHandle ELEMENT_VOLUME )
{
  Real* volume = Variable(ELEMENT_VOLUME);

  Real* n0 = Node(0)->Variable(NODE_POSITION);
  Real* n1 = Node(1)->Variable(NODE_POSITION);
  Real* n2 = Node(2)->Variable(NODE_POSITION);
  Real* n3 = Node(3)->Variable(NODE_POSITION);
  Real* n4 = Node(4)->Variable(NODE_POSITION);
  Real* n5 = Node(5)->Variable(NODE_POSITION);
  Real* n6 = Node(6)->Variable(NODE_POSITION);
  Real* n7 = Node(7)->Variable(NODE_POSITION);
  
  Real x1 = n0[0];
  Real x2 = n1[0];
  Real x3 = n2[0];
  Real x4 = n3[0];
  Real x5 = n4[0];
  Real x6 = n5[0];
  Real x7 = n6[0];
  Real x8 = n7[0];

  Real y1 = n0[1];
  Real y2 = n1[1];
  Real y3 = n2[1];
  Real y4 = n3[1];
  Real y5 = n4[1];
  Real y6 = n5[1];
  Real y7 = n6[1];
  Real y8 = n7[1];

  Real z1 = n0[2];
  Real z2 = n1[2];
  Real z3 = n2[2];
  Real z4 = n3[2];
  Real z5 = n4[2];
  Real z6 = n5[2];
  Real z7 = n6[2];
  Real z8 = n7[2];

  Real rx0 = (y2*((z6-z3)-(z4-z5))+y3*(z2-z4)+y4*((z3-z8)-
             (z5-z2))+y5*((z8-z6)-(z2-z4))+y6*(z5-z2)+
             y8*(z4-z5));
  Real rx1 = (y3*((z7-z4)-(z1-z6))+y4*(z3-z1)+y1*((z4-z5)-
             (z6-z3))+y6*((z5-z7)-(z3-z1))+y7*(z6-z3)+
             y5*(z1-z6));
  Real rx2 = (y4*((z8-z1)-(z2-z7))+y1*(z4-z2)+y2*((z1-z6)-
             (z7-z4))+y7*((z6-z8)-(z4-z2))+y8*(z7-z4)+
             y6*(z2-z7));
  Real rx3 = (y1*((z5-z2)-(z3-z8))+y2*(z1-z3)+y3*((z2-z7)-
             (z8-z1))+y8*((z7-z5)-(z1-z3))+y5*(z8-z1)+
             y7*(z3-z8));
  Real rx4 = (y8*((z4-z7)-(z6-z1))+y7*(z8-z6)+y6*((z7-z2)-
             (z1-z8))+y1*((z2-z4)-(z8-z6))+y4*(z1-z8)+
             y2*(z6-z1));
  Real rx5 = (y5*((z1-z8)-(z7-z2))+y8*(z5-z7)+y7*((z8-z3)-
             (z2-z5))+y2*((z3-z1)-(z5-z7))+y1*(z2-z5)+
             y3*(z7-z2));
  Real rx6 = (y6*((z2-z5)-(z8-z3))+y5*(z6-z8)+y8*((z5-z4)-
             (z3-z6))+y3*((z4-z2)-(z6-z8))+y2*(z3-z6)+
             y4*(z8-z3));
  Real rx7 = (y7*((z3-z6)-(z5-z4))+y6*(z7-z5)+y5*((z6-z1)-
             (z4-z7))+y4*((z1-z3)-(z7-z5))+y3*(z4-z7)+
             y1*(z5-z4));

  *volume = (x1 * rx0 +
	     x2 * rx1 +
	     x3 * rx2 +
	     x4 * rx3 +
	     x5 * rx4 +
	     x6 * rx5 +
	     x7 * rx6 +
	     x8 * rx7) / 12.0;
}

void
ContactHexElementL8::Compute_Centroid( VariableHandle NODE_POSITION,
				       VariableHandle CENTROID )
{
  Real* centroid = Variable(CENTROID);
  centroid[0] = 0.0;
  centroid[1] = 0.0;
  centroid[2] = 0.0;
  for (int i=0; i<8; ++i) {
    Real* position = Node(i)->Variable(NODE_POSITION);
    centroid[0] += position[0];
    centroid[1] += position[1];
    centroid[2] += position[2];
  }
  centroid[0] /= 8.0;
  centroid[1] /= 8.0;
  centroid[2] /= 8.0;
}

bool
ContactHexElementL8::Is_Local_Coordinates_Inside_Element( Real* local_coords )
{
  if( local_coords[0] >= -1.0 && local_coords[0] <= 1.0 &&
      local_coords[1] >= -1.0 && local_coords[1] <= 1.0 &&
      local_coords[2] >= -1.0 && local_coords[2] <= 1.0 )
    return true;
  return false;
}

bool
ContactHexElementL8::Is_Local_Coordinates_Near_Element( Real* local_coords, Real tolerance )
{
  Real low_coord  = -(1.+tolerance);
  Real high_coord = 1.+tolerance;
  if( local_coords[0] >= low_coord && local_coords[0] <= high_coord &&
      local_coords[1] >= low_coord && local_coords[1] <= high_coord &&
      local_coords[2] >= low_coord && local_coords[2] <= high_coord )
    return true;
  return false;
}

void ContactHexElementL8::Evaluate_Shape_Functions( Real* local_coords,
						    Real* shape_functions )
{
  Compute_Shape_Functions(local_coords, shape_functions);
}

void ContactHexElementL8::Compute_Global_Coordinates( VariableHandle POSITION,
						      Real* local_coords,
						      Real* global_coords )
{
  Real node_positions[8][3];
  for(int i=0; i<Nodes_Per_Element(); ++i ){
    Real* node_position = Node(i)->Variable(POSITION);
    for (int j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void ContactHexElementL8::Compute_Local_Coordinates( VariableHandle POSITION,
						     Real* global_coords,
						     Real* local_coords )
{
  int i, j;
  Real node_positions[8][3];
  for (i=0; i<Nodes_Per_Element(); ++i) {
    Real* node_position = Node(i)->Variable(POSITION);
    for (j=0; j<3; ++j) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Local_Coords(node_positions, global_coords, local_coords);
}


void ContactHexElementL8::Interpolate_Scalar_Value( Real* local_coords,
						    Real* node_scalars,
						    Real& interpolated_scalar )
{
  Interpolate_Scalar( local_coords, node_scalars, interpolated_scalar );
  return;
}

/*************************************************************************/
/*************************************************************************/
/*                                                                       */
/* The following functions are supplied for doing generic computations   */
/* on a linear H8 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary H8 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

void ContactHexElementL8::Compute_Shape_Functions( Real* local_coords,
						   Real* shape_functions )
{
  shape_functions[0] = 0.125*(1.0-local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[1] = 0.125*(1.0+local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[2] = 0.125*(1.0+local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[3] = 0.125*(1.0-local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0-local_coords[2]);
  shape_functions[4] = 0.125*(1.0-local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[5] = 0.125*(1.0+local_coords[0])*
                             (1.0-local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[6] = 0.125*(1.0+local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0+local_coords[2]);
  shape_functions[7] = 0.125*(1.0-local_coords[0])*
                             (1.0+local_coords[1])*
                             (1.0+local_coords[2]);
}

void ContactHexElementL8::Compute_Shape_Derivatives( Real* local_coords,
						     Real shape_derivs[3][8] )
{
  shape_derivs[0][0] = -0.125*(1.0-local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][1] =  0.125*(1.0-local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][2] =  0.125*(1.0+local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][3] = -0.125*(1.0+local_coords[1])*(1.0-local_coords[2]);
  shape_derivs[0][4] = -0.125*(1.0-local_coords[1])*(1.0+local_coords[2]);
  shape_derivs[0][5] =  0.125*(1.0-local_coords[1])*(1.0+local_coords[2]);
  shape_derivs[0][6] =  0.125*(1.0+local_coords[1])*(1.0+local_coords[2]);
  shape_derivs[0][7] = -0.125*(1.0+local_coords[1])*(1.0+local_coords[2]);

  shape_derivs[1][0] = -0.125*(1.0-local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][1] = -0.125*(1.0+local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][2] =  0.125*(1.0+local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][3] =  0.125*(1.0-local_coords[0])*(1.0-local_coords[2]);
  shape_derivs[1][4] = -0.125*(1.0-local_coords[0])*(1.0+local_coords[2]);
  shape_derivs[1][5] = -0.125*(1.0+local_coords[0])*(1.0+local_coords[2]);
  shape_derivs[1][6] =  0.125*(1.0+local_coords[0])*(1.0+local_coords[2]);
  shape_derivs[1][7] =  0.125*(1.0-local_coords[0])*(1.0+local_coords[2]);

  shape_derivs[2][0] = -0.125*(1.0-local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][1] = -0.125*(1.0+local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][2] = -0.125*(1.0+local_coords[0])*(1.0+local_coords[1]);
  shape_derivs[2][3] = -0.125*(1.0-local_coords[0])*(1.0+local_coords[1]);
  shape_derivs[2][4] =  0.125*(1.0-local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][5] =  0.125*(1.0+local_coords[0])*(1.0-local_coords[1]);
  shape_derivs[2][6] =  0.125*(1.0+local_coords[0])*(1.0+local_coords[1]);
  shape_derivs[2][7] =  0.125*(1.0-local_coords[0])*(1.0+local_coords[1]);
}

void 
ContactHexElementL8::Compute_Local_Coords( Real node_positions[8][3], 
					   Real global_coords[3],
					   Real local_coords[3] )
{
  using std::abs;
  using std::min;
  using std::max;

  int  i, j;
  int  nnodes=8;
  Real spatial_tolerance = 1.0e-10;

  // are we on a node?
  for (i=0; i<nnodes; ++i) {
    Real dx = node_positions[i][0]-global_coords[0];
    Real dy = node_positions[i][1]-global_coords[1];
    Real dz = node_positions[i][2]-global_coords[2];
    Real d  = std::sqrt(dx*dx+dy*dy+dz*dz);
    if (d<spatial_tolerance) break;
  }
  switch (i) {
  case 0:
    local_coords[0] = -1.0;
    local_coords[1] = -1.0;
    local_coords[2] = -1.0;
    break;
  case 1:
    local_coords[0] =  1.0;
    local_coords[1] = -1.0;
    local_coords[2] = -1.0;
    break;
  case 2:
    local_coords[0] =  1.0;
    local_coords[1] =  1.0;
    local_coords[2] = -1.0;
    break;
  case 3:
    local_coords[0] = -1.0;
    local_coords[1] =  1.0;
    local_coords[2] = -1.0;
    break;
  case 4:
    local_coords[0] = -1.0;
    local_coords[1] = -1.0;
    local_coords[2] =  1.0;
    break;
  case 5:
    local_coords[0] =  1.0;
    local_coords[1] = -1.0;
    local_coords[2] =  1.0;
    break;
  case 6:
    local_coords[0] =  1.0;
    local_coords[1] =  1.0;
    local_coords[2] =  1.0;
    break;
  case 7:
    local_coords[0] = -1.0;
    local_coords[1] =  1.0;
    local_coords[2] =  1.0;
    break;
  }
  if (i<nnodes) return;
  //
  // else use newton's method to iterate
  //
  int  iterations=0;
  int  max_iterations=400;
  bool converged = false;
  Real tolerance = 1.0e-12;
  Real u, u0=0.0, u1, du;
  Real v, v0=0.0, v1, dv;
  Real w, w0=0.0, w1, dw;
  Real f[3], J[3][3], invJ[3][3];
  Real shape_derivatives[3][8];
  while (!converged && iterations<max_iterations) {
    local_coords[0] = u0;
    local_coords[1] = v0;
    local_coords[2] = w0;
    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives( local_coords, shape_derivatives );
    for (i=0; i<3; ++i) {
      J[0][i] = 0.0;
      J[1][i] = 0.0;
      J[2][i] = 0.0;
      for (j=0; j<8; ++j) {
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

    Compute_Global_Coords( node_positions, local_coords, f );
    u  = f[0]-global_coords[0];
    v  = f[1]-global_coords[1];
    w  = f[2]-global_coords[2];
    u1 = u0-(invJ[0][0]*u+invJ[0][1]*v+invJ[0][2]*w);
    v1 = v0-(invJ[1][0]*u+invJ[1][1]*v+invJ[1][2]*w);
    w1 = w0-(invJ[2][0]*u+invJ[2][1]*v+invJ[2][2]*w);
    du = abs(u1-u0);
    dv = abs(v1-v0);
    dw = abs(w1-w0);
    u0 = u1;
    v0 = v1;
    w0 = w1;
    if (du<tolerance && dv<tolerance && dw<tolerance) converged = true;
    ++iterations;
  }
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  if (!converged) {
    std::cerr << "ContactHexElementL8::Compute_Local_Coordinates() did not converge" 
	 << std::endl;
    std::cerr << "                     after "<<max_iterations
         << " iterations:  du = "<<du
         <<";  dv = "<<dv
         <<";  dw = "<<dw<<std::endl;
  }
#endif
  POSTCONDITION(converged);
  // If it's close to any of the edges, snap to it
  if (abs(u0)<1.0+spatial_tolerance) {
    u0 = min(u0, 1.0);
    u0 = max(u0,-1.0);
  }
  if (abs(v0)<1.0+spatial_tolerance) {
    v0 = min(v0, 1.0);
    v0 = max(v0,-1.0);
  }
  if (abs(w0)<1.0+spatial_tolerance) {
    w0 = min(w0, 1.0);
    w0 = max(w0,-1.0);
  }
  local_coords[0] = u0;
  local_coords[1] = v0;
  local_coords[2] = w0;
}

void ContactHexElementL8::Compute_Global_Coords( Real node_positions[8][3],
						 Real local_coords[3],
						 Real global_coords[3] )
{
  Real N[8];
  int  nnodes=8;
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

void  ContactHexElementL8::Interpolate_Scalar( Real  local_coords[3],
					       Real  node_scalars[8],
					       Real& interpolated_scalar )
{
  Real N[8];
  int  nnodes=8;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; ++i ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

#endif  // #define ContactHexElementL8_C_
