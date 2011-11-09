// $Id: AcmeTriFaceL3.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeTriFaceL3.h"
#include "AcmeNode.h"
#include "contact_assert.h"
#include <iostream.h>
#include <math.h>
#include <new.h>

AcmeTriFaceL3::AcmeTriFaceL3( int Block_Index, 
		              int Index_in_Block,
                              int Exo_ID ) 
  : AcmeFace( AcmeEntity::TRIFACEL3, 
              Block_Index, Index_in_Block,
              Exo_ID, (AcmeNode**) &nodes, &Node_Ids[0] ) 
{
  number_node_connections = 3;
}

AcmeTriFaceL3::~AcmeTriFaceL3() {}

void AcmeTriFaceL3::Evaluate_Shape_Functions( Real* local_coords,
						 Real* shape_functions )
{
  
  PRECONDITION( fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  Compute_Shape_Functions( local_coords, shape_functions );
}

void AcmeTriFaceL3::Compute_Global_Coordinates( Real* local_coords,
						Real* global_coords )
{
  Real node_positions[3][3];
  for(int i=0; i<3; i++ ){
    Real* node_position = Node(i)->Coordinates();
    for (int j=0; j<3; j++) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void AcmeTriFaceL3::Compute_Local_Coordinates( Real* global_coords,
					       Real* local_coords )
{
  int  i, j;
  Real node_positions[3][3];
  for (i=0; i<number_node_connections; i++) {
    Real* node_position = Node(i)->Coordinates();
    for (j=0; j<3; j++) {
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

void AcmeTriFaceL3::Compute_Shape_Functions( Real local_coords[3],
						Real shape_functions[3] )
{
  
  PRECONDITION( fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  shape_functions[0] = local_coords[0];
  shape_functions[1] = local_coords[1];
  shape_functions[2] = local_coords[2];
}

void AcmeTriFaceL3::Compute_Shape_Derivatives( Real local_coords[3],
					          Real shape_derivs[2][3] )
{
  
  PRECONDITION( fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  shape_derivs[0][0] =  1.0;
  shape_derivs[0][1] =  0.0;
  shape_derivs[0][2] = -1.0;
  shape_derivs[1][0] =  0.0;
  shape_derivs[1][1] =  1.0;
  shape_derivs[1][2] = -1.0;
}

void AcmeTriFaceL3::Compute_Local_Coords( Real node_positions[3][3],
					     Real global_coords[3],
					     Real local_coords[3] )
{
  int  i;
  Real spatial_tolerance = 1.0e-10;
  //
  // check for coincidence with one of the face nodes
  //
  for (i=0; i<3; i++) {
    Real dx = node_positions[i][0]-global_coords[0];
    Real dy = node_positions[i][1]-global_coords[1];
    Real dz = node_positions[i][2]-global_coords[2];
    Real d  = sqrt(dx*dx+dy*dy+dz*dz);
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
  Real v_c_0_x = node_positions[0][0] - global_coords[0];
  Real v_c_0_y = node_positions[0][1] - global_coords[1];
  Real v_c_0_z = node_positions[0][2] - global_coords[2];
  Real v_c_1_x = node_positions[1][0] - global_coords[0];
  Real v_c_1_y = node_positions[1][1] - global_coords[1];
  Real v_c_1_z = node_positions[1][2] - global_coords[2];
  Real v_c_2_x = node_positions[2][0] - global_coords[0];
  Real v_c_2_y = node_positions[2][1] - global_coords[1];
  Real v_c_2_z = node_positions[2][2] - global_coords[2];

  // Compute the normal (as v_0_1 x v_0_2)
  Real Normal[3];
  Real v_0_1_x = node_positions[1][0] - node_positions[0][0];
  Real v_0_1_y = node_positions[1][1] - node_positions[0][1];
  Real v_0_1_z = node_positions[1][2] - node_positions[0][2];
  Real v_0_2_x = node_positions[2][0] - node_positions[0][0];
  Real v_0_2_y = node_positions[2][1] - node_positions[0][1];
  Real v_0_2_z = node_positions[2][2] - node_positions[0][2];
  Normal[0] = v_0_1_y*v_0_2_z - v_0_1_z*v_0_2_y;
  Normal[1] = v_0_1_z*v_0_2_x - v_0_1_x*v_0_2_z;
  Normal[2] = v_0_1_x*v_0_2_y - v_0_1_y*v_0_2_x;

  // Compute the areas
  Real a0 = ( v_c_1_y*v_c_2_z-v_c_2_y*v_c_1_z)*Normal[0] +
            (-v_c_1_x*v_c_2_z+v_c_2_x*v_c_1_z)*Normal[1] +
            ( v_c_1_x*v_c_2_y-v_c_2_x*v_c_1_y)*Normal[2] ;
  Real a1 = ( v_c_2_y*v_c_0_z-v_c_0_y*v_c_2_z)*Normal[0] +
            (-v_c_2_x*v_c_0_z+v_c_0_x*v_c_2_z)*Normal[1] +
            ( v_c_2_x*v_c_0_y-v_c_0_x*v_c_2_y)*Normal[2] ;
  Real a2 = ( v_c_0_y*v_c_1_z-v_c_1_y*v_c_0_z)*Normal[0] +
            (-v_c_0_x*v_c_1_z+v_c_1_x*v_c_0_z)*Normal[1] +
            ( v_c_0_x*v_c_1_y-v_c_1_x*v_c_0_y)*Normal[2] ;

  // Compute the local coordinates
  Real area = a0 + a1 + a2;
  Real L1   = a0/area;
  Real L2   = a1/area;
  // If it's close to any of the edges, snap to it
  if (L1<1.0+spatial_tolerance) {
    L1 = MIN(L1, 1.0);
  }
  if (L1>-spatial_tolerance) {
    L1 = MAX(L1, 0.0);
  }
  if (L2<1.0+spatial_tolerance) {
    L2 = MIN(L2, 1.0);
  }
  if (L2>-spatial_tolerance) {
    L2 = MAX(L2, 0.0);
  }
  local_coords[0] = L1;
  local_coords[1] = L2;
  local_coords[2] = 1.0-L1-L2;
}

void AcmeTriFaceL3::Compute_Global_Coords( Real node_positions[3][3],
					      Real local_coords[3],
					      Real global_coords[3] )
{
  Real N[3];
  int  nnodes=3;
  global_coords[0] = 0.0;
  global_coords[1] = 0.0;
  global_coords[2] = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; i++ ){
    for (int j=0; j<3; j++) {
      global_coords[j] += N[i]*node_positions[i][j];
    }
  }
}

void AcmeTriFaceL3::Interpolate_Scalar( Real  local_coords[3],
					   Real  node_scalars[3],
					   Real& interpolated_scalar )
{
  Real N[3];
  int  nnodes=3;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; i++ ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

void AcmeTriFaceL3::Interpolate_Vector( Real local_coords[3],
					   Real node_vectors[3][3],
					   Real interpolated_vector[3] )
{
  Real N[3];
  int  nnodes=3;
  interpolated_vector[0] = 0.0;
  interpolated_vector[1] = 0.0;
  interpolated_vector[2] = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; i++ ){
    for (int j=0; j<3; j++) {
      interpolated_vector[j] += N[i]*node_vectors[i][j];
    }
  }
}
