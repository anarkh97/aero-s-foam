// $Id: AcmeTriFaceQ6.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeTriFaceQ6.h"
#include "AcmeNode.h"
#include "contact_assert.h"
#include <iostream.h>
#include <stdlib.h>
#include <math.h>
#include <new.h>

AcmeTriFaceQ6::AcmeTriFaceQ6( int Block_Index, 
			      int Index_in_Block,
                              int Exo_ID ) 
  : AcmeFace( AcmeEntity::TRIFACEQ6, 
              Block_Index, Index_in_Block,
              Exo_ID, (AcmeNode**) &nodes, &Node_Ids[0] ) 
{
  number_node_connections = 6;
}

AcmeTriFaceQ6::~AcmeTriFaceQ6() {}

void AcmeTriFaceQ6::Evaluate_Shape_Functions( Real* local_coords,
						 Real* shape_functions )
{
  Compute_Shape_Functions( local_coords, shape_functions );
}

void AcmeTriFaceQ6::Compute_Global_Coordinates( Real* local_coords,
                                                Real* global_coords )
{
  Real node_positions[6][3];
  for(int i=0; i<6; i++ ){
    Real* node_position = Node(i)->Coordinates();
    for (int j=0; j<3; j++) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void AcmeTriFaceQ6::Compute_Local_Coordinates( Real* global_coords,
					       Real* local_coords )
{
  int  i, j;
  Real node_positions[6][3];
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

void AcmeTriFaceQ6::Compute_Shape_Functions( Real local_coords[3],
						Real shape_functions[6] )
{
  PRECONDITION( fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  Real a0 = local_coords[0];
  Real a1 = local_coords[1];
  Real a2 = local_coords[2];
  shape_functions[0] = a0*(2.0*a0-1.0);
  shape_functions[1] = a1*(2.0*a1-1.0);
  shape_functions[2] = a2*(2.0*a2-1.0);
  shape_functions[3] = 4.0*a0*a1;
  shape_functions[4] = 4.0*a1*a2;
  shape_functions[5] = 4.0*a2*a0;
}

void AcmeTriFaceQ6::Compute_Shape_Derivatives( Real local_coords[3],
                                                  Real shape_derivatives[2][6])
{
  PRECONDITION( fabs(local_coords[0]+local_coords[1]+local_coords[2]-1.0)<1.e-13 );
  Real a0 = local_coords[0];
  Real a1 = local_coords[1];
  Real a2 = local_coords[2];
  shape_derivatives[0][0] =  4.0*a0 - 1.0;
  shape_derivatives[0][1] =  0.0;
  shape_derivatives[0][2] = -4.0*a2 + 1.0;
  shape_derivatives[0][3] =  4.0*a1;
  shape_derivatives[0][4] = -4.0*a1;
  shape_derivatives[0][5] =  4.0*(a2-a0);

  shape_derivatives[1][0] =  0.0;
  shape_derivatives[1][1] =  4.0*a1 - 1.0;
  shape_derivatives[1][2] = -4.0*a2 + 1.0;
  shape_derivatives[1][3] =  4.0*a0;
  shape_derivatives[1][4] =  4.0*(a2-a1);
  shape_derivatives[1][5] = -4.0*a0;
}

void AcmeTriFaceQ6::Compute_Local_Coords( Real node_positions[6][3],
                                             Real global_coords[3],
                                             Real local_coords[3] )
{
  int  i, j;
  int  nnodes=6;
  Real spatial_tolerance = 1.0e-10;
  //
  // check for coincidence with one of the face nodes
  //
  for (i=0; i<nnodes; i++) {
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
  case 3:
    local_coords[0] = 0.5;
    local_coords[1] = 0.5;
    local_coords[2] = 0.0;
    break;
  case 4:
    local_coords[0] = 0.0;
    local_coords[1] = 0.5;
    local_coords[2] = 0.5;
    break;
  case 5:
    local_coords[0] = 0.5;
    local_coords[1] = 0.0;
    local_coords[2] = 0.5;
    break;
  }
  if (i<nnodes) return;
  /*
     else use newton's method to iterate
    
     Assume straight-sided triangle for first guess
    
                      2
                     /|\
                    / | \
                   /a1|a0\
                  /  /c\  \
                 / /     \ \
                //    a2   \\ 
               0-------------1
  */

  // Compute vectors from the contact point to each node
  Real v_c_0_x = node_positions[0][0] - global_coords[0];
  Real v_c_0_y = node_positions[0][1] - global_coords[1];
  Real v_c_0_z = node_positions[0][2] - global_coords[2];
  Real v_c_1_x = node_positions[1][0] - global_coords[0];
  Real v_c_1_y = node_positions[1][1] - global_coords[1];
  Real v_c_1_z = node_positions[1][2] - global_coords[2];
  Real v_c_2_x = node_positions[2][0] - global_coords[0];
  Real v_c_2_y = node_positions[2][1] - global_coords[1];
  Real v_c_2_z = node_positions[2][2] - global_coords[2];

  // Compute the areas
  Real shape_derivatives[2][6];
  Real midpnt[3] = {1.0/3.0, 1.0/3.0, 1.0 - 1.0/3.0 - 1.0/3.0};
  Compute_Shape_Derivatives( midpnt, shape_derivatives );
  Real e1[3] = {0.0, 0.0, 0.0};
  Real e2[3] = {0.0, 0.0, 0.0};
  for (i=0; i<nnodes; i++) {
    e1[0] += shape_derivatives[0][i]*node_positions[i][0];
    e1[1] += shape_derivatives[0][i]*node_positions[i][1];
    e1[2] += shape_derivatives[0][i]*node_positions[i][2];
    e2[0] += shape_derivatives[1][i]*node_positions[i][0];
    e2[1] += shape_derivatives[1][i]*node_positions[i][1];
    e2[2] += shape_derivatives[1][i]*node_positions[i][2];
  }
  Real a=0.0, b=0.0, c=0.0;
  for (i=0; i<3; i++) {
    a += e1[i]*e1[i];
    b += e2[i]*e2[i];
    c += e1[i]*e2[i];
  }
  Real face_normal[3];
  Real detJ      = sqrt(a*b-c*c);
  face_normal[0] = (e1[1]*e2[2]-e1[2]*e2[1])/detJ;
  face_normal[1] = (e2[0]*e1[2]-e1[0]*e2[2])/detJ;
  face_normal[2] = (e1[0]*e2[1]-e2[0]*e1[1])/detJ;
  Real a0        = ( v_c_1_y*v_c_2_z-v_c_2_y*v_c_1_z)*face_normal[0] +
                   (-v_c_1_x*v_c_2_z+v_c_2_x*v_c_1_z)*face_normal[1] +
                   ( v_c_1_x*v_c_2_y-v_c_2_x*v_c_1_y)*face_normal[2] ;
  Real a1        = ( v_c_2_y*v_c_0_z-v_c_0_y*v_c_2_z)*face_normal[0] +
                   (-v_c_2_x*v_c_0_z+v_c_0_x*v_c_2_z)*face_normal[1] +
                   ( v_c_2_x*v_c_0_y-v_c_0_x*v_c_2_y)*face_normal[2] ;
  Real a2        = ( v_c_0_y*v_c_1_z-v_c_1_y*v_c_0_z)*face_normal[0] +
                   (-v_c_0_x*v_c_1_z+v_c_1_x*v_c_0_z)*face_normal[1] +
                   ( v_c_0_x*v_c_1_y-v_c_1_x*v_c_0_y)*face_normal[2] ;

  // Compute the local coordinates
  Real area = a0 + a1 + a2;
  a0 /= area;
  a1 /= area;
  a2 /= area;

  int  iterations=0;
  int  max_iterations=200;
  bool converged = false;
  Real tolerance = 1.0e-12;
  Real s, s0=a0, s1, ds=0.0;
  Real t, t0=a1, t1, dt=0.0;
  Real coords[3];
  Real J[3][2], f[3];
  while (!converged && iterations<max_iterations) {
    coords[0] = s0;
    coords[1] = t0;
    coords[2] = 1.0-s0-t0;
    Compute_Global_Coords(node_positions, coords, f );
    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives( coords, shape_derivatives );
    for (i=0; i<2; i++) {
      J[0][i] = 0.0;
      J[1][i] = 0.0;
      J[2][i] = 0.0;
      for (j=0; j<nnodes; j++) {
        J[0][i] += shape_derivatives[i][j]*node_positions[j][0];
        J[1][i] += shape_derivatives[i][j]*node_positions[j][1];
        J[2][i] += shape_derivatives[i][j]*node_positions[j][2];
      }
    }
    Real JT[2][3];
    JT[0][0] = J[0][0];
    JT[0][1] = J[1][0];
    JT[0][2] = J[2][0];
    JT[1][0] = J[0][1];
    JT[1][1] = J[1][1];
    JT[1][2] = J[2][1];
    
    Real JTJ[2][2];
    JTJ[0][0] = JT[0][0]*J[0][0] + JT[0][1]*J[1][0] + JT[0][2]*J[2][0];
    JTJ[0][1] = JT[0][0]*J[0][1] + JT[0][1]*J[1][1] + JT[0][2]*J[2][1];
    JTJ[1][0] = JT[1][0]*J[0][0] + JT[1][1]*J[1][0] + JT[1][2]*J[2][0];
    JTJ[1][1] = JT[1][0]*J[0][1] + JT[1][1]*J[1][1] + JT[1][2]*J[2][1];
    
    Real invJTJ[2][2];
    Real detJTJ  = JTJ[0][0]*JTJ[1][1]-JTJ[0][1]*JTJ[1][0];
    invJTJ[0][0] =  JTJ[1][1]/detJTJ;
    invJTJ[0][1] = -JTJ[0][1]/detJTJ;
    invJTJ[1][0] = -JTJ[1][0]/detJTJ;
    invJTJ[1][1] =  JTJ[0][0]/detJTJ;

    // APPLY NEWTON ALGORITHM

    Real dx = f[0]-global_coords[0];
    Real dy = f[1]-global_coords[1];
    Real dz = f[2]-global_coords[2];
    s = JT[0][0]*dx + JT[0][1]*dy + JT[0][2]*dz;
    t = JT[1][0]*dx + JT[1][1]*dy + JT[1][2]*dz;
    
    s1 = s0-(invJTJ[0][0]*s+invJTJ[0][1]*t);
    t1 = t0-(invJTJ[1][0]*s+invJTJ[1][1]*t);
    ds = fabs(s1-s0);
    dt = fabs(t1-t0);
    s0 = s1;
    t0 = t1;
    if (ds<tolerance && dt<tolerance) converged = true;
    iterations++;
  }
  // If it's close to any of the edges, snap to it
  if (s0<1.0+spatial_tolerance) {
    s0 = MIN(s0, 1.0);
  }
  if (s0>-spatial_tolerance) {
    s0 = MAX(s0, 0.0);
  }
  if (t0<1.0+spatial_tolerance) {
    t0 = MIN(t0, 1.0);
  }
  if (t0>-spatial_tolerance) {
    t0 = MAX(t0, 0.0);
  }
  local_coords[0] = s0;
  local_coords[1] = t0;
  local_coords[2] = 1.0-s0-t0;
}

void AcmeTriFaceQ6::Compute_Global_Coords( Real node_positions[6][3],
					      Real local_coords[3],
					      Real global_coords[3] )
{
  Real N[6];
  int  nnodes=6;
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

void AcmeTriFaceQ6::Interpolate_Scalar( Real  local_coords[3],
					   Real  node_scalars[6],
					   Real& interpolated_scalar )
{
  Real N[6];
  int  nnodes=6;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions( local_coords, N );
  for( int i=0 ; i<nnodes ; i++ ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

void AcmeTriFaceQ6::Interpolate_Vector( Real local_coords[3],
					   Real node_vectors[6][3],
					   Real interpolated_vector[3] )
{
  Real N[6];
  int  nnodes=6;
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


