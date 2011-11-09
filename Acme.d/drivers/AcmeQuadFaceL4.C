// $Id: AcmeQuadFaceL4.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeQuadFaceL4.h"
#include "AcmeNode.h"
#include "contact_assert.h"
#include <iostream.h>
#include <stdlib.h>
#include <math.h>
#include <new.h>

AcmeQuadFaceL4::AcmeQuadFaceL4( int Block_Index, 
				int Index_in_Block, 
                                int Exo_ID ) 
  : AcmeFace( AcmeEntity::QUADFACEL4, 
              Block_Index, Index_in_Block,
              Exo_ID, (AcmeNode**) &nodes, &Node_Ids[0] )
{
  number_node_connections = 4;
}

AcmeQuadFaceL4::~AcmeQuadFaceL4() {}

void AcmeQuadFaceL4::Evaluate_Shape_Functions( Real* local_coords,
					       Real* shape_functions )
{
  Compute_Shape_Functions(local_coords, shape_functions);
}

void AcmeQuadFaceL4::Compute_Global_Coordinates( Real* local_coords,
						 Real* global_coords )
{
  Real node_positions[4][3];
  for(int i=0; i<4; i++ ){
    Real* node_position = Node(i)->Coordinates();
    for (int j=0; j<3; j++) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void AcmeQuadFaceL4::Compute_Local_Coordinates( Real* global_coords,
						Real* local_coords )
{
  int i, j;
  Real node_positions[4][3];
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
/* on a linear Q4 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary Q4 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

void AcmeQuadFaceL4::Compute_Shape_Functions( Real* local_coords,
						 Real* shape_functions )
{
  shape_functions[0] = 0.25*(1.0-local_coords[0])*(1.0-local_coords[1]);
  shape_functions[1] = 0.25*(1.0+local_coords[0])*(1.0-local_coords[1]);
  shape_functions[2] = 0.25*(1.0+local_coords[0])*(1.0+local_coords[1]);
  shape_functions[3] = 0.25*(1.0-local_coords[0])*(1.0+local_coords[1]);
}

void AcmeQuadFaceL4::Compute_Shape_Derivatives( Real* local_coords,
						   Real shape_derivs[2][4] )
{
  shape_derivs[0][0] = -0.25*(1.0-local_coords[1]);
  shape_derivs[0][1] =  0.25*(1.0-local_coords[1]);
  shape_derivs[0][2] =  0.25*(1.0+local_coords[1]);
  shape_derivs[0][3] = -0.25*(1.0+local_coords[1]);
  
  shape_derivs[1][0] = -0.25*(1.0-local_coords[0]);
  shape_derivs[1][1] = -0.25*(1.0+local_coords[0]);
  shape_derivs[1][2] =  0.25*(1.0+local_coords[0]);
  shape_derivs[1][3] =  0.25*(1.0-local_coords[0]);
}

void 
AcmeQuadFaceL4::Compute_Local_Coords( Real node_positions[4][3], 
					 Real global_coords[3],
					 Real local_coords[3] )
{
  int  i, j;
  int  nnodes=4;
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
    local_coords[0] = -1.0;
    local_coords[1] = -1.0;
    break;
  case 1:
    local_coords[0] =  1.0;
    local_coords[1] = -1.0;
    break;
  case 2:
    local_coords[0] =  1.0;
    local_coords[1] =  1.0;
    break;
  case 3:
    local_coords[0] = -1.0;
    local_coords[1] =  1.0;
    break;
  }
  if (i<nnodes) {
    local_coords[2] = 0.0;
    return;
  }
  //
  // else use newton's method to iterate
  //
  
  int  iterations=0;
  int  max_iterations=500;
  bool converged = false;
  Real tolerance = 1.0e-12;
  Real s, s0=0.0, s1, ds=0.0; 
  Real t, t0=0.0, t1, dt=0.0;
  Real coords[3];
  Real J[3][2], f[3], shape_derivatives[2][4];
  
  while (!converged && iterations<max_iterations) {
    coords[0] = s0;
    coords[1] = t0;
    Compute_Global_Coords( node_positions, coords, f );
    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives(coords, shape_derivatives);
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
  if (fabs(s0)<1.0+spatial_tolerance) {
    s0 = MIN(s0, 1.0);
    s0 = MAX(s0,-1.0);
  }
  if (fabs(t0)<1.0+spatial_tolerance) {
    t0 = MIN(t0, 1.0);
    t0 = MAX(t0,-1.0);
  }
  local_coords[0] = s0;
  local_coords[1] = t0;
  local_coords[2] = 0.0;
  
}

void AcmeQuadFaceL4::Compute_Global_Coords( Real node_positions[4][3],
					       Real local_coords[2],
					       Real global_coords[3] )
{
  Real N[4];
  int  nnodes=4;
  global_coords[0] = 0.0;
  global_coords[1] = 0.0;
  global_coords[2] = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; i++ ){
    for (int j=0; j<3; j++) {
      global_coords[j] += N[i]*node_positions[i][j];
    }
  }
}

void  AcmeQuadFaceL4::Interpolate_Scalar( Real  local_coords[2],
					     Real  node_scalars[4],
					     Real& interpolated_scalar )
{
  Real N[4];
  int  nnodes=4;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; i++ ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

void  AcmeQuadFaceL4::Interpolate_Vector( Real local_coords[2],
					     Real node_vectors[4][3],
					     Real interpolated_vector[3] )
{
  Real N[4];
  int  nnodes=4;
  interpolated_vector[0] = 0.0;
  interpolated_vector[1] = 0.0;
  interpolated_vector[2] = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; i++ ){
    for (int j=0; j<3; j++) {
      interpolated_vector[j] += N[i]*node_vectors[i][j];
    }
  }
}
