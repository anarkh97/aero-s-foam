// $Id: AcmeHexElemL8.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeHexElemL8.h"
#include "AcmeNode.h"
#include "contact_assert.h"
#include <iostream.h>
#include <stdlib.h>
#include <math.h>
#include <new.h>

AcmeHexElemL8::AcmeHexElemL8( int Block_Index, 
			      int Host_Index_in_Block, 
                              int Exo_ID ) 
  : AcmeElem( AcmeEntity::HEXELEML8, 
              Block_Index, Host_Index_in_Block,
              Exo_ID, (AcmeNode**) &nodes, &Node_Ids[0] ) 
{
  number_node_connections = 8;
}

AcmeHexElemL8::~AcmeHexElemL8() {}

void AcmeHexElemL8::Evaluate_Shape_Functions( Real* local_coords,
					      Real* shape_functions )
{
  Compute_Shape_Functions(local_coords, shape_functions);
}

void AcmeHexElemL8::Compute_Global_Coordinates( Real* local_coords,
						Real* global_coords )
{
  Real node_positions[8][3];
  for(int i=0; i<number_node_connections; i++ ){
    Real* node_position = Node(i)->Coordinates();
    for (int j=0; j<3; j++) {
      node_positions[i][j] = node_position[j];
    }
  }
  Compute_Global_Coords(node_positions, local_coords, global_coords);
}

void AcmeHexElemL8::Compute_Local_Coordinates( Real* global_coords,
					       Real* local_coords )
{
  int i, j;
  Real node_positions[8][3];
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
/* on a linear H8 that isn't created as an object.  This is useful for   */
/* for normal smoothing and other situations where you want to compute   */
/* on a temporary H8 with out the necessity of creating nodes and edges. */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/

void AcmeHexElemL8::Compute_Shape_Functions( Real* local_coords,
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

void AcmeHexElemL8::Compute_Shape_Derivatives( Real* local_coords,
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
AcmeHexElemL8::Compute_Local_Coords( Real node_positions[8][3], 
					Real global_coords[3],
					Real local_coords[3] )
{
  //
  // 1st check for coincidence with one of the face nodes
  //
  int  i, j;
  int  nnodes=8;
  Real spatial_tolerance = 1.0e-10;

  // are we on a node?
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
  Real f[3], J[3][3], invJ[3][3], detJ;
  Real shape_derivatives[3][8];
  while (!converged && iterations<max_iterations) {
    local_coords[0] = u0;
    local_coords[1] = v0;
    local_coords[2] = w0;
    // BUILD JACOBIAN AND INVERT
    Compute_Shape_Derivatives( local_coords, shape_derivatives );
    for (i=0; i<3; i++) {
      J[0][i] = 0.0;
      J[1][i] = 0.0;
      J[2][i] = 0.0;
      for (j=0; j<8; j++) {
        J[0][i] += shape_derivatives[i][j]*node_positions[j][0];
        J[1][i] += shape_derivatives[i][j]*node_positions[j][1];
        J[2][i] += shape_derivatives[i][j]*node_positions[j][2];
      }
    }
    
    detJ       =  J[0][0]*(J[1][1]*J[2][2]-J[1][2]*J[2][1])-
                  J[0][1]*(J[1][0]*J[2][2]-J[2][0]*J[1][2])+
                  J[0][2]*(J[1][0]*J[2][1]-J[2][0]*J[1][1]);
    
    POSTCONDITION(fabs(detJ) > TINYNUM);
    
    invJ[0][0] =  (J[1][1]*J[2][2]-J[1][2]*J[2][1])/detJ;
    invJ[0][1] = -(J[0][1]*J[2][2]-J[2][1]*J[0][2])/detJ;
    invJ[0][2] =  (J[1][2]*J[0][1]-J[0][2]*J[1][1])/detJ;
    invJ[1][0] = -(J[1][0]*J[2][2]-J[2][0]*J[1][2])/detJ;
    invJ[1][1] =  (J[0][0]*J[2][2]-J[0][2]*J[2][0])/detJ;
    invJ[1][2] = -(J[0][0]*J[1][2]-J[1][0]*J[0][2])/detJ;
    invJ[2][0] =  (J[1][0]*J[2][1]-J[2][0]*J[1][1])/detJ;
    invJ[2][1] = -(J[0][0]*J[2][1]-J[2][0]*J[0][1])/detJ;
    invJ[2][2] =  (J[0][0]*J[1][1]-J[0][1]*J[1][0])/detJ;

    // APPLY NEWTON ALGORITHM

    Compute_Global_Coords( node_positions, local_coords, f );
    u  = f[0]-global_coords[0];
    v  = f[1]-global_coords[1];
    w  = f[2]-global_coords[2];
    u1 = u0-(invJ[0][0]*u+invJ[0][1]*v+invJ[0][2]*w);
    v1 = v0-(invJ[1][0]*u+invJ[1][1]*v+invJ[1][2]*w);
    w1 = w0-(invJ[2][0]*u+invJ[2][1]*v+invJ[2][2]*w);
    du = fabs(u1-u0);
    dv = fabs(v1-v0);
    dw = fabs(w1-w0);
    u0 = u1;
    v0 = v1;
    w0 = w1;
    if (du<tolerance && dv<tolerance && dw<tolerance) converged = true;
    iterations++;
  }
  POSTCONDITION(converged);
  // If it's close to any of the edges, snap to it
  if (fabs(u0)<1.0+spatial_tolerance) {
    u0 = MIN(u0, 1.0);
    u0 = MAX(u0,-1.0);
  }
  if (fabs(v0)<1.0+spatial_tolerance) {
    v0 = MIN(v0, 1.0);
    v0 = MAX(v0,-1.0);
  }
  if (fabs(w0)<1.0+spatial_tolerance) {
    w0 = MIN(w0, 1.0);
    w0 = MAX(w0,-1.0);
  }
  local_coords[0] = u0;
  local_coords[1] = v0;
  local_coords[2] = w0;
}

void AcmeHexElemL8::Compute_Global_Coords( Real node_positions[8][3],
					      Real local_coords[3],
					      Real global_coords[3] )
{
  Real N[8];
  int  nnodes=8;
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

void  AcmeHexElemL8::Interpolate_Scalar( Real  local_coords[3],
					    Real  node_scalars[8],
					    Real& interpolated_scalar )
{
  Real N[8];
  int  nnodes=8;
  interpolated_scalar = 0.0;
  Compute_Shape_Functions(local_coords, N);
  for( int i=0 ; i<nnodes ; i++ ){
    interpolated_scalar += N[i]*node_scalars[i];
  }
}

void  AcmeHexElemL8::Interpolate_Vector( Real local_coords[3],
					    Real node_vectors[8][3],
					    Real interpolated_vector[3] )
{
  Real N[8];
  int  nnodes=8;
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
