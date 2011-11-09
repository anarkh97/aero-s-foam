// $Id: AcmeTriFaceQ6.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeTriFaceQ6_h_
#define _AcmeTriFaceQ6_h_

#include "AcmeFace.h"

class AcmeNode;

/* This class represents the quadratic six node triangle face with the
   following fortran numbering convention.

                              2
                              /\
                             /  \
                            /    \
                        E2 5      4 E1
                          /        \
                         /          \
                        /            \
                       0------3-------1
                              E0
*/

class AcmeTriFaceQ6 : public AcmeFace {
 public:
  AcmeTriFaceQ6( int blk_indx=-1, 
                    int indx_in_block=-1, int key=-1 );
  ~AcmeTriFaceQ6( );
  int Nodes_Per_Face() { return 6; };
  void Compute_Local_Coordinates( Real*, Real* );
  void Compute_Global_Coordinates( Real*, Real* );
  void Evaluate_Shape_Functions( Real* local_coords, Real* shape_funcs );
  
  static void Compute_Shape_Functions( Real local_coords[3], 
                                       Real shape_funcs[6] );
  
  static void Compute_Shape_Derivatives( Real local_coords[3],
                                         Real shape_derivatives[2][6] );
                                   
  static void Compute_Local_Coords( Real node_positions[6][3],
				    Real global_coords[3],
				    Real local_coords[3] );
  
  static void Compute_Global_Coords( Real node_positions[6][3],
				     Real local_coords[3],
				     Real global_coords[3] );

  static void Interpolate_Scalar( Real  local_coords[3],
				  Real  node_scalars[6],
				  Real& interpolated_scalar );

  static void Interpolate_Vector( Real local_coords[3],
				  Real node_vectors[6][3],
				  Real interpolated_vector[3] );

 protected:
 private:
  AcmeNode* nodes[6];
  int Node_Ids[6];

};

#endif // ifndef _AcmeTriFaceQ6_h_
