// $Id: AcmeQuadFaceL4.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeQuadFaceL4_h_
#define _AcmeQuadFaceL4_h_

#include "AcmeFace.h"

class AcmeNode;

/* This class represents the bilinear four node quadrilateral face with the
   following fortran numbering convention.

                             E2
                       3-------------2
                       |             |
                       |             |
                    E3 |             | E1
                       |             |
                       0-------------1
                             E0
*/

class AcmeQuadFaceL4 : public AcmeFace {

 public:
  AcmeQuadFaceL4(int blk_indx=-1, int indx_in_block=-1, int key=-1);
  ~AcmeQuadFaceL4();
  int Nodes_Per_Face() { return 4; };

  void Compute_Local_Coordinates( Real*, Real* );
  void Compute_Global_Coordinates( Real*, Real* );
  void Evaluate_Shape_Functions( Real* Local_Coords, Real* Shape_Funcs );
                      
  static void Compute_Shape_Functions( Real local_coords[3], 
                                       Real Shape_Funcs[4] );
  
  static void Compute_Shape_Derivatives( Real local_coords[3],
                                         Real shape_derivatives[2][4] );
  
  static void Compute_Local_Coords( Real node_positions[4][3],
				    Real global_coords[3],
				    Real local_coords[3] );
  
  static void Compute_Global_Coords( Real node_positions[4][3],
				     Real local_coords[3],
				     Real global_coords[3] );

  static void Interpolate_Scalar( Real  local_coords[3],
				  Real  node_scalars[4],
				  Real& interpolated_scalar );

  static void Interpolate_Vector( Real local_coords[3],
				  Real node_vectors[4][3],
				  Real interpolated_vector[3] );

 protected:
 private:
  AcmeNode* nodes[4];
  int Node_Ids[4];

};


#endif // ifdef _AcmeQuadFaceL4_h_
