// $Id: AcmeQuadFaceQ8.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeQuadFaceQ8_h_
#define _AcmeQuadFaceQ8_h_

#include "AcmeFace.h"

class AcmeNode;

/* This class represents the quadratic eight node quadrilateral face with the
   following fortran numbering convention.

                             E2
                       3------6------2
                       |             |
                       |             |
                       7             5
                    E3 |             | E1
                       |             |
                       0------4------1
                             E0
*/

class AcmeQuadFaceQ8 : public AcmeFace {

 public:
  AcmeQuadFaceQ8(int blk_indx=-1, 
                    int indx_in_block=-1, int key=-1);
  ~AcmeQuadFaceQ8();
  int Nodes_Per_Face() { return 8; };

  void Compute_Local_Coordinates( Real*, Real* );
  void Compute_Global_Coordinates( Real*, Real* );
  void Evaluate_Shape_Functions( Real* local_coords, Real* shape_funcs );
  
  static void Compute_Shape_Functions( Real local_coords[2], 
                                       Real shape_funcs[8] );
  
  static void Compute_Shape_Derivatives( Real local_coords[2],
                                         Real shape_derivatives[2][8] );
  
  static void Compute_Local_Coords( Real node_positions[8][3],
				    Real global_coords[3],
				    Real local_coords[2] );
  
  static void Compute_Global_Coords( Real node_positions[8][3],
				     Real local_coords[2],
				     Real global_coords[3] );

  static void Interpolate_Scalar( Real  local_coords[2],
				  Real  node_scalars[8],
				  Real& interpolated_scalar );

  static void Interpolate_Vector( Real local_coords[2],
				  Real node_vectors[8][3],
				  Real interpolated_vector[3] );

 protected:
 private:
  AcmeNode* nodes[8];
  int Node_Ids[8];

};

#endif // ifdef _AcmeQuadFaceQ8_h_

