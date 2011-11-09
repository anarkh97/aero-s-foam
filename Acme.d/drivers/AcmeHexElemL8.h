// $Id: AcmeHexElemL8.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeHexElemL8_h_
#define _AcmeHexElemL8_h_

#include "AcmeElem.h"

class AcmeNode;


class AcmeHexElemL8 : public AcmeElem {

 public:
  AcmeHexElemL8(int blk_indx=-1, int indx_in_block=-1, int key=-1);
  ~AcmeHexElemL8();
  
  int Nodes_Per_Element() { return 8; };
  void Evaluate_Shape_Functions( Real*, Real* );
  void Compute_Local_Coordinates( Real*, Real* );
  void Compute_Global_Coordinates( Real*, Real* );
                      
  static void Compute_Shape_Functions( Real local_coords[3], 
                                       Real Shape_Funcs[8] );
  
  static void Compute_Shape_Derivatives( Real local_coords[3],
                                         Real shape_derivatives[3][8] );
  
  static void Compute_Local_Coords( Real node_positions[8][3],
				    Real global_coords[3],
				    Real local_coords[3] );
  
  static void Compute_Global_Coords( Real node_positions[8][3],
				     Real local_coords[3],
				     Real global_coords[3] );

  static void Interpolate_Scalar( Real  local_coords[3],
				  Real  node_scalars[8],
				  Real& interpolated_scalar );

  static void Interpolate_Vector( Real local_coords[3],
				  Real node_vectors[8][3],
				  Real interpolated_vector[3] );
 private:
  AcmeNode* nodes[8];
  int Node_Ids[16];

};


#endif // ifdef _AcmeHexElemL8_h_

