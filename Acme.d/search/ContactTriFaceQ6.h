// $Id$

#ifndef ContactTriFaceQ6_h_
#define ContactTriFaceQ6_h_

#include "ContactFace.h"
#include "ContactEdge.h"

class ContactFixedSizeAllocator;
class ContactNode;
class ContactEdge;

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

class ContactTriFaceQ6 : public ContactFace {
 public:
  ContactTriFaceQ6( ContactFixedSizeAllocator*, int blk_indx=-1, 
                    int indx_in_block=-1, int key=-1 );
  static ContactTriFaceQ6* new_ContactTriFaceQ6(ContactFixedSizeAllocator*,
                                                int blk_indx=-1, 
                                                int indx_in_block=-1, 
                                                int key=-1);
  ~ContactTriFaceQ6( );
  ContactSearch::ContactEdge_Type Edge_Type() 
    {return ContactSearch::LINEEDGEQ3;};
  void Get_Edge_Nodes( int, ContactNode**);
  int Get_Edge_Number( ContactNode** );
  int Get_Edge_Number( Real* );

  void Compute_Normal(VariableHandle, VariableHandle );
  void Compute_Normal(VariableHandle, Real*, Real* );
  void Compute_Normal(Real**, Real*, Real* );
  void Compute_CharacteristicLength(VariableHandle, VariableHandle );
  void Compute_Centroid(VariableHandle, VariableHandle );
  void Compute_Edge_Normal( VariableHandle, VariableHandle,
					int , Real*);
  void Compute_Local_Coordinates( Real, VariableHandle, VariableHandle,
				  VariableHandle, Real*, Real* );
  void Compute_Local_Coordinates( VariableHandle, Real*, Real* );
  void Compute_Global_Coordinates( VariableHandle, Real*, Real* );
  void Evaluate_Shape_Functions( Real* local_coords, Real* shape_funcs );
  bool Is_Inside_Face( Real* local_coords );
  inline bool IsPlanar(VariableHandle) {return false;};
  ContactFace* Neighbor( Real* local_coords );
  void Get_Close_Edges( Real*, int&, int&, int& );
  void FacetDecomposition(int &, 
                          Real*, Real*, VariableHandle,
                          Real*, Real*, VariableHandle,
                          Real*, Real*, VariableHandle);
  void FacetStaticRestriction(int, Real*, Real*, Real*, Real*);
  void FacetDynamicRestriction(int, Real*, Real*);

  void Smooth_Normal( VariableHandle, VariableHandle, VariableHandle, 
		      VariableHandle, ContactSearch::Smoothing_Resolution,
		      Real, Real*, Real*, Real );

  void Compute_Node_Areas( VariableHandle, VariableHandle, Real* );
                      
  int FaceEdge_Intersection(VariableHandle, ContactEdge*, Real*);
  
  static void Compute_Shape_Functions( Real local_coords[3], 
                                       Real shape_funcs[6] );
  
  static void Compute_Shape_Derivatives( Real local_coords[3],
                                         Real shape_derivatives[2][6] );
                                   
  void Compute_Local_Coords( Real node_positions[MAX_NODES_PER_FACE][3],
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
  ContactNode* nodes[6];
  ContactEdge* edges[3];
  connection_data Node_Info[6];
  connection_data Edge_Info[3];
};




#endif // ifndef ContactTriFaceQ6_h_
