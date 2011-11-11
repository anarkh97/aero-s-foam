// $Id$

#include "ContactTriFaceL3.C"

template
ContactTriFaceL3<Real>::ContactTriFaceL3( ContactFixedSizeAllocator* alloc,
                                    int Block_Index, 
				    int Index_in_Block,int key );

template
ContactTriFaceL3<Real>* ContactTriFaceL3<Real>::new_ContactTriFaceL3(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key);

template
void ContactTriFaceL3_SizeAllocator<Real>(ContactFixedSizeAllocator& alloc);

template
ContactTriFaceL3<Real>::~ContactTriFaceL3();

template
void ContactTriFaceL3<Real>::Compute_Normal( VariableHandle CURRENT_POSITION,
				       VariableHandle FACE_NORMAL );

template
void ContactTriFaceL3<Real>::Compute_Normal(VariableHandle POSITION,
				       Real* normal, Real* local_coords );
template
void ContactTriFaceL3<Real>::Compute_Normal(Real** nodal_positions,
				      Real* local_coords, 
				      Real* normal );

template
void ContactTriFaceL3<Real>::Compute_CharacteristicLength( VariableHandle CURRENT_POSITION,
						     VariableHandle CHARACTERISTIC_LENGTH );

template
void ContactTriFaceL3<Real>::Compute_Centroid( VariableHandle CURRENT_POSITION,
					 VariableHandle CENTROID );

template
void ContactTriFaceL3<Real>::Get_Edge_Nodes( int i, ContactNode<Real>** node );

template
int ContactTriFaceL3<Real>::Get_Edge_Number( ContactNode<Real>** edge_nodes );

template
int ContactTriFaceL3<Real>::Get_Edge_Number( Real* local_coords );

template
void
ContactTriFaceL3<Real>::Compute_Edge_Normal( VariableHandle POSITION,
				       VariableHandle FACE_NORMAL,
				       int Edge,
				       Real* edge_normal);

template
void ContactTriFaceL3<Real>::Get_Close_Edges( Real* local_coords, int& number,
					int& edge_1_id, int& edge_2_id );

template
bool ContactTriFaceL3<Real>::Is_Inside_Face( Real* local_coords );

template
ContactFace<Real>* ContactTriFaceL3<Real>::Neighbor( Real* local_coords );

template
void ContactTriFaceL3<Real>::FacetDecomposition(int& nfacets,
          Real* coordinates0, Real* normals0, VariableHandle POSITION0,
          Real* coordinates1, Real* normals1, VariableHandle POSITION1,
          Real* coordinates2, Real* normals2, VariableHandle POSITION2);

template
void ContactTriFaceL3<Real>::FacetStaticRestriction(int nfacets, Real* coordinates, 
					      Real* normals, Real* ctrcl_facets,
					      Real* ctrcl);

template
void ContactTriFaceL3<Real>::FacetDynamicRestriction(int nfacets, 
                                               Real* ctrcl_facets, 
                                               Real* ctrcl);

template
void ContactTriFaceL3<Real>::Smooth_Normal( VariableHandle CURRENT_POSITION,
				      VariableHandle NODE_NORMAL,
				      VariableHandle FACE_NORMAL,
				      VariableHandle CURVATURE,
    		       ContactSearch::Smoothing_Resolution resolution,
				      Real percentage,
				      Real* coordinates,
				      Real* smooth_normal,
                                      Real critical_curvature );

template
void ContactTriFaceL3<Real>::Compute_Node_Areas( VariableHandle POSITION,
	                                   VariableHandle /*FACE_NORMAL*/,
                                           Real* node_areas );

template
int ContactTriFaceL3<Real>::FaceEdge_Intersection(VariableHandle POSITION,
                                            ContactEdge<Real>* edge,Real* coords);

template
void ContactTriFaceL3<Real>::Evaluate_Shape_Functions( Real* local_coords,
						 Real* shape_functions );

template
void ContactTriFaceL3<Real>::Compute_Global_Coordinates( VariableHandle POSITION,
						   Real* local_coords,
						   Real* global_coords );

template
void ContactTriFaceL3<Real>::Compute_Local_Coordinates( Real Config_Param,
						  VariableHandle POSITION0,
						  VariableHandle POSITION1,
						  VariableHandle FACE_NORMAL,
						  Real* global_coords,
						  Real* local_coords );

template
void ContactTriFaceL3<Real>::Compute_Local_Coordinates( VariableHandle POSITION,
						  Real* global_coords,
						  Real* local_coords );

template
void ContactTriFaceL3<Real>::Compute_Shape_Functions( Real local_coords[3],
						Real shape_functions[3] );

template
void ContactTriFaceL3<Real>::Compute_Shape_Derivatives( Real local_coords[3],
					          Real shape_derivs[2][3] );

template
void ContactTriFaceL3<Real>::Compute_Local_Coords( Real node_positions[MAX_NODES_PER_FACE][3],
					     Real global_coords[3],
					     Real local_coords[3] );

template
void ContactTriFaceL3<Real>::Compute_Global_Coords( Real node_positions[3][3],
					      Real local_coords[3],
					      Real global_coords[3] );

template
void ContactTriFaceL3<Real>::Interpolate_Scalar( Real  local_coords[3],
					   Real  node_scalars[3],
					   Real& interpolated_scalar );

template
void ContactTriFaceL3<Real>::Interpolate_Vector( Real local_coords[3],
					   Real node_vectors[3][3],
					   Real interpolated_vector[3] );

