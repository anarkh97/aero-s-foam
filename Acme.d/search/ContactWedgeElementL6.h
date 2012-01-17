// This file is part of a modified version of ACME: Algorithms for Contact in
// a Multiphysics Environment, derived from version 2.7f
//
// Copyright (C) 2007 Sandia Corporation
// Copyright (C) 2011 Stanford University
//
// ACME is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ACME.  If not, see <http://www.gnu.org/licenses/>.


#ifndef ContactWedgeElementL6_h_
#define ContactWedgeElementL6_h_

#include "ContactElement.h"

class ContactNode;
class ContactEdge;
class ContactFace;
class ContactFixedSizeAllocator;


class ContactWedgeElemL6 : public ContactElem {

 public:
  ContactWedgeElemL6(int blk_indx=-1, int host_indx_in_block=-1, int key=-1);
  static ContactWedgeElemL6* new_ContactWedgeElemL6(ContactFixedSizeAllocator&,
                    int blk_indx=-1, int host_indx_in_block=-1, int key=-1);
  ~ContactWedgeElemL6();
  void BuildTopology(int, int, int, ContactFixedSizeAllocator*);
  void DeleteTopology(ContactFixedSizeAllocator*);
  void UpdateTopology(ContactFace*, VariableHandle, VariableHandle,
                      VariableHandle, Real, bool use_node_normals=false);
  int Nodes_Per_Element() { return 6; };
  int Edges_Per_Element() { return 9; };
  int Faces_Per_Element() { return 5; };
  void Evaluate_Shape_Functions( Real*, Real* );
  void Compute_Local_Coordinates( Real, VariableHandle, VariableHandle,
				  VariableHandle, Real*, Real* );
  void Compute_Local_Coordinates( VariableHandle, Real*, Real* );
  void Compute_Global_Coordinates( VariableHandle, Real*, Real* );
  bool Is_Local_Coordinates_Inside_Element( Real* );
  bool Is_Local_Coordinates_Near_Element( Real*, Real );
  ContactSearch::ContactNode_Type Node_Type() 
    {return ContactSearch::NODE;};
  ContactSearch::ContactEdge_Type Edge_Type() 
    {return ContactSearch::LINEEDGEL2;};
  ContactSearch::ContactFace_Type Face_Type(int i) 
    {return faces[i]->FaceType();};
                      
  static void Compute_Shape_Functions( Real local_coords[4], 
                                       Real Shape_Funcs[6] );
  
  static void Compute_Shape_Derivatives( Real local_coords[4],
                                         Real shape_derivatives[3][6] );
  
  static void Compute_Local_Coords( Real node_positions[6][3],
				    Real global_coords[3],
				    Real local_coords[4] );
  
  static void Compute_Global_Coords( Real node_positions[6][3],
				     Real local_coords[4],
				     Real global_coords[3] );

  static void Interpolate_Scalar( Real  local_coords[4],
				  Real  node_scalars[6],
				  Real& interpolated_scalar );

  static void Interpolate_Vector( Real local_coords[4],
				  Real node_vectors[6][3],
				  Real interpolated_vector[3] );


  inline ContactNode** Nodes() {return nodes;};
  inline ContactEdge** Edges() {return edges;};
  inline ContactFace** Faces() {return faces;};
  inline int* Node_Ids() {return node_ids;};
  inline int* Edge_Ids() {return edge_ids;};
  inline int* Face_Ids() {return face_ids;};

 protected:
 private:
  ContactNode* nodes[6];
  ContactEdge* edges[9];
  ContactFace* faces[5];
  int node_ids[6];
  int edge_ids[9];
  int face_ids[5];

};


#endif // ifdef ContactWedgeElementL6_h_

