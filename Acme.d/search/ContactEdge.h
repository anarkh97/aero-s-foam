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


#ifndef ContactEdge_h_
#define ContactEdge_h_

#include "ContactTopologyEntity.h"
#include "contact_assert.h"
#include "ContactSearch.h"
#include "ContactFace.h"

class ContactNode;

class ContactEdge : public ContactTopologyEntity {

 public:
  
  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#undef  EDGE_SCALAR_VAR
#define EDGE_SCALAR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#undef  EDGE_VECTOR_VAR
#define EDGE_VECTOR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactEdge( ContactSearch::ContactEdge_Type, int, int, ContactNode **nodes);
  virtual ~ContactEdge();

  static void Initialize_Lookup_Arrays();

  int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};
  inline int Nodes_Per_Edge() {
    PRECONDITION(array_init);
    return NODES_PER_EDGE[edge_type];
  };
  inline int Number_Face_Connections() { return number_face_connections; };
  inline int EdgeType() { return(edge_type); };

  inline ContactNode* Node( const int i ) {
    PRECONDITION( i>=0 && i<Nodes_Per_Edge() );
    PRECONDITION( node_list );    
    return node_list[i];
  };
  inline ContactNode** Nodes() {return node_list;};
  inline void ConnectNode(const int num, ContactNode* node ) {
    PRECONDITION( num<Nodes_Per_Edge() );
    PRECONDITION( node_list );
    node_list[num] = node;
  };


  ContactFace* Face( int i );
  ContactFace* Neighbor_Face( ContactFace* face );

  void ConnectFace( ContactFace* );

  // Packing/Unpacking Functions
  inline int  Size();
  void Pack( char* );
  inline void Unpack( char* );
  void Copy( ContactEdge* );

  // Restart Pack/Unpack Functions

  inline int Restart_Size() {return DataArray_Length();}
  inline void Restart_Pack( Real* buffer ) {
    std::memcpy( buffer, DataArray_Buffer(), DataArray_Length()*sizeof(Real) );
  }
  inline void Restart_Unpack( Real* buffer ){
    std::memcpy( DataArray_Buffer(), buffer, DataArray_Length()*sizeof(Real) );
  }

  inline connection_data* FaceInfo() { return &face_info; };

  void Smooth_Normal( VariableHandle, Real* coords, Real* smooth_normal );

  void Compute_Smoothed_Normal( VariableHandle );
  
  inline void Initialize_Memory() {std::memset(DataArray, 0, DataArray_Length()*sizeof(Real));};

 protected:
  int number_face_connections;
  ContactSearch::ContactEdge_Type edge_type;

 private:
  static bool array_init;
  static int                             NODES_PER_EDGE[ContactSearch::NEDGE_TYPES];

  ContactFace* faces[2];    // An edge can only be connected to at most 2 faces

  ContactNode **node_list;

  // The following is used only in the secondary decomposition to store
  // entity information until the actual pointer connections can be made.
  connection_data  face_info;

  Real DataArray[NUMBER_SCALAR_VARS + 3*NUMBER_VECTOR_VARS];
};

inline ContactFace* ContactEdge::Face( int i )
{
  PRECONDITION( i>=0 && i<number_face_connections );
  PRECONDITION( faces );
  return( faces[i] );
}

inline ContactFace* ContactEdge::Neighbor_Face( ContactFace* face )
{
  PRECONDITION( face );
  PRECONDITION( face == faces[0] || face == faces[1] );
  ContactFace* null_neighbor = NULL;
  if( number_face_connections != 2 ) return null_neighbor;
  if( faces[0] == face ) return faces[1];
  if( faces[1] == face ) return faces[0];
  POSTCONDITION( 0 ); // We should never get here
  return null_neighbor;
}

inline int ContactEdge::Size()
{
  return( ContactTopologyEntity::Size(DataArray_Length()) + 
          Nodes_Per_Edge()*sizeof(connection_data) );
}

inline void ContactEdge::Unpack( char* buffer )
{
  ContactTopologyEntity::Unpack( buffer, DataArray_Length() );
  PRECONDITION( ((int*)buffer)[1] == edge_type );
}

#endif // #ifndef ContactEdge_h_
