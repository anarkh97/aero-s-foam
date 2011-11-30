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


#ifndef ContactFace_h_
#define ContactFace_h_

#include "ContactTopologyEntity.h"
#include "contact_assert.h"
#include "ContactEdge.h"
#include "Contact_Defines.h"
#include "ContactSearch.h"
#include "ContactDoublyLinkedList.h"
#include "ContactNode.h"

class ContactNode;
class ContactEdge;
class ContactFaceFaceInteraction;
class ContactFaceFaceInteraction;
class ContactFaceCoverageInteraction;

class ContactFace : public ContactTopologyEntity {

 public:

  enum Ctrcl_Index{MSPARAM = 0,  ICPOINTX = 1,  ICPOINTY = 2,  ICPOINTZ  = 3,
                   IPENMAG = 4,  IPUSHX   = 5,  IPUSHY   = 6,  IPUSHZ    = 7,
                   INORMX  = 8,  INORMY   = 9,  INORMZ   = 10, ILOCATION = 11,
                   ICTIMC  = 12, LENGTH  = 13};
  
  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#undef  FACE_SCALAR_VAR
#define FACE_SCALAR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#undef  FACE_VECTOR_VAR
#define FACE_VECTOR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactFace( ContactFixedSizeAllocator*,
               ContactSearch::ContactFace_Type, 
               int Block_Index,
	       int Host_Index_in_Block, 
               int key,
               ContactNode **node_list_,
               ContactEdge **edge_list_,
               connection_data *node_info_list_,
               connection_data *edge_info_list_);

  virtual ~ContactFace();

  static inline int Nodes_Per_Face(ContactSearch::ContactFace_Type check_face_type) {
    PRECONDITION(array_init);
    return NODES_PER_FACE[check_face_type];
  }



  inline int Nodes_Per_Face() {
    PRECONDITION(array_init);
    return NODES_PER_FACE[face_type];
  };

  inline int Edges_Per_Face() {
    PRECONDITION(array_init);
    return EDGES_PER_FACE[face_type];
  };

  inline ContactNode**     Nodes()    {return node_list;};
  inline connection_data*  NodeInfo() {return node_info_list;};
  inline ContactEdge**     Edges()    {return edge_list;};
  inline connection_data*  EdgeInfo() {return edge_info_list;};


  virtual ContactSearch::ContactEdge_Type Edge_Type() = 0;
  inline ContactNode* Node( const int i );
  static void Initialize_Lookup_Arrays();


  // This function get the two nodes that terminate edge i
  virtual void Get_Edge_Nodes( int, ContactNode**) = 0;
  virtual int  Get_Edge_Number( ContactNode** ) = 0;
  virtual int  Get_Edge_Number( Real* ) = 0;
  virtual void Compute_Normal(VariableHandle, VariableHandle ) = 0;
  virtual void Compute_Normal(VariableHandle, Real*, Real* ) = 0;
  virtual void Compute_Normal(Real**, Real*, Real* ) = 0;
  virtual void Compute_CharacteristicLength(VariableHandle, VariableHandle ) = 0;
  virtual void Compute_Centroid( VariableHandle, VariableHandle ) = 0;
  virtual void Compute_Edge_Normal( VariableHandle, VariableHandle,
				    int , Real*) = 0;
  virtual void Compute_Local_Coordinates( Real, VariableHandle, 
                                          VariableHandle, VariableHandle, 
		                          Real*, Real* ) = 0;

  virtual void Compute_Local_Coords(Real node_positions[MAX_NODES_PER_FACE][3], 
				    Real global_coords[3],
			            Real local_coords[3]) = 0;


  virtual void Compute_Local_Coordinates( VariableHandle, Real*, Real* ) = 0;
  virtual void Compute_Global_Coordinates( VariableHandle, Real*, Real* ) = 0;
  virtual void Evaluate_Shape_Functions( Real* local_coord, 
					 Real* shape_fnc) = 0;
  virtual bool Is_Inside_Face( Real* Local_Coordinates ) = 0;
  virtual ContactFace* Neighbor( Real* Local_Coordinates ) = 0;

  // this is provided for the curvature modification 
  virtual void Get_Close_Edges( Real*, int&, int&, int& ) = 0;
  inline void Initialize_Memory() {std::memset(DataArray, 0, DataArray_Length()*sizeof(Real));};

  virtual void FacetDecomposition(int& nfacets, 
                                  Real* coords0       , Real* normal0       , 
				  VariableHandle POSITION0,
                                  Real* coords1 = NULL, Real* normal1 = NULL, 
				  VariableHandle POSITION1 = 0,
                                  Real* coords2 = NULL, Real* normal2 = NULL, 
				  VariableHandle POSITION2 = 0) = 0;
  virtual void FacetStaticRestriction(int, Real*, Real*, Real*, Real*) = 0;
  virtual void FacetDynamicRestriction(int, Real*, Real*) = 0;
  virtual int  FaceEdge_Intersection(VariableHandle, ContactEdge*, Real*) = 0;
  virtual bool IsPlanar(VariableHandle) = 0;
 
#ifndef CONTACT_NO_MPI
  Real MaxSize(VariableHandle POSITION);
#endif

  int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};
  inline ContactEdge* Edge( const int i );
  void ConnectNode( const int, ContactNode* );
  void ConnectEdge( const int, ContactEdge* );

  inline ContactEdge* Clockwise_Edge( ContactEdge* );
  inline void Clockwise_EdgeNode( int, ContactNode** );

  inline ContactSearch::ContactFace_Type FaceType() {return face_type;};

  // Packing/Unpacking Functions
  inline int  Size(int flag=0);
  inline void Pack( char*, int flag=0 );
  inline void Unpack( char* );
  inline void Copy( ContactFace* face, int neighbors=0);
  
  inline int  Size_ForSecondary(int include_neighbors=0, int include_edgeinfo=0);
  inline void Pack_ForSecondary( char*, int include_neighbors=0, int include_edgeinfo=0 );
  inline void Unpack_ForSecondary( char* );
  inline void Copy_ForSecondary( ContactFace* face, int include_neighbors=0, int include_edgeinfo=0);
  
  inline int  Size_ForDataUpdate();
  inline void Pack_ForDataUpdate( char* );
  inline void Unpack_ForDataUpdate( char* );
  
  int  Size_Interactions( int state=0 );
  void Pack_Interactions( char*, int state=0 );
  void Unpack_Interactions( char*, int state=0 );
  void Copy_Interactions( ContactFace*, int state=0 );
  
  int  Size_Interactions_ForSecondary( int state=0 );
  void Pack_Interactions_ForSecondary( char*, int state=0 ); 
  void Unpack_Interactions_ForSecondary( char*, int state=0 );
  void Copy_Interactions_ForSecondary( ContactFace*, int state=0 );
  
  // Restart Pack/Unpack Functions

  inline int Restart_Size() {return DataArray_Length();}
  inline void Restart_Pack( Real* buffer ) {
    std::memcpy( buffer, DataArray_Buffer(), DataArray_Length()*sizeof(Real) );
  }
  inline void Restart_Unpack( Real* buffer ){
    std::memcpy( DataArray_Buffer(), buffer, DataArray_Length()*sizeof(Real) );
  }

  virtual void Smooth_Normal( VariableHandle, VariableHandle, VariableHandle,
			      VariableHandle, 
			      ContactSearch::Smoothing_Resolution,
			      Real, Real*, Real*, Real ) = 0;
  
  virtual void Compute_Node_Areas( VariableHandle, VariableHandle, Real* ) = 0;

  inline int Number_Interactions(int state = 0) {
    int num_ffi = 0;
    if((int)FaceFaceInteractions.size() > state) num_ffi = FaceFaceInteractions[state].NumEntities();
    int num_fci = 0;
    if((int)FaceCoverageInteractions.size() > state) num_fci = FaceCoverageInteractions[state].NumEntities();
    return num_ffi + num_fci;
  };
  
  inline int Number_FaceFace_Interactions (int state = 0) { 
    if((int)FaceFaceInteractions.size() <= state) return 0;
    return FaceFaceInteractions[state].NumEntities(); 
  };
                                         
  inline ContactInteractionDLL* Get_FaceFace_Interactions( int state = 0 ) { 
    if((int)FaceFaceInteractions.size() <= state) return NULL;
    return &(FaceFaceInteractions[state]); 
  };
      
  ContactFaceFaceInteraction* Get_FaceFace_Interaction(int interaction_number,
						       int state = 0 );
                                         
  void Store_FaceFace_Interaction( ContactFaceFaceInteraction*, 
                                   int state = 0 );
                                         
  void Delete_FaceFace_Interaction( ContactFaceFaceInteraction*, 
                                    int state = 0 );
                                   
  void Display_FaceFace_Interactions( ContactParOStream&, int state = 0 );
  
  inline int Number_FaceCoverage_Interactions (int state = 0) { 
    if((int)FaceCoverageInteractions.size() <= state) return 0;
    return FaceCoverageInteractions[state].NumEntities(); 
  };
                                         
  inline ContactInteractionDLL* Get_FaceCoverage_Interactions( int state = 0 ) { 
    if((int)FaceCoverageInteractions.size() <= state) return NULL;
    return &(FaceCoverageInteractions[state]); 
  };
                                         
  void Store_FaceCoverage_Interaction( ContactFaceCoverageInteraction*, 
                                       int state = 0 );
                                       
  void Display_FaceCoverage_Interactions( ContactParOStream&, int state = 0 );
  
  void Update_Interactions( );
  
  void SetNeighborFacesInfo( );
  
  int NumberOfNeighbors() { return number_of_neighbors; };
  connection_data* NeighborInfo() { return neighbor_face_info; };

  void SetEdgeCurvature(VariableHandle);

  void SetEdgeCurvature(VariableHandle &var, ContactEdge *edge);


  Real GetEdgeCurvature(int);
  
  void GetEdgeInfo(ContactNode* node, ContactNode** edge_nodes, 
                   int* edge_nums);

  void SetEdgeSmoothedNormal(VariableHandle);
  void GetEdgeSmoothedNormal(int, Real*);
  
  void ComputeBoundingBoxForSearch(const int num_configs,
                                   const VariableHandle &NODE_COORD_START,
                                   const VariableHandle &NODE_COORD_END,
                                   const int  auto_tol,
                                   const Real box_inflation,
                                   const Real user_tol,
                                   ContactBoundingBox &box_c,
                                   ContactBoundingBox &box_p,
                                   ContactBoundingBox &box_s);
  
  void ComputeBoundingBoxForSearch(const int num_configs,
                                   const VariableHandle &NODE_COORD_START,
                                   const VariableHandle &NODE_COORD_END,
                                   const int  auto_tol,
                                   const Real box_inflation,
                                   const Real* max_node_motion,
                                   const Real max_remaining_gap_mag,
                                   const Real user_search_tol,
                                   ContactBoundingBox &box_c,
                                   ContactBoundingBox &box_p,
                                   ContactBoundingBox &box_s);

  

 protected:
  int number_of_neighbors;
  ContactSearch::ContactFace_Type face_type;
  ContactFixedSizeAllocator* allocators;

 private:
  //
  //  Arrays for fast lookup of face info based on face type
  //
  static bool array_init;
  static int  NODES_PER_FACE[ContactSearch::NFACE_TYPES];
  static int  EDGES_PER_FACE[ContactSearch::NFACE_TYPES];

  // The edges and node arrays are actually owned in the derived class
  // but we hold a pointer to them to provide access through the base class.
  
  ContactNode **node_list;
  ContactEdge **edge_list;
  connection_data *node_info_list;   
  connection_data *edge_info_list;

  connection_data* neighbor_face_info;
  
//  int entity_key;

  Real DataArray[NUMBER_SCALAR_VARS + 3*NUMBER_VECTOR_VARS];
  
  inline int NumberOfStates() const {return 2;};

  std::vector<ContactInteractionDLL> FaceFaceInteractions;
  std::vector<ContactInteractionDLL> FaceCoverageInteractions;
};

inline ContactNode* ContactFace::Node( const int i ) { 
  PRECONDITION( i>=0 && i<Nodes_Per_Face() );
  return( Nodes()[i] );
}

inline ContactEdge* ContactFace::Edge( const int i ) { 
  PRECONDITION( i>=0 && i<Edges_Per_Face() );
  return( Edges()[i] );
}

inline ContactEdge* ContactFace::Clockwise_Edge( ContactEdge* edge ) {
  int i;
  // Find this edge in the edge list
  for( i=0 ; i<Edges_Per_Face() ; ++i){
    if( Edges()[i] == edge ){
      if( i == Edges_Per_Face()-1 )
	return Edges()[0];
      else 
	return Edges()[i+1];
    }
  }
  // We didn't find this edge which is an error
  POSTCONDITION( 0 );
  return NULL;
}

inline void ContactFace::Clockwise_EdgeNode( int edge, ContactNode** edge_nodes)
{
  int n1 = edge==Edges_Per_Face()-1?0:edge+1;
  int n2 =   n1==Edges_Per_Face()-1?0:n1+1;
  edge_nodes[0] = node_list[n1];
  edge_nodes[1] = node_list[n2];
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for DLB
//--------------------------------------------------------------------
inline int ContactFace::Size(int include_neighbors)
{
  if (include_neighbors) {
    return( ContactTopologyEntity::Size(DataArray_Length()) + sizeof(int) +
	    Nodes_Per_Face()*sizeof(connection_data) +
	    Edges_Per_Face()*sizeof(connection_data) +
            Edges_Per_Face()*sizeof(connection_data));
  } else {
    return( ContactTopologyEntity::Size(DataArray_Length()) + sizeof(int) +
            Nodes_Per_Face()*sizeof(connection_data) +
	    Edges_Per_Face()*sizeof(connection_data));
  }
}

inline void ContactFace::Pack( char* buffer, int include_neighbors )
{
  int* i_buf = reinterpret_cast<int*>(buffer);
  // ContactTopologyEntity packs in location 0 as ContactFace and here we pack
  // in the derived type in location 1.
  i_buf[1] = face_type;
  ContactTopologyEntity::Pack( buffer, DataArray_Length() );
  // Add the entity data for the nodes and edges
  i_buf = reinterpret_cast<int*>(buffer+ContactTopologyEntity::Size(DataArray_Length()));
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    cnt += PackConnection(reinterpret_cast<ContactTopologyEntity*>(Node(i)), &i_buf[cnt]);
  }
  for( int i=0 ; i<Edges_Per_Face() ; ++i ){
    cnt += PackConnection(reinterpret_cast<ContactTopologyEntity*>(Edge(i)), &i_buf[cnt]);
  }
  if (include_neighbors) {
    i_buf[cnt++] = number_of_neighbors;
    for( int i=0 ; i<Edges_Per_Face() ; ++i ){
      cnt += PackConnection(&neighbor_face_info[i], &i_buf[cnt]);
    }
  } else {
    i_buf[cnt] = 0;
  }
}

inline void ContactFace::Unpack( char* buffer )
{
  ContactTopologyEntity::Unpack( buffer, DataArray_Length() );
  entity_key = block_id;

  PRECONDITION( ((int*)buffer)[1] == face_type );
  // Store off the entity data for the nodes and edge
  int* i_buf = reinterpret_cast<int*> ( buffer + ContactTopologyEntity::Size(DataArray_Length()) );
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    cnt += UnPackConnection(&NodeInfo()[i], &i_buf[cnt]);
  }
  for( int i=0 ; i<Edges_Per_Face() ; ++i ){
    cnt += UnPackConnection(&EdgeInfo()[i], &i_buf[cnt]);
  }
  number_of_neighbors = i_buf[cnt++];
  if (number_of_neighbors>0) {
    neighbor_face_info = new connection_data[Edges_Per_Face()];
    for( int i=0 ; i<Edges_Per_Face() ; ++i ){
      cnt += UnPackConnection(&neighbor_face_info[i], &i_buf[cnt]);
    }
  }
}

inline void ContactFace::Copy( ContactFace* src, int include_neighbors )
{
  ContactTopologyEntity::Copy( src, DataArray_Length() );
  entity_key = src->entity_key;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    PackConnection(reinterpret_cast<ContactTopologyEntity*>(src->Node(i)), &NodeInfo()[i]);
  }
  for( int i=0 ; i<Edges_Per_Face() ; ++i ){
    PackConnection(reinterpret_cast<ContactTopologyEntity*>(src->Edge(i)), &EdgeInfo()[i]);
  }
  if (include_neighbors) {
    if (neighbor_face_info) delete [] neighbor_face_info;
    neighbor_face_info  = NULL;
    number_of_neighbors = src->number_of_neighbors;
    if (number_of_neighbors>0) {
      neighbor_face_info = new connection_data[Edges_Per_Face()];
      for (int i=0; i<Edges_Per_Face(); ++i) {
        neighbor_face_info[i] = src->neighbor_face_info[i];
      }
    }
  } else {
    if (neighbor_face_info) delete [] neighbor_face_info;
    neighbor_face_info  = NULL;
    number_of_neighbors = 0;
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for 
// transferring entities from the primary to secondary decomposition
//--------------------------------------------------------------------
inline int ContactFace::Size_ForSecondary(int include_neighbors, 
                                          int include_edgeinfo)
{
  int cnt = sizeof(int) + Nodes_Per_Face()*sizeof(connection_data);
  if (include_edgeinfo) {
    cnt += ContactTopologyEntity::Size_ForSecondary(DataArray_Length());
    cnt += Edges_Per_Face()*sizeof(connection_data);
  } else {
    cnt += ContactTopologyEntity::Size_ForSecondary(DataArray_Length());
  }
  if (include_neighbors) {
    cnt += number_of_neighbors*(sizeof(connection_data)+sizeof(int));
  }
  return cnt;
}

inline void ContactFace::Pack_ForSecondary( char* buffer, 
                                            int include_neighbors, 
                                            int include_edgeinfo )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  // ContactTopologyEntity packs in location 0 as ContactFace and here we pack
  // in the derived type in location 1.
  i_buf[1] = face_type;
  
  if (include_edgeinfo) {
    ContactTopologyEntity::Pack_ForSecondary( buffer, DataArray_Length() );
  } else {
    ContactTopologyEntity::Pack_ForSecondary( buffer, &DataArray[Edge0_Curvature], DataArray_Length() );
  }
  i_buf[SEC_OWNER] |= include_edgeinfo<<24;
  i_buf = reinterpret_cast<int*>(buffer+ContactTopologyEntity::Size_ForSecondary(DataArray_Length()));
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    cnt += PackConnection(reinterpret_cast<ContactTopologyEntity*>(Node(i)), &i_buf[cnt]);
  }
  if (include_edgeinfo) {
    for( int i=0 ; i<Edges_Per_Face() ; ++i ){
      cnt += PackConnection(reinterpret_cast<ContactTopologyEntity*>(Edge(i)), &i_buf[cnt]);
    }
  }
  if (include_neighbors) {
    int nn=0;
    i_buf[cnt++] = number_of_neighbors;
    for( int i=0 ; i<Edges_Per_Face() ; ++i ){
      if (neighbor_face_info[i].owner>=0) {
        ++nn;
        i_buf[cnt++] = i;
        cnt += PackConnection(&neighbor_face_info[i], &i_buf[cnt]);
      }
    }
    POSTCONDITION(nn==number_of_neighbors);
  } else {
    i_buf[cnt++] = 0;
  }
}

inline void ContactFace::Unpack_ForSecondary( char* buffer )
{
  int* i_buf = reinterpret_cast<int*>( buffer );
  PRECONDITION( i_buf[1] == face_type );
  int include_edgeinfo = i_buf[SEC_OWNER]>>24;
  i_buf[SEC_OWNER] &= 0xFFFFFF;
  if (include_edgeinfo) {
    ContactTopologyEntity::Unpack_ForSecondary( buffer, DataArray_Length() );
  } else {
    ContactTopologyEntity::Unpack_ForSecondary( buffer, &DataArray[Edge0_Curvature], DataArray_Length() );
  }
  entity_key = block_id;

  // Store off the entity data for the nodes and edge
  i_buf = reinterpret_cast<int*>( buffer + ContactTopologyEntity::Size_ForSecondary(DataArray_Length()) );
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i){
    cnt += UnPackConnection(&NodeInfo()[i], &i_buf[cnt]);
  }
  if (include_edgeinfo) {
    for( int i=0 ; i<Edges_Per_Face() ; ++i){
      cnt += UnPackConnection(&EdgeInfo()[i], &i_buf[cnt]);
    }
  }
  neighbor_face_info = new connection_data[Edges_Per_Face()];
  for (int i=0; i<Edges_Per_Face(); ++i) {
    neighbor_face_info[i].owner = -1;;
  }
  number_of_neighbors = i_buf[cnt++];
  if (number_of_neighbors>0) {
    int nn = 0;
    while (nn<number_of_neighbors) {
      int ii = i_buf[cnt++];
      cnt   += UnPackConnection(&neighbor_face_info[ii], &i_buf[cnt]);
      ++nn;
    }
  }
}

inline void ContactFace::Copy_ForSecondary( ContactFace* src, 
                                            int include_neighbors, 
                                            int include_edgeinfo )
{
  ContactTopologyEntity::Copy_ForSecondary( src, DataArray_Length() );
  entity_key = src->entity_key;
  for( int i=0 ; i<Nodes_Per_Face() ; ++i ){
    PackConnection(reinterpret_cast<ContactTopologyEntity*>(src->Node(i)), &NodeInfo()[i]);
  }
  if (include_edgeinfo) {
    for( int i=0 ; i<Edges_Per_Face() ; ++i ){
      PackConnection(reinterpret_cast<ContactTopologyEntity*>(src->Edge(i)), &EdgeInfo()[i]);
    }
  }
  if (include_neighbors) {
    if (neighbor_face_info) delete [] neighbor_face_info;
    neighbor_face_info  = NULL;
    number_of_neighbors = src->number_of_neighbors;
    if (number_of_neighbors>0) {
      neighbor_face_info = new connection_data[Edges_Per_Face()];
      for (int i=0; i<Edges_Per_Face(); ++i) {
        neighbor_face_info[i] = src->neighbor_face_info[i];
      }
    }
  } else {
    if (neighbor_face_info) delete [] neighbor_face_info;
    neighbor_face_info = new connection_data[Edges_Per_Face()];
    for (int i=0; i<Edges_Per_Face(); ++i) {
      neighbor_face_info[i].owner = -1;;
    }
    number_of_neighbors = 0;
  }
}

inline int ContactFace::Size_ForDataUpdate()
{
  return( ContactTopologyEntity::Size_ForDataUpdate(DataArray_Length()) );
}

inline void ContactFace::Pack_ForDataUpdate( char* buffer )
{
  ContactTopologyEntity::Pack_ForDataUpdate( buffer, DataArray_Length() );
}

inline void ContactFace::Unpack_ForDataUpdate( char* buffer )
{
  ContactTopologyEntity::Unpack_ForDataUpdate( buffer, DataArray_Length() );
}

#endif // #ifdef ContactFace_h_

