// $Id$

// We currently have two kinds of "elements" in ACME std::right now. 
// The first is the "extrusion" of a face to an element (i.e., a T3 to Wedge6
//     or a Q4 to Hex8) used for face-face searches
// The second is an true element that only has node connections (nodes do not
//     have back pointers to the elements std::right now).  These are used for the
//     volume overlap searches.

#ifndef ContactElement_h_
#define ContactElement_h_

#include "ContactEntity.h"
#include "contact_assert.h"
#include "ContactEdge.h"
#include "Contact_Defines.h"
#include "ContactSearch.h"
#include "ContactDoublyLinkedList.h"

class ContactNode;
class ContactEdge;
class ContactFace;
class ContactElementElementInteraction;

class ContactElem : public ContactTopologyEntity {

 public:
  
  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#undef  ELEM_SCALAR_VAR
#define ELEM_SCALAR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#undef  ELEM_VECTOR_VAR
#define ELEM_VECTOR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactElem( ContactSearch::ContactElem_Type, 
	       int Block_Index, int Host_Index_in_Block, int key);
  virtual ~ContactElem();
  
  virtual void BuildTopology(int, int, int, ContactFixedSizeAllocator*) = 0;
  virtual void DeleteTopology(ContactFixedSizeAllocator*) = 0;
  virtual void UpdateTopology(ContactFace*, VariableHandle, VariableHandle,
                              VariableHandle, Real, bool use_node_normals=false) = 0;

  int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};
  virtual int Nodes_Per_Element() = 0;
  virtual int Edges_Per_Element() = 0;
  virtual int Faces_Per_Element() = 0;
  virtual void Evaluate_Shape_Functions( Real*, Real* ) = 0;
  virtual void Compute_Global_Coordinates( VariableHandle, Real*, Real* ) = 0;
  virtual void Compute_Local_Coordinates( VariableHandle, Real*, Real* ) = 0;
  virtual bool Is_Local_Coordinates_Inside_Element( Real* ) = 0;
  virtual bool Is_Local_Coordinates_Near_Element( Real*, Real ) = 0;
  
  virtual ContactSearch::ContactNode_Type Node_Type() = 0;
  virtual ContactSearch::ContactEdge_Type Edge_Type() = 0;
  virtual ContactSearch::ContactFace_Type Face_Type(int) = 0;
  inline  ContactSearch::ContactElem_Type Elem_Type() {return type;};

  inline void Initialize_Memory() {std::memset(DataArray, 0, DataArray_Length()*sizeof(Real));};

  virtual ContactNode** Nodes() = 0;
  virtual ContactEdge** Edges() = 0;
  virtual ContactFace** Faces() = 0;

  ContactNode* Node( int i );
  ContactEdge* Edge( int i );
  ContactFace* Face( int i );
  
  void ConnectNode( int,ContactNode* );
  void ConnectEdge( int,ContactEdge* );
  void ConnectFace( int,ContactFace* );

  Real MaxSize(VariableHandle POSITION);
  // Packing/Unpacking Functions
  int  Size();
  void Pack( char* );
  void Unpack( char* );

  // Restart Pack/Unpack Functions

  inline int Restart_Size() {return DataArray_Length();}
  inline void Restart_Pack( Real* buffer ) {
    std::memcpy( buffer, DataArray_Buffer(), DataArray_Length()*sizeof(Real) );
  }
  inline void Restart_Unpack( Real* buffer ){
    std::memcpy( DataArray_Buffer(), buffer, DataArray_Length()*sizeof(Real) );
  }

  virtual int* Node_Ids() = 0;
  virtual int* Edge_Ids() = 0;
  virtual int* Face_Ids() = 0;

 protected:
  ContactSearch::ContactElem_Type type;

 private:
  Real DataArray[NUMBER_SCALAR_VARS + 3*NUMBER_VECTOR_VARS];
};

inline ContactNode* ContactElem::Node( int i )
{ 
  PRECONDITION( i>=0 && i<Nodes_Per_Element() );
  return( Nodes()[i] );
}

inline ContactEdge* ContactElem::Edge( int i )
{ 
  PRECONDITION( i>=0 && i<Edges_Per_Element() );
  return( Edges()[i] );
}

inline ContactFace* ContactElem::Face( int i )
{ 
  PRECONDITION( i>=0 && i<Faces_Per_Element() );
  return( Faces()[i] );
}




class ContactElement : public ContactTopologyEntity {

 public:

  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#undef  ELEMENT_SCALAR_VAR
#define ELEMENT_SCALAR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#undef  ELEMENT_VECTOR_VAR
#define ELEMENT_VECTOR_VAR( a,b ) a,
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactElement( ContactFixedSizeAllocator*,
                  ContactSearch::ContactElement_Type, 
		  int Block_Index, int Host_Index_in_Block, int key);
  virtual ~ContactElement();
  
  virtual int Nodes_Per_Element() = 0;
  virtual void Compute_Volume(VariableHandle, VariableHandle ) = 0;
  virtual void Compute_Centroid(VariableHandle, VariableHandle ) = 0;
  virtual void TetDice(int &ntets, Real thex[][4][3], VariableHandle) = 0;

  int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};

  virtual void Evaluate_Shape_Functions( Real*, Real* ) = 0;
  virtual void Compute_Global_Coordinates( VariableHandle, Real*, Real* ) = 0;
  virtual void Compute_Local_Coordinates( VariableHandle, Real*, Real* ) = 0;
  virtual bool Is_Local_Coordinates_Inside_Element( Real* ) = 0;
  virtual bool Is_Local_Coordinates_Near_Element( Real*, Real ) = 0;
  virtual void Interpolate_Scalar_Value( Real *, Real *, Real& ) { return; };

  inline ContactSearch::ContactElement_Type ElementType() {return type;};
  virtual ContactNode** Nodes() = 0;
  ContactNode* Node( int i );  
  void ConnectNode( int, ContactNode* );

  // Restart Pack/Unpack Functions

  inline int Restart_Size() {return DataArray_Length();}
  inline void Restart_Pack( Real* buffer ) {
    std::memcpy( buffer, DataArray_Buffer(), DataArray_Length()*sizeof(Real) );
  }
  inline void Restart_Unpack( Real* buffer ){
    std::memcpy( DataArray_Buffer(), buffer, DataArray_Length()*sizeof(Real) );
  }

#ifndef CONTACT_NO_MPI

  Real MaxSize(VariableHandle POSITION);
  
  // Packing/Unpacking Functions
  inline int  Size();
  inline void Pack( char* );
  inline void Unpack( char* );
  inline void Copy( ContactElement* );
  
  inline int  Size_ForSecondary();
  inline void Pack_ForSecondary( char* );
  inline void Unpack_ForSecondary( char* );
  inline void Copy_ForSecondary( ContactElement* );
  
  inline int  Size_ForDataUpdate();
  inline void Pack_ForDataUpdate( char* );
  inline void Unpack_ForDataUpdate( char* );
  
  virtual int  Size_Interactions( int state=0 );
  virtual void Pack_Interactions( char*, int state=0 );
  virtual void Unpack_Interactions( char*, int state=0 );
  virtual void Copy_Interactions( ContactElement*, int state=0 );
  
  virtual int  Size_Interactions_ForSecondary( int state=0 );
  virtual void Pack_Interactions_ForSecondary( char*, int state=0 );
  virtual void Unpack_Interactions_ForSecondary( char*, int state=0 );
  virtual void Copy_Interactions_ForSecondary( ContactElement*, int state=0 );

#endif

  virtual connection_data* NodeInfo() = 0;
  
  inline int Number_Interactions(int state = 0) 
      { int n = ElementElementInteractions[state]->NumEntities();
        return n; };
  
  inline int Number_ElementElement_Interactions (int state = 0) 
      { return ElementElementInteractions[state]->NumEntities(); };
                                         
  inline ContactInteractionDLL* Get_ElementElement_Interactions( int state = 0 )
       { return ElementElementInteractions[state]; };
      
  ContactElementElementInteraction* Get_ElementElement_Interaction(int interaction_number,
						                   int state = 0 );
       
  void Store_ElementElement_Interaction( ContactElementElementInteraction*, 
					 int state = 0 );
                                         
  void Delete_ElementElement_Interaction( ContactElementElementInteraction*, 
			 		  int state = 0 );
                                         
  void Display_ElementElement_Interactions( ContactParOStream&, int state = 0 );
  
  void Update_Interactions( );

  inline void Initialize_Memory() {std::memset(DataArray, 0, DataArray_Length()*sizeof(Real));};
  
  void ComputeBoundingBoxForSearch(const int num_configs,
                                   const VariableHandle &NODE_COORD_START,
                                   const VariableHandle &NODE_COORD_END,
                                   const int  auto_tol,
                                   const Real box_inflation,
                                   const Real user_tol,
                                   ContactBoundingBox &box);

 protected:
  ContactSearch::ContactElement_Type type;
  ContactFixedSizeAllocator* allocators;

 private:
  Real DataArray[NUMBER_SCALAR_VARS + 3*NUMBER_VECTOR_VARS];
  
  int number_of_states;
  ContactInteractionDLL** ElementElementInteractions;
};

inline ContactNode* ContactElement::Node( int i )
{ 
  PRECONDITION( i>=0 && i<Nodes_Per_Element() );
  return( Nodes()[i] );
}

inline void ContactElement::ConnectNode( int i, ContactNode* node ) {
  PRECONDITION( i>=0 && i<Nodes_Per_Element() );
  Nodes()[i] = node;
  PackConnection((ContactTopologyEntity*)node, &NodeInfo()[i]);
}

#ifndef CONTACT_NO_MPI
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used 
// to just update the data for ghosted elements
//--------------------------------------------------------------------
inline int ContactElement::Size_ForDataUpdate()
{
  return( ContactTopologyEntity::Size_ForDataUpdate(DataArray_Length()) );
}

inline void ContactElement::Pack_ForDataUpdate( char* buffer )
{
  ContactTopologyEntity::Pack_ForDataUpdate( buffer, DataArray_Length() );
}

inline void ContactElement::Unpack_ForDataUpdate( char* buffer )
{
  ContactTopologyEntity::Unpack_ForDataUpdate( buffer, DataArray_Length() );
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for DLB
//--------------------------------------------------------------------
inline int ContactElement::Size()
{
  return( ContactTopologyEntity::Size(DataArray_Length()) + 
          Nodes_Per_Element()*sizeof(connection_data) );
}

inline void ContactElement::Pack( char* buffer )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  // ContactTopologyEntity packs in location 0 as ContactElement 
  // and here we pack in the derived type in location 1.
  i_buf[1] = type;
  ContactTopologyEntity::Pack( buffer, DataArray_Length() );
  // Add the global ids of the nodes
  i_buf = reinterpret_cast<int*> (buffer+ContactTopologyEntity::Size(DataArray_Length()));
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Element() ; ++i ){
    cnt += PackConnection(reinterpret_cast<ContactTopologyEntity*>(Node(i)), &i_buf[cnt]);
  }
}

inline void ContactElement::Unpack( char* buffer )
{
  ContactTopologyEntity::Unpack( buffer, DataArray_Length() );
  entity_key = block_id;

  PRECONDITION( ((int*)buffer)[1] == type );
  // Store off the global ids of the nodes
  int* i_buf = reinterpret_cast<int*> ( buffer + ContactTopologyEntity::Size(DataArray_Length()) );
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Element() ; ++i ){
    cnt += UnPackConnection(&NodeInfo()[i], &i_buf[cnt]);
  }
}

inline void ContactElement::Copy( ContactElement* src )
{
  ContactTopologyEntity::Copy( src, DataArray_Length() );
  entity_key = src->entity_key;
  for( int i=0 ; i<Nodes_Per_Element() ; ++i ){
    PackConnection(reinterpret_cast<ContactTopologyEntity*>(src->Node(i)), &NodeInfo()[i]);
  }
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for 
// transferring entities from the primary to secondary decomposition
//--------------------------------------------------------------------
inline int ContactElement::Size_ForSecondary()
{
  return( ContactTopologyEntity::Size_ForSecondary(DataArray_Length()) + 
          Nodes_Per_Element()*sizeof(connection_data) );
}

inline void ContactElement::Pack_ForSecondary( char* buffer )
{
  int* i_buf = reinterpret_cast<int*> (buffer);
  // ContactTopologyEntity packs in location 0 as ContactElement 
  // and here we pack in the derived type in location 1.
  i_buf[1] = type;
  ContactTopologyEntity::Pack_ForSecondary( buffer, DataArray_Length() );
  // Add the global ids of the nodes
  i_buf = reinterpret_cast<int*> (buffer+ContactTopologyEntity::Size_ForSecondary(DataArray_Length()));
  int cnt = 0;
  for( int i=0 ; i<Nodes_Per_Element() ; ++i ){
    cnt += PackConnection(reinterpret_cast<ContactTopologyEntity*>(Node(i)), &i_buf[cnt]);
  }
}

inline void ContactElement::Unpack_ForSecondary( char* buffer )
{
  ContactTopologyEntity::Unpack_ForSecondary( buffer, DataArray_Length() );
  entity_key = block_id;

  PRECONDITION( ((int*)buffer)[1] == type );
  // Store off the global ids of the nodes
  int* i_buf = reinterpret_cast<int*> ( buffer + ContactTopologyEntity::Size_ForSecondary(DataArray_Length()) );
  int cnt = 0;
  for( int i=0; i<Nodes_Per_Element(); ++i ){
    cnt += UnPackConnection(&NodeInfo()[i], &i_buf[cnt]);
  }
}

inline void ContactElement::Copy_ForSecondary( ContactElement* src )
{
  ContactTopologyEntity::Copy_ForSecondary( src, DataArray_Length() );
  entity_key = src->entity_key;
  for( int i=0; i<Nodes_Per_Element(); ++i ){
    PackConnection(reinterpret_cast<ContactTopologyEntity*>(src->Node(i)), &NodeInfo()[i]);
  }
}
#endif

#endif // #ifdef ContactElement_h_
