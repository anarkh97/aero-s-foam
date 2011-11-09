// $Id: AcmeEntityBlock.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeEntityBlock_h_
#define _AcmeEntityBlock_h_

#include "AcmeEntity.h"
#include "AcmeBlockEntityList.h"
#include "ContactSearch.h"

class AcmeEntityBlock {

 public:

   AcmeEntityBlock(AcmeEntity::BaseEntityType, 
                   AcmeEntity::DerivedEntityType, 
                   int, int, int*);
  ~AcmeEntityBlock();

  void UpdateEntityBlock( int, int, int& );

  inline AcmeEntity::BaseEntityType BaseType() { return base_type;};
  inline AcmeEntity::DerivedEntityType DerivedType() { return derived_type;};
  
  inline void NodeType(ContactSearch::ContactNode_Type n) { node_type=n; };
  inline void FaceType(ContactSearch::ContactFace_Type f) { face_type=f; };
  inline void ElemType(ContactSearch::ContactElement_Type e) { elem_type=e; };
  
  inline ContactSearch::ContactNode_Type    NodeType() { return node_type; };
  inline ContactSearch::ContactFace_Type    FaceType() { return face_type; };
  inline ContactSearch::ContactElement_Type ElemType() { return elem_type; };
  
  inline void NumAttributes(int n) { num_attributes=n;};
  inline int NumAttributes() { return num_attributes;};
  inline void ShellOffset(Real s) { shell_layer_offset=s;};
  inline Real ShellOffset() { return shell_layer_offset;};
  inline void IsShell(bool s) { is_a_shell=s;};
  inline bool IsShell() { return is_a_shell;};
  inline int ID() { return id; };
  inline int Number_of_Entities() { return number_of_entities; };
  AcmeBlockEntityList* EntityList() { return entity_list; };
  
  void Delete_Entity_List( );
  void Delete_Entities( );
  void Add_Entities( int, int* );
  void Insert_Entity( AcmeEntity* );
  void Delete_Entity( AcmeEntity* );
#ifndef CONTACT_NO_MPI
  void Insert_Entity( char* );
#endif

 private:

  int id;
  AcmeEntity::BaseEntityType    base_type;
  AcmeEntity::DerivedEntityType derived_type;
  ContactSearch::ContactNode_Type    node_type;
  ContactSearch::ContactFace_Type    face_type;
  ContactSearch::ContactElement_Type elem_type;
  
  int num_attributes;
  int number_of_entities;
  AcmeBlockEntityList* entity_list;

  bool is_a_shell;
  Real shell_layer_offset;
  int num_entities_added; // This is used to place new nodes in Add_Node

};

#endif  // #ifdef _AcmeEntityBlock_h_
