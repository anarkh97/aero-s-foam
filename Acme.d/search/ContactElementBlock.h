// $Id$

#ifndef ContactElementBlock_h_
#define ContactElementBlock_h_

#include "ContactElement.h"
#include "ContactBlockEntityList.h"
#include "ContactBoundingBox.h"

class ContactElementBlock {

 public:

  ContactElementBlock( ContactSearch::ContactElement_Type, int, int, int&, int*,
		    ContactTopology* );
  ContactElementBlock( ContactSearch::ContactElement_Type, int,
		    ContactTopology* );
  ~ContactElementBlock();
  
  void Delete_Element_List( );
  
  void Delete_Elements( );
  
  void Add_Elements( int, int* );

  inline ContactSearch::ContactElement_Type Type() { return type;};
  inline int Entity_Key() { return entity_key; };
  inline int ID() { return id; };
  inline int Number_of_Elements() { return number_of_elements; };
  ContactBlockEntityList* ElemList() { return elem_list; };
  ContactBoundingBox* LocalBoundingBox() { return &local_bounding_box; };
  ContactBoundingBox* GlobalBoundingBox() { return &global_bounding_box; };
  void ComputeBoundingBox(int, VariableHandle, VariableHandle, MPI_Comm&);

  void Insert_Element( ContactElement* );
  void Delete_Element( ContactElement* );
#ifndef CONTACT_NO_MPI
  void Insert_Element( char* );
  void Insert_Element_ForSecondary( char* );
#endif

 void SetMaster() {master=true;};
 void SetSlave()  {slave=true;};
 void SetMaster(bool val) {master=val;};
 void SetSlave(bool val)  {slave=val;};
 bool IsMaster()          {return master;};
 bool IsSlave()           {return slave;};

 private:
  
  ContactElementBlock(ContactElementBlock&);
  ContactElementBlock& operator=(ContactElementBlock&);
  
  ContactTopology* topology;

  int number_of_elements;
  int num_elements_added;
  ContactSearch::ContactElement_Type type;
  bool master;
  bool slave;
  int entity_key;
  int id;

  ContactBlockEntityList* elem_list;
  
  ContactBoundingBox local_bounding_box;
  ContactBoundingBox global_bounding_box;

};

#endif //_ContactElementBlock_h_
