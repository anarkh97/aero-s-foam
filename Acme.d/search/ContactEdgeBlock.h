// $Id$

#ifndef ContactEdgeBlock_h_
#define ContactEdgeBlock_h_

#include "ContactEdge.h"
#include "ContactBlockEntityList.h"

class ContactEdgeBlock {

 public:

  ContactEdgeBlock( ContactSearch::ContactEdge_Type, int, int, int, int&,
		    ContactTopology* );
  ContactEdgeBlock( ContactSearch::ContactEdge_Type, int,
		    ContactTopology* );
  ~ContactEdgeBlock();
  
  void Delete_Edge_List( );
  
  void Delete_Edges( );

  inline ContactSearch::ContactEdge_Type Type() { return type;};
  inline int Entity_Key() { return entity_key; };
  inline int ID() { return id; };
  inline int Number_of_Edges() { return number_of_edges; };
  ContactBlockEntityList* EdgeList() { return edge_list; };

  void Insert_Edge( ContactEdge<Real>* );
  void Delete_Edge( ContactEdge<Real>* );
#ifndef CONTACT_NO_MPI
  void Insert_Edge( char* );
#endif

 private:
  
  ContactEdgeBlock(ContactEdgeBlock&);
  ContactEdgeBlock& operator=(ContactEdgeBlock&);
  
  ContactTopology* topology;

  int number_of_edges;
  int num_edges_added;
  ContactSearch::ContactEdge_Type type;
  int entity_key;
  int id;

  ContactBlockEntityList* edge_list;

};

#endif //_ContactEdgeBlock_h_
