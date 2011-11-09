// $Id: AcmeElem.h,v 2002.2 2004/06/23 17:08:55 mwglass Exp $

// We currently have two kinds of "elements" in ACME right now. 
// The first is the "extrusion" of a face to an element (i.e., a T3 to Wedge6
//     or a Q4 to Hex8) used for face-face searches
// The second is an true element that only has node connections (nodes do not
//     have back pointers to the elements right now).  These are used for the
//     volume overlap searches.

#ifndef _AcmeElem_h_
#define _AcmeElem_h_

#include "AcmeEntity.h"
#include "AcmeDoublyLinkedList.h"

class AcmeNode;
class AcmeTopologyEntityList;

class AcmeElem : public AcmeEntity {

 public:

  AcmeElem( AcmeEntity::AcmeElemType, 
            int Block_Index,
	    int Host_Index_in_Block, 
	    int Exo_ID, 
            AcmeNode** NodeList, 
            int* NodeIds );
            
  virtual ~AcmeElem();
  
  inline AcmeEntity::AcmeElemType ElemType() { return type; };
  
  virtual int Nodes_Per_Element() = 0;

  virtual void Evaluate_Shape_Functions( Real*, Real* ) = 0;
  virtual void Compute_Global_Coordinates( Real*, Real* ) = 0;
  virtual void Compute_Local_Coordinates( Real*, Real* ) = 0;

  inline AcmeNode** Nodes() { return node_list; };
  inline AcmeNode* Node( int i ) { return( node_list[i] ); }
  inline void NodeIds( int* n ) { for (int i=0; i<number_node_connections; i++) node_ids[i]=n[i]; };
  inline int* NodeIds() { return node_ids; };
  inline int  NodeId(int n) { return node_ids[n]; };
  inline void ConnectNode( int i, AcmeNode* node) { node_list[i] = node; };
  void ConnectNodes(AcmeTopologyEntityList*);

  // Packing/Unpacking Functions
  int  Size();
  void Pack( char* );
  void Unpack( char* );


 protected:
  int number_node_connections;
  AcmeEntity::AcmeElemType type;

 private:
  AcmeNode** node_list;
  int*       node_ids;
};

#endif // #ifdef _AcmeElem_h_
