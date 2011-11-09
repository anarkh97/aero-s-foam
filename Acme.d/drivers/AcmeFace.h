// $Id: AcmeFace.h,v 2002.2 2004/06/23 17:08:55 mwglass Exp $

#ifndef _AcmeFace_h_
#define _AcmeFace_h_

#include "AcmeEntity.h"
#include "AcmeDoublyLinkedList.h"

class AcmeNode;
class AcmeTopologyEntityList;

class AcmeFace : public AcmeEntity {

 public:

  AcmeFace( AcmeEntity::AcmeFaceType, 
            int Block_Index,
	    int Host_Index_in_Block, 
	    int Exo_ID, 
            AcmeNode** NodeList, 
            int* NodeIds );
  virtual ~AcmeFace();

  inline AcmeEntity::AcmeFaceType FaceType() { return type; };

  virtual int Nodes_Per_Face() = 0;
  virtual void Compute_Global_Coordinates( Real*, Real* ) = 0;
  virtual void Compute_Local_Coordinates( Real*, Real* ) = 0;
  virtual void Evaluate_Shape_Functions( Real* local_coord, 
					 Real* shape_fnc) = 0;


  inline AcmeNode* Node( int i ) { return( node_list[i] ); }
  inline AcmeNode** Nodes() { return node_list; };
  inline void NodeIds( int* n ) { for (int i=0; i<number_node_connections; i++) node_ids[i]=n[i]; };
  inline int* NodeIds() { return node_ids; };
  inline int  NodeId(int n) { return node_ids[n]; };

  inline void ConnectNode( int i, AcmeNode* node) { node_list[i] = node; };
  void ConnectNodes(AcmeTopologyEntityList*);

#ifndef CONTACT_NO_MPI
  // Packing/Unpacking Functions
  int  Size();
  void Pack( char* );
  void Unpack( char* );
#endif

 protected:
  int number_node_connections;
  AcmeEntity::AcmeFaceType type;

 private:
  AcmeNode** node_list;
  int*       node_ids;
  
};

#endif // #ifdef _AcmeFace_h_

