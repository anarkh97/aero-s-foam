// $Id$

#ifndef ContactShellNode_h_
#define ContactShellNode_h_

#include "ContactNode.h"
#include <cstring>

class ContactShellNode : public ContactNode<Real> {

 public:

  ContactShellNode( ContactFixedSizeAllocator*,
                    ContactSearch::ContactNode_Type, 
                    int Block_Index=-1,
		    int Index_in_Block=-1 );
  static ContactShellNode* new_ContactShellNode( ContactFixedSizeAllocator*,
					    ContactSearch::ContactNode_Type, 
					    int Block_Index=-1,
					    int Index_in_Block=-1 );
  ~ContactShellNode();

  void  Shell_Node_Base_ID( int id ) { shell_node_base_id = id; };
  int   Shell_Node_Base_ID() { return shell_node_base_id; };
  Real* Previous_Lofting() { return &previous_lofting[0]; };

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Size/Pack/Unpack functions that are to be used for DLB
  //--------------------------------------------------------------------
 virtual int Size() {
    return Size(-2);
  };

  virtual int Size(int state) {
    int s = ContactNode<Real>::Size(state);
    // want interactions
    s += sizeof(int) + 3*sizeof(Real);
    return s;
  };

  virtual void Pack( char* buffer) {
    Pack(buffer, -2);
  };

  virtual void Pack( char* buffer, int state ){
    ContactNode<Real>::Pack(buffer, state);
    std::memcpy( buffer+ContactNode<Real>::Size(state), 
            &shell_node_base_id, sizeof(int) );
    std::memcpy( buffer+ContactNode<Real>::Size(state)+sizeof(int), 
            &previous_lofting[0], 3*sizeof(Real) );
  };
  
  virtual void Unpack( char* buffer ){
    int state = (reinterpret_cast<int*>
                 (buffer+ContactTopologyEntity<Real>::Size(DataArray_Length())+5*sizeof(int)))[0];
    ContactNode<Real>::Unpack(buffer);
    std::memcpy( &shell_node_base_id, 
                 buffer+ContactNode<Real>::Size(state), 
                 sizeof(int) );
    std::memcpy( &previous_lofting[0], 
                 buffer+ContactNode<Real>::Size(state)+sizeof(int),
	         3*sizeof(Real) );
  };
  
  virtual void Copy( ContactNode<Real>* node ) { Copy(node, -2); }
  
  virtual void Copy( ContactNode<Real>* src, int state=-1 ){
    ContactNode<Real>::Copy( src, state );
    ContactShellNode* s = static_cast<ContactShellNode*>(src);
    shell_node_base_id  = s->shell_node_base_id;
    previous_lofting[0] = s->previous_lofting[0];
    previous_lofting[1] = s->previous_lofting[1];
    previous_lofting[2] = s->previous_lofting[2];
  };

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Size/Pack/Unpack/Copy functions that are to be used for 
  // transferring entities from the primary to secondary decomposition
  //--------------------------------------------------------------------
   virtual int Size_ForSecondary() {
    return Size(-2);
  }

  virtual int Size_ForSecondary(int state) {
    int s = ContactNode<Real>::Size_ForSecondary(state);
    // want interactions
    s += sizeof(int) + 3*sizeof(Real);
    return s;
  }

  virtual void Pack_ForSecondary( char* buffer) {
    Pack_ForSecondary(buffer, -2);
  }

  virtual void Pack_ForSecondary( char* buffer, int state ){
    ContactNode<Real>::Pack_ForSecondary(buffer, state);
    std::memcpy( buffer+ContactNode<Real>::Size_ForSecondary(state), 
                 &shell_node_base_id, sizeof(int) );
    std::memcpy( buffer+ContactNode<Real>::Size_ForSecondary(state)+sizeof(int), 
                 &previous_lofting[0], 3*sizeof(Real) );
  };
  
  virtual void Unpack_ForSecondary( char* buffer ){
    int* i_buf = reinterpret_cast<int*> (buffer);
    int state  = ((i_buf[1]>>12) & 0xf)-2;
    ContactNode<Real>::Unpack_ForSecondary(buffer);
    std::memcpy( &shell_node_base_id, 
                 buffer+ContactNode<Real>::Size_ForSecondary(state), 
                 sizeof(int) );
    std::memcpy( &previous_lofting[0], 
                 buffer+ContactNode<Real>::Size_ForSecondary(state)+sizeof(int),
	         3*sizeof(Real) );
  };
  
  virtual void Copy_ForSecondary( ContactNode<Real>* node) 
           { Copy_ForSecondary(node, -2); }
  
  virtual void Copy_ForSecondary( ContactNode<Real>* src, int state=-1 ){
    ContactNode<Real>::Copy_ForSecondary( src, state );
    ContactShellNode* s = static_cast<ContactShellNode*>(src);
    shell_node_base_id  = s->shell_node_base_id;
    previous_lofting[0] = s->previous_lofting[0];
    previous_lofting[1] = s->previous_lofting[1];
    previous_lofting[2] = s->previous_lofting[2];
  };

  bool ConnectedToAllShellFaces();
  
 private:

  Real previous_lofting[3];

};

#endif
