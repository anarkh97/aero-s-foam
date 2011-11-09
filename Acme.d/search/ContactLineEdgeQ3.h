// $Id$

#ifndef ContactLineEdgeQ3_h_
#define ContactLineEdgeQ3_h_

#include "ContactEdge.h"

class ContactFixedSizeAllocator;
class ContactNode;

class ContactLineEdgeQ3 : public ContactEdge {

 public:
  ContactLineEdgeQ3( int block_index=-1, int host_index_in_block=-1 );
  static ContactLineEdgeQ3* new_ContactLineEdgeQ3( ContactFixedSizeAllocator&,
						   int block_index=-1, 
						   int host_index_in_block = -1 );
  ~ContactLineEdgeQ3();

 protected:
 private:
  ContactNode* nodes[3];
};

#endif // ContactLineEdgeQ3_h_
