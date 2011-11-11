// $Id$

#ifndef ContactLineEdgeL2_h_
#define ContactLineEdgeL2_h_

#include "ContactEdge.h"

class ContactFixedSizeAllocator;
template<typename DataType> class ContactNode;

template<typename DataType>
class ContactLineEdgeL2 : public ContactEdge<DataType> {

 public:
  ContactLineEdgeL2( int block_index = -1, int index_in_block = -1 );
  static ContactLineEdgeL2<DataType>* new_ContactLineEdgeL2( ContactFixedSizeAllocator&,
						   int block_index=-1, 
						   int index_in_block = -1 );
  ~ContactLineEdgeL2();

 protected:
 private:
  ContactNode<DataType>* nodes[2];
};

#endif // ContactLineEdgeL2_h_
