// $Id$

#include "ContactEntity.h"
#include "ContactEdge.h"
#include "ContactNode.h"
#include "ContactDoublyLinkedList.h"
#include <new>

using namespace std;

ContactTopologyDLL::ContactTopologyDLL() {
  current = data.begin();
  num_entities = 0;
}

ContactTopologyDLL::~ContactTopologyDLL() {
}

int ContactTopologyDLL::Append( ContactTopologyEntity* entity ) {
  data.push_back(entity);
  num_entities++;
  return data.size() - 1;
}

void ContactTopologyDLL::Clear() {
  data.clear();
  current = data.begin();
  num_entities = 0;
}

void ContactTopologyDLL::Display(ContactParOStream& postream)
{
  IteratorStart();
  while( ContactTopologyEntity* entity=IteratorForward() ){
    entity->Display(postream);
  }
}

ContactInteractionDLL::ContactInteractionDLL() {
  current = data.begin();
}

ContactInteractionDLL::~ContactInteractionDLL()
{
  Clear();
}

void ContactInteractionDLL::Append( ContactInteractionEntity* entity )
{
  data.push_back(entity);
}

void ContactInteractionDLL::Clear()
{
  data.clear();
  current = data.begin();
}

void ContactInteractionDLL::Display(ContactParOStream& postream)
{
  int cur_size = data.size();
  for(int i = 0; i < cur_size; ++i) {
    data[i]->Display(postream);
  }
}


    
