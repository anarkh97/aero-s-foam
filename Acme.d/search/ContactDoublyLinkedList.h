// $Id$

#ifndef ContactDoublyLinkedList_h_
#define ContactDoublyLinkedList_h_

#include <vector>
#include "ContactTopologyEntity.h"
#include "ContactInteractionEntity.h"
#include "contact_assert.h"
#include "Contact_Defines.h"

class ContactTopologyDLL {

 public:
  ContactTopologyDLL();
  ~ContactTopologyDLL();
  int Append( ContactTopologyEntity<Real>* );
  
  inline void IteratorStart() {
    current = data.begin(); 
    while(current < data.end()) {
      if(*current != NULL) return;
      ++current;
    }
  };
  inline ContactTopologyEntity<Real>* IteratorForward() { 
    while(current < data.end()) {
      ContactTopologyEntity<Real> *cur_object = *current;
      current++;
      if(cur_object != NULL) return cur_object;
    }
    return NULL;
  };

  void SetEntity(const int pos, ContactTopologyEntity<Real>* entity) {
    data[pos] = entity;
  }

  void RemoveEntity(const int pos) {
    data[pos] = NULL;
    num_entities--;
  }

  inline int NumEntities() {return num_entities;};
  
  void Display(ContactParOStream&);
  void Clear();
  
  std::vector< ContactTopologyEntity<Real>* > *EntityList() { return &data; };

 private:
  std::vector< ContactTopologyEntity<Real>* > data;
  std::vector< ContactTopologyEntity<Real>* >::iterator current;
  int num_entities;
};

//=================================================================================================================

class ContactInteractionDLL {

 public:
  ContactInteractionDLL();
  ~ContactInteractionDLL();
  void Append( ContactInteractionEntity* );
  
  void DeletePrev() {
    std::vector< ContactInteractionEntity* >::iterator data_to_erase = current; 
    data_to_erase--;   
    data.erase(data_to_erase);
  }
  
  inline void IteratorStart() {current = data.begin();};

  inline ContactInteractionEntity* IteratorForward() { 
    if(current < data.end()) {
      return *(current++);
    } else {
      return NULL;
    }
  };

  inline ContactInteractionEntity* HeadEntity() {return data[0];};

  inline ContactInteractionEntity* CurrentEntity() {return *current;};

  inline int NumEntities() {return data.size();};
  
  void Display(ContactParOStream&);
  void Reset();
  void Clear();

 private:
  std::vector< ContactInteractionEntity* > data;
  std::vector< ContactInteractionEntity* >::iterator current;
};

#endif
