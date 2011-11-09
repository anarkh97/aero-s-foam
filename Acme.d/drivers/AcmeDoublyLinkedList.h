// $Id: AcmeDoublyLinkedList.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeDoublyLinkedList_h_
#define _AcmeDoublyLinkedList_h_

#include "AcmeEntity.h"

class AcmeDLLnode {

 public:
  AcmeDLLnode( AcmeEntity* );
  ~AcmeDLLnode();
  inline void Entity(AcmeEntity* e) { entity=e; };
  inline AcmeEntity*  Entity() { return entity; };
  inline AcmeDLLnode* Prev()   { return prev; };
  inline AcmeDLLnode* Next()   { return next; };
  inline void Prev( AcmeDLLnode* link ) { prev = link; };
  inline void Next( AcmeDLLnode* link ) { next = link; };
  
 private:
  AcmeEntity*  entity;
  AcmeDLLnode* prev;
  AcmeDLLnode* next;
  AcmeDLLnode* hash_prev;
  AcmeDLLnode* hash_next;
};

class AcmeDLL {

 public:
  AcmeDLL();
  ~AcmeDLL();
  void Append( AcmeEntity* );
  void Append( AcmeDLLnode* );
  void Prepend( AcmeDLLnode* );
  void AppendFirstLink( AcmeDLL& );
  void InsertAfter(AcmeDLLnode*, AcmeEntity*);
  void InsertAfter(AcmeDLLnode*, AcmeDLLnode*);
  
  void Delete( AcmeDLLnode* );
  
  inline void IteratorStart() {current=head; direction=0;};
  inline void IteratorEnd  () {current=tail; direction=0;};
  inline AcmeEntity* IteratorForward()
    { 
      direction = 1;
      AcmeDLLnode* link=current;
      if (current) {
        current = current->Next();
      }
      if (link) return (link->Entity());
      return NULL;
    };
  inline AcmeEntity* IteratorBackward()
    { 
      direction = -1;
      AcmeDLLnode* link=current;
      if (current) {
        current = current->Prev();
      }
      if (link) return (link->Entity());
      return NULL;
    };
  inline AcmeDLLnode* LinkIteratorForward()
    { 
      direction = 1;
      AcmeDLLnode* link=current;
      if (current) {
        current = current->Next();
      }
      return link;
    };
  inline AcmeDLLnode* LinkIteratorBackward()
    { 
      direction = -1;
      AcmeDLLnode* link=current;
      if (current) {
        current = current->Prev();
      }
      return link;
    };

  inline AcmeDLLnode* Head() {return head;};
  inline AcmeDLLnode* Tail() {return tail;};
  inline AcmeDLLnode* Current() {return current;};
  inline int NumEntities() {return num_entities;};
  
  void DisplayList();
  void Reset();
  void Clear();
  void Sort();
  int  CompareEntities( AcmeEntity*, AcmeEntity* );

 protected:
 
 private:
  AcmeDLLnode* head;
  AcmeDLLnode* tail;
  AcmeDLLnode* current;
  int num_entities;
  int direction;
};

#endif
