// $Id: AcmeLinkList.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeLinkList_h_
#define _AcmeLinkList_h_

#include "AcmeEntity.h"

class AcmeLLnode;

class AcmeLLnode {

 public:
  AcmeLLnode( AcmeEntity* );
  ~AcmeLLnode();
  inline AcmeEntity* Entity() { return entity; };
  inline AcmeLLnode* Next() { return next; };
  inline void Next( AcmeLLnode* llnode ) { next = llnode; };
  
 private:
  AcmeEntity* entity;
  AcmeLLnode* next;
};

class AcmeLinkList {

 public:
  AcmeLinkList();
  ~AcmeLinkList();
  void Clear();
  void Reset();
  void Append( AcmeEntity* );
  void Append( AcmeLLnode* );
  void AppendFirstLink( AcmeLinkList& );
  inline void IteratorStart(AcmeLLnode* link=NULL) {current=(link?link:head);};
  inline AcmeLLnode* Iterator()
    {AcmeLLnode* link=current;
     if (current) current = current->Next();
     return link;};
  inline AcmeLLnode* Head() {return head;};
  inline AcmeLLnode* Tail() {return tail;};
  inline int NumEntities() {return num_entities;};

 private:
  AcmeLLnode* head;
  AcmeLLnode* tail;
  AcmeLLnode* current;
  int num_entities;
};

#endif
