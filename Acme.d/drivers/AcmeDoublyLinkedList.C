// $Id: AcmeDoublyLinkedList.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeEntity.h"
#include "AcmeDoublyLinkedList.h"
#include "contact_assert.h"

#include <iostream.h>
#include <new.h>

AcmeDLL::AcmeDLL( )
{
  head         = NULL;
  tail         = NULL;
  current      = NULL;
  num_entities = 0;
  direction    = 0;
}

AcmeDLL::~AcmeDLL()
{
  Clear();
}

void AcmeDLL::Prepend( AcmeDLLnode* link )
{
  if( head ){
    head->Prev(link);
    link->Next(head);
    head = link;
  } else {
    head = link;
    tail = link;
    link->Prev(NULL);
    link->Next(NULL);
  }
  num_entities++;
}

void AcmeDLL::Append( AcmeEntity* entity )
{
  AcmeDLLnode* link = new AcmeDLLnode( entity );
  if( tail ){
    tail->Next(link);
    link->Prev(tail);
    tail = link;
  } else {
    head = link;
    tail = link;
    link->Prev(NULL);
    link->Next(NULL);
  }
  num_entities++;
}

void AcmeDLL::Append( AcmeDLLnode* link )
{
  if( head ){
    tail->Next(link);
    link->Prev(tail);
    tail = link;
  } else {
    head = link;
    tail = link;
    link->Prev(NULL);
    link->Next(NULL);
  }
  num_entities++;
}

void AcmeDLL::AppendFirstLink( AcmeDLL& list )
{
  AcmeDLLnode* link = head;
  head = head->Next();
  link->Next(NULL);
  list.Append(link);
  num_entities--;
}

void
AcmeDLL::InsertAfter( AcmeDLLnode* after, 
                      AcmeEntity* entity )
{
  if (after==NULL || after==tail) {
    Append(entity);
  } else {
    AcmeDLLnode* link = new AcmeDLLnode( entity );
    AcmeDLLnode* prev = after;
    AcmeDLLnode* next = after->Next();
    prev->Next(link);
    link->Prev(prev);
    next->Prev(link);
    link->Next(next);
    num_entities++;
  }
}

void
AcmeDLL::InsertAfter( AcmeDLLnode* after, 
                      AcmeDLLnode* link )
{
  if (after==NULL || after==tail) {
    Append(link);
  } else {
    AcmeDLLnode* prev = after;
    AcmeDLLnode* next = after->Next();
    prev->Next(link);
    link->Prev(prev);
    next->Prev(link);
    link->Next(next);
    num_entities++;
  }
}

void AcmeDLL::Delete( AcmeDLLnode* link )
{
  if( num_entities==1 ){
    Clear();
    current = head;
  } else {
    if (head==link) {
      head = link->Next();
      head->Prev(NULL);
    }
    if (tail==link) {
      tail = link->Prev();
      tail->Next(NULL);
    }
    if (link->Prev()) link->Prev()->Next(link->Next());
    if (link->Next()) link->Next()->Prev(link->Prev());
    if (current==link) {
      if (direction>0) current = link->Next();
      if (direction<0) current = link->Prev();
    }
    link->~AcmeDLLnode();
    num_entities--;
  }
}

void AcmeDLL::DisplayList()
{
  cout<<"    Displaying Linked List..."<<endl;
  if (num_entities==0) {
    cout<<"      The list is empty"<<endl;
  } else {
    cout<<"      num_entities = "<<num_entities<<endl;
    cout<<"      head = "<<head<<endl;
    cout<<"        entity = "<<head->Entity()<<endl;
    cout<<"        prev   = "<<head->Prev()<<endl;
    cout<<"        next   = "<<head->Next()<<endl;
    cout<<"      tail = "<<tail<<endl;
    cout<<"        entity = "<<tail->Entity()<<endl;
    cout<<"        prev   = "<<tail->Prev()<<endl;
    cout<<"        next   = "<<tail->Next()<<endl;
    AcmeDLLnode* link = head;
    while( link ){
      cout<<"      link = "<<link<<endl;
      cout<<"        entity = "<<link->Entity()<<endl;
      cout<<"        prev   = "<<link->Prev()<<endl;
      cout<<"        next   = "<<link->Next()<<endl;
      link = link->Next();
    }
    cout<<flush;
  }
}

void AcmeDLL::Reset()
{
  head         = NULL;
  tail         = NULL;
  current      = NULL;
  num_entities = 0;
  direction    = 0;
}

void AcmeDLL::Clear()
{
  AcmeDLLnode* link = head;
  while( link ){
    AcmeDLLnode* next = link->Next();
    link->~AcmeDLLnode();
    link = next;
  }
  head         = NULL;
  tail         = NULL;
  current      = NULL;
  num_entities = 0;
  direction    = 0;
}

void
AcmeDLL::Sort() 
{
  if (num_entities<=1) return;
  AcmeDLL *L1 = new AcmeDLL();
  AcmeDLL *L2 = new AcmeDLL();
  IteratorStart();
  while (head) {
    AppendFirstLink(*L1);
    if (head) AppendFirstLink(*L2);
  }
  L1->Sort();
  L2->Sort();
  while (L1->Head() && L2->Head()) {
    if (CompareEntities(L1->Head()->Entity(),L2->Head()->Entity())<0) {
      L1->AppendFirstLink(*this);
    } else {
      L2->AppendFirstLink(*this);
    }
  }
  while (L1->Head()) L1->AppendFirstLink(*this);
  while (L2->Head()) L2->AppendFirstLink(*this);
  delete L1;
  delete L2;
}

int 
AcmeDLL::CompareEntities( AcmeEntity* Entity1, 
			  AcmeEntity* Entity2 )
{
  if (Entity1->BlockID()<Entity2->BlockID()) return -1;
  if (Entity1->BlockID()>Entity2->BlockID()) return  1;
  
  if (Entity1->Exodus_ID()  < Entity2->Exodus_ID()) return -1;
  if (Entity1->Exodus_ID() == Entity2->Exodus_ID()) return  0;
  if (Entity1->Exodus_ID()  > Entity2->Exodus_ID()) return  1;
  // Should never get here but must have a return value
  POSTCONDITION( 0 );
  return -1024;
}

AcmeDLLnode::AcmeDLLnode( AcmeEntity* Entity )
{
  PRECONDITION( Entity );
  entity = Entity;
  prev   = NULL;
  next   = NULL;
}

AcmeDLLnode::~AcmeDLLnode( ) { }
