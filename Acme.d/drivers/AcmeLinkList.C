// $Id: AcmeLinkList.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeLinkList.h"
#include "AcmeEntity.h"

#include <new.h>

AcmeLinkList::AcmeLinkList()
{
  head         = NULL;
  tail         = NULL;
  current      = NULL;
  num_entities = 0;
}

AcmeLinkList::~AcmeLinkList()
{
  Clear();
}

void AcmeLinkList::Append( AcmeEntity* entity )
{
  AcmeLLnode* llnode = new AcmeLLnode( entity );
  if( head ){
    tail->Next( llnode );
    tail = llnode;
  } else {
    head = llnode;
    tail = llnode;
  }
  num_entities++;
}

void AcmeLinkList::Append( AcmeLLnode* llnode )
{
  if( head ){
    tail->Next( llnode );
    tail = llnode;
  } else {
    head = llnode;
    tail = llnode;
  }
  num_entities++;
}

void AcmeLinkList::AppendFirstLink( AcmeLinkList& list )
{
  AcmeLLnode* link = head;
  head = head->Next();
  link->Next(NULL);
  list.Append(link);
  num_entities--;
}

void AcmeLinkList::Clear()
{
  AcmeLLnode* llnode = head;
  while( llnode ){
    AcmeLLnode* next = llnode->Next();
    delete llnode ;
    llnode = next;
  }
  head         = NULL;
  tail         = NULL;
  current      = NULL;
  num_entities = 0;
}

void AcmeLinkList::Reset()
{
  head         = NULL;
  tail         = NULL;
  current      = NULL;
  num_entities = 0;
}

AcmeLLnode::AcmeLLnode( AcmeEntity* Entity )
{
  entity = Entity;
  next = NULL;
}

AcmeLLnode::~AcmeLLnode() { }
