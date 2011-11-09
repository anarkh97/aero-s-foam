// $Id: AcmeBlockEntityList.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeBlockEntityList_
#define _AcmeBlockEntityList_

#include "AcmeEntity.h"
#include "AcmeDoublyLinkedList.h"
#include "AcmeLinkList.h"
#include <iostream.h>

class AcmeBlockEntityList : public AcmeDLL {
  
 public:
  AcmeBlockEntityList( int );
  AcmeBlockEntityList( );
  ~AcmeBlockEntityList( );

  struct hash {
    int id;
    AcmeEntity* entity;
    AcmeDLLnode* link;
    struct hash *prev;
    struct hash *next;
  };

  struct hash_bin {
    struct hash* head;
    struct hash* tail;
  };
  
  void CleanUp();
  void Insert( AcmeEntity* );
  void Delete( AcmeEntity* );

  void SetupHash( int );
  void Rehash();
  
  AcmeEntity*  Find( AcmeEntity* );
  AcmeEntity*  Find( int );
#ifndef CONTACT_NO_MPI
  void Insert( char* );
  AcmeEntity*  CreateEntity(char*);
  AcmeDLLnode* CreateLink(char*);
#endif
  AcmeDLLnode* FindLink( AcmeEntity* );
  AcmeDLLnode* FindLink( int );
  
  void DisplayHash();

 private:
  int       hash_size;
  int       hash_nbins;
  hash_bin* hash_bins;
  
  void PrependHash( hash_bin* bin, hash* hlink );
  void AppendHash( hash_bin* bin, hash* hlink );
  void InsertBeforeHash( hash_bin* bin, hash* hlink, hash* ptr );
  
};

#endif // _AcmeBlockEntityList_
