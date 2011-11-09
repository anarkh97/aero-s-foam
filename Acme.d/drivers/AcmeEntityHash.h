// $Id: AcmeEntityHash.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeEntityHash_
#define _AcmeEntityHash_

#include "AcmeEntity.h"
#include "AcmeLinkList.h"
#include <iostream.h>

class AcmeEntityHash {
  
 public:
  AcmeEntityHash( );
  AcmeEntityHash( int, AcmeEntity** );
  AcmeEntityHash( AcmeLinkList* );
  ~AcmeEntityHash();

  struct hash {
    int global_id;
    AcmeEntity* entity;
    struct hash *next;
  };

  struct hash_list {
    struct hash* ptr;
    struct hash_list *next;
  };

  void SetupHash( int, AcmeEntity** );
  void ReHash( int, AcmeEntity** );
  void ClearHash();
  AcmeEntity* find( int, int, hash*, AcmeEntity* );

 private:
  int nbins;
  int number_of_entities;

  hash* hash_space;
  hash** bins;
  hash_list* hash_space_list;

  void create_space();
};

#endif // _AcmeEntityHash_
