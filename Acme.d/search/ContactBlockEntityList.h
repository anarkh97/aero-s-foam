// $Id$

#ifndef ContactBlockEntityList_
#define ContactBlockEntityList_

#include "ContactTopologyEntity.h"
#include "ContactDoublyLinkedList.h"
#include "ContactParOStream.h"
#include "ContactFixedSizeAllocator.h"
#include <iostream>

class ContactBlockEntityList : public ContactTopologyDLL {
  
 public:
  ContactBlockEntityList( ContactFixedSizeAllocator*, int );
  ContactBlockEntityList( ContactFixedSizeAllocator* );
  ~ContactBlockEntityList();

  struct hash {
    ContactHostGlobalID id;
    ContactTopologyEntity<Real>* entity;
    int list_index;
    struct hash *prev;
    struct hash *next;
  };

  struct hash_bin {
    struct hash* head;
    struct hash* tail;
  };
  
  void CleanUp();
  void Insert( char* );
  void Insert_ForSecondary( char* );
  void Insert( ContactTopologyEntity<Real>* );
  void Delete( ContactTopologyEntity<Real>* );

  void SetupHash( int );
  void Rehash();
  //void UseHash(int s) {do_hash=s;};
  void UseHash(int s) {SetupHash(s);};
  
  ContactTopologyEntity<Real>*  Find( ContactHostGlobalID& );
  ContactTopologyEntity<Real>*  Find( ContactTopologyEntity<Real>::connection_data* );
  ContactTopologyEntity<Real>*  Find( ContactInteractionEntity::entity_data* );
#ifndef CONTACT_NO_MPI
  ContactTopologyEntity<Real>*  CreateEntity(char*);
  ContactTopologyEntity<Real>*  CreateEntity_ForSecondary(char*);
#endif
  
  friend std::ostream& operator<<( std::ostream& os, const ContactBlockEntityList& hash );
  friend ContactParOStream& operator<<( ContactParOStream&  os, const ContactBlockEntityList& hash );

 private:
  int       do_hash;
  int       hash_size;
  int       hash_size_orig;
  int       hash_nbins;
  int       hash_nbins_orig;
  hash_bin* hash_bins;
#ifdef CONTACT_ANALYZE_HASH
  int       hash_collisions;
#endif
  
  ContactFixedSizeAllocator* allocators;
  
  int  hash_func(int);
  void ComputeNbins(int);
  void PrependHash( hash_bin* bin, hash* hlink );
  void AppendHash( hash_bin* bin, hash* hlink );
  void InsertBeforeHash( hash_bin* bin, hash* hlink, hash* ptr );
  
};

inline int
ContactBlockEntityList::hash_func(int lo_int)
{
#if defined(JENKINS96_HASH_FUNC)
  unsigned int mask = hash_nbins-1;
  unsigned int c    = (unsigned int)lo_int;
  unsigned int a    = 0x9e3779b9;
  unsigned int b    = 0x7b6552e3;
  a=a-b;  a=a-c;  a=a^(c>>13);
  b=b-c;  b=b-a;  b=b^(a<<8);
  c=c-a;  c=c-b;  c=c^(b>>13);
  a=a-b;  a=a-c;  a=a^(c>>12);
  b=b-c;  b=b-a;  b=b^(a<<16);
  c=c-a;  c=c-b;  c=c^(b>>5);
  a=a-b;  a=a-c;  a=a^(c>>3);
  b=b-c;  b=b-a;  b=b^(a<<10);
  c=c-a;  c=c-b;  c=c^(b>>15);
  c=(c & mask);
  return (int)c;
#elif defined(JENKINS32_HASH_FUNC)
  unsigned int mask = hash_nbins-1;
  unsigned int key  = (unsigned int)lo_int;
  key += (key << 12);
  key ^= (key >> 22);
  key += (key <<  4);
  key ^= (key >>  9);
  key += (key << 10);
  key ^= (key >>  2);
  key += (key <<  7);
  key ^= (key >> 12);
  key  = (key & mask);
  return (int)key;
#elif defined(WANG32_HASH_FUNC)
  unsigned int mask = hash_nbins-1;
  unsigned int key  = (unsigned int)lo_int;
  key += ~(key << 15);
  key ^=  (key >> 10);
  key +=  (key <<  3);
  key ^=  (key >>  6);
  key += ~(key << 11);
  key ^=  (key >> 16);
  key  =  (key & mask);
  return (int)key;
#elif defined(SEDGEWICK_HASH_FUNC)
  return ((16161*lo_int)%hash_nbins);
#elif defined(KNUTH_HASH_FUNC)
  return ((2654435761U*(unsigned)lo_int)%hash_nbins);
#else
  return (lo_int % hash_nbins);
#endif
}

#endif // ContactBlockEntityList_
