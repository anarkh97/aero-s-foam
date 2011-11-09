// $Id$

#ifndef ContactTopologyEntityList_
#define ContactTopologyEntityList_

#include "ContactTopologyEntity.h"
#include "ContactInteractionEntity.h"
#include "ContactTopologyEntityHash.h"
#include "ContactHostGlobalID.h"
#include <iostream>

class ContactNodeBlock;
class ContactEdgeBlock;
class ContactFaceBlock;
class ContactElementBlock;

class ContactTopologyEntityList : public ContactTopologyEntityHash {
  
 public:
  ContactTopologyEntityList( );
  ~ContactTopologyEntityList( );

  void CleanUp();
  void Rehash();
  void BuildList( ContactNodeBlock** , int, int );
  void BuildList( ContactEdgeBlock** , int, int );
  void BuildList( ContactFaceBlock** , int, int );
  void BuildList( ContactElementBlock** , int, int );
  void BuildList( ContactNodeBlock** , ContactNodeBlock** , int, int );
  void BuildList( ContactEdgeBlock** , ContactEdgeBlock** , int, int );
  void BuildList( ContactFaceBlock** , ContactFaceBlock** , int, int );
  void BuildList( ContactElementBlock** , ContactElementBlock** , int, int );
  
  void SortByNodeGID();
  
  inline ContactTopologyEntity* Find( ContactTopologyEntity::connection_data* data )
    { ContactHostGlobalID GID( data->host_gid[0], data->host_gid[1] );
      return find(GID); };
  inline ContactTopologyEntity* Find( ContactInteractionEntity::entity_data* data )
    { ContactHostGlobalID GID( data->host_gid[0], data->host_gid[1] );
      return find(GID); };
  inline ContactTopologyEntity* Find( ContactHostGlobalID& GID ) 
    {return find(GID);};
  inline ContactTopologyEntity* Find( int index) {return entity_list[index];};
  inline void IteratorStart() {current=0;};
  inline void IteratorEnd  () {current=num_entities-1;};
  inline ContactTopologyEntity* IteratorForward() 
    { if (current==num_entities) return NULL; 
      else return entity_list[current++]; };
  inline ContactTopologyEntity* IteratorBackward()
    { if (current<0) return NULL; 
      else return entity_list[current--]; };
  inline int NumEntities() {return num_entities;};
  inline ContactTopologyEntity** EntityList() {return entity_list;};
  inline ContactTopologyEntity** BlockEntityList(int i) 
    {return num_entities>0?&entity_list[block_offset[i]]:NULL;};
  inline int Num_Blocks() {return nblocks;};
  inline int BlockNumEntities(int i) {return num_entities>0?block_cnt[i]:0;};
  
  
  void Display(ContactParOStream&);

 private:
  int  current;
  int  num_entities;
  int  nblocks;
  int* block_cnt;
  int* block_offset;
  ContactTopologyEntity** entity_list;
};

#endif // ContactTopologyEntityList_
