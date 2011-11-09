// $Id: AcmeTopologyEntityList.h,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#ifndef _AcmeTopologyEntityList_
#define _AcmeTopologyEntityList_

#include "AcmeEntity.h"
#include "AcmeEntityHash.h"
#include <iostream.h>

class AcmeEntityBlock;

class AcmeTopologyEntityList : public AcmeEntityHash {
  
 public:
  AcmeTopologyEntityList( );
  ~AcmeTopologyEntityList( );

  void CleanUp();
  void Rehash();
  void BuildList( AcmeEntityBlock** blocks, int count, int active_only=0);
  void BuildSortedList( AcmeEntityBlock** blocks, int count, int active_only=0);
  
  inline AcmeEntity* Find( int GID ) {return find(GID,0,NULL,0);};
  inline AcmeEntity* Entity( int index) {return entity_list[index];};
  inline void IteratorStart() {current=0;};
  inline void IteratorEnd  () {current=num_entities-1;};
  inline AcmeEntity* IteratorForward() 
    { if (current==num_entities) return NULL; 
      else return entity_list[current++]; };
  inline AcmeEntity* IteratorBackward()
    { if (current<0) return NULL; 
      else return entity_list[current--]; };
  inline int NumEntities() {return num_entities;};
  inline AcmeEntity** EntityList() {return entity_list;};
  inline AcmeEntity** BlockEntityList(int i) 
    {return num_entities>0?&entity_list[block_offset[i]]:NULL;};
  inline int Num_Blocks() {return nblocks;};
  inline int BlockNumEntities(int i) {return num_entities>0?block_cnt[i]:0;};

 private:
  int  current;
  int  num_entities;
  int  nblocks;
  int* block_cnt;
  int* block_offset;
  AcmeEntity** entity_list;
  void SortEntityList(int, AcmeEntity**);
  int  CompareEntities(AcmeEntity*, AcmeEntity*);
  
};

#endif // _AcmeTopologyEntityList_
