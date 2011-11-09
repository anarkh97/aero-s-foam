// $Id: AcmeTopologyEntityList.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeTopologyEntityList.h"
#include "AcmeEntityBlock.h"
#include "contact_assert.h"

#define BIN_FRACTION 0.25
#define BIN_MINIMUM 100

AcmeTopologyEntityList::AcmeTopologyEntityList( )
  : AcmeEntityHash( )
{ 
  current      = 0;
  num_entities = 0;
  nblocks      = 0;
  block_cnt    = NULL;
  block_offset = NULL;
  entity_list  = NULL;
}

AcmeTopologyEntityList::~AcmeTopologyEntityList()
{
  CleanUp();
}

void
AcmeTopologyEntityList::CleanUp()
{
  ClearHash();
  if (block_cnt)    delete [] block_cnt;
  if (block_offset) delete [] block_offset;
  if (entity_list)  delete [] entity_list;
  block_cnt    = NULL;
  block_offset = NULL;
  entity_list  = NULL;
  num_entities = 0;
  nblocks      = 0;
}

void
AcmeTopologyEntityList::BuildList( AcmeEntityBlock** blocks, 
                                   int num_blocks, int active_only)
{
  if (nblocks>0) CleanUp();
  nblocks      = num_blocks;
  block_cnt    = new int[nblocks];
  block_offset = new int[nblocks];
  
  int i;
  int cnt=0;
  for (i=0; i<num_blocks; i++) {
    int knt=0;
    if (active_only) {
      AcmeBlockEntityList* block_entities = blocks[i]->EntityList();
      block_entities->IteratorStart();
      while (AcmeEntity* entity=block_entities->IteratorForward()) {
        if (entity->State()==AcmeEntity::ACTIVE) knt++;
      }
    } else {
      knt = blocks[i]->Number_of_Entities();
    }
    block_cnt[i]    = knt;
    block_offset[i] = cnt;
    cnt            += knt;
  }
  entity_list = new AcmeEntity*[cnt];
  
  num_entities = 0;
  for (i=0; i<num_blocks; i++) {
    AcmeBlockEntityList* block_entities = blocks[i]->EntityList();
    block_entities->IteratorStart();
    while (AcmeEntity* entity=block_entities->IteratorForward()) {
      bool add_entity = true;
      if (active_only) {
        if (entity->State()==AcmeEntity::INACTIVE) add_entity = false;
      }
      if (add_entity) {
        entity->ProcArrayIndex(num_entities);
        entity_list[num_entities++] = entity;
      }
    }
  }
  SetupHash(num_entities, entity_list);
}

void
AcmeTopologyEntityList::BuildSortedList( AcmeEntityBlock** blocks, 
                                         int num_blocks, int active_only)
{
  if (nblocks>0) CleanUp();
  nblocks      = num_blocks;
  block_cnt    = new int[nblocks];
  block_offset = new int[nblocks];
  
  int i;
  int cnt=0;
  for (i=0; i<num_blocks; i++) {
    int knt=0;
    if (active_only) {
      AcmeBlockEntityList* block_entities = blocks[i]->EntityList();
      block_entities->IteratorStart();
      while (AcmeEntity* entity=block_entities->IteratorForward()) {
        if (entity->State()==AcmeEntity::ACTIVE) knt++;
      }
    } else {
      knt = blocks[i]->Number_of_Entities();
    }
    block_cnt[i]    = knt;
    block_offset[i] = cnt;
    cnt            += knt;
  }
  entity_list = new AcmeEntity*[cnt];
  
  num_entities = 0;
  for (i=0; i<num_blocks; i++) {
    AcmeEntity** list = &entity_list[num_entities];
    AcmeBlockEntityList* block_entities = blocks[i]->EntityList();
    block_entities->IteratorStart();
    while (AcmeEntity* entity=block_entities->IteratorForward()) {
      bool add_entity = true;
      if (active_only) {
        if (entity->State()==AcmeEntity::INACTIVE) add_entity = false;
      }
      if (add_entity) {
        entity_list[num_entities] = entity;
        num_entities++;
      }
    }
    SortEntityList(cnt, list);
    for (int j=0; j<block_cnt[i]; j++) {
      list[j]->ProcArrayIndex(block_offset[i]+j);
    }
  }
  SetupHash(num_entities, entity_list);
}

void 
AcmeTopologyEntityList::Rehash( )
{
  ReHash(num_entities, entity_list);
}

void AcmeTopologyEntityList::SortEntityList(int cnt, 
                                            AcmeEntity** list)
{
  int i, j, k, n;
  AcmeEntity* entity;

  if (cnt>1) {
    k = (cnt>>1)+1;
    n = cnt;
    for (;;) {
      if (k>1) {
        entity = list[--k-1];
      } else {
        entity    = list[n-1];
        list[n-1] = list[0];
        if (--n == 1) {
          list[0] = entity;
          break;
        }
      }
      i = k;
      j = k<<1;
      while (j<=n) {
        if ((j<n) && (CompareEntities(list[j-1],list[j])==-1)) ++j;
        if (CompareEntities(entity,list[j-1])==-1) {
          list[i-1] = list[j-1];
          j += (i=j);
        } else {
          j = n+1;
        }
      }
      list[i-1] = entity;
    }
  }
}

int AcmeTopologyEntityList::CompareEntities( AcmeEntity* Entity1, 
				             AcmeEntity* Entity2 )
{
  if (Entity1->BlockID()  < Entity2->BlockID())   return -1;
  if (Entity1->BlockID()  > Entity2->BlockID())   return  1;
  if (Entity1->Exodus_ID()< Entity2->Exodus_ID()) return -1;
  if (Entity1->Exodus_ID()==Entity2->Exodus_ID()) return  0;
  if (Entity1->Exodus_ID()> Entity2->Exodus_ID()) return  1;

  // Should never get here but must have a return value
  POSTCONDITION( 0 );
  return -1024;
}
