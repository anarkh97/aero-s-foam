// $Id: AcmeBlockEntityList.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeBlockEntityList.h"
#include "AcmeNode.h"
#include "AcmeFace.h"
#include "AcmeElem.h"
#include "AcmeQuadFaceL4.h"
#include "AcmeQuadFaceQ8.h"
#include "AcmeTriFaceL3.h"
#include "AcmeTriFaceQ6.h"
#ifdef ACME_2D
#include "AcmeLineFaceL2.h"
#include "AcmeLineFaceL3.h"
#endif
#include "AcmeHexElemL8.h"

#include "contact_assert.h"

#define BIN_FRACTION 0.25
#define BIN_MINIMUM 100

#define HashID_FindPtr(a,b)     while((a)!=NULL && b>(a)->id   ) (a)=(a)->next

#define HashID_InvalidPtr(a,b)     ((a)==NULL)||((b)<(a)->id)

#define Hash_AtTheEnd(a)           (a)==NULL
#define HashID_InTheMiddle(a,b)    (b)<(a)->id

AcmeBlockEntityList::AcmeBlockEntityList(int n)
{ 
  hash_size  = 0;
  hash_nbins = 0;
  hash_bins  = NULL;
  SetupHash(n);
}

AcmeBlockEntityList::AcmeBlockEntityList()
{ 
  hash_size  = 0;
  hash_nbins = 0;
  hash_bins  = NULL;
}

AcmeBlockEntityList::~AcmeBlockEntityList()
{
  CleanUp();
}

void
AcmeBlockEntityList::CleanUp()
{
  int i;
  if (hash_nbins>0) {
    for (i=0; i<hash_nbins; i++) {
      hash* hlink = hash_bins[i].head;
      while (hlink) {
        hash* ptr = hlink;
        hlink = hlink->next;
        delete ptr;
      }
    }
    delete [] hash_bins;
    hash_nbins = 0;
    hash_size  = 0;
  }
  Clear();
}

void AcmeBlockEntityList::Insert( AcmeEntity* entity )
{
  int       ID   = entity->Exodus_ID();
  int       ibin = ID % hash_nbins;
  hash_bin* bin  = &hash_bins[ibin];
  hash*     ptr  = bin->head;
  
  if (Hash_AtTheEnd(ptr)) {
    //=============================================
    // HASH BIN IS EMPTY SO WE KNOW THIS ENTITY 
    // IS NOT YET IN THE HASH TABLE SO JUST ADD IT
    //=============================================
    AcmeDLLnode* link = new AcmeDLLnode(entity);
    hash* hlink   = new hash;
    hlink->id     = ID;
    hlink->entity = entity;
    hlink->link   = link;
    PrependHash(bin, hlink);
    Append(link);
  } else {
    //==================================================
    // HASH BIN IS NOT EMPTY SO SEE IF IT'S IN THIS BIN
    //==================================================
    HashID_FindPtr(ptr,ID);
    if( Hash_AtTheEnd(ptr) ) {
      //===========================================
      // ENTITY IS NOT IN THE HASH TABLE SO ADD IT
      //===========================================
      AcmeDLLnode* link = new AcmeDLLnode(entity);
      hash* hlink   = new hash;
      hlink->id     = ID;
      hlink->entity = entity;
      hlink->link   = link;
      AppendHash(bin, hlink);
      Append(link);
    } else if (HashID_InTheMiddle(ptr,ID)) {
      //===========================================
      // ENTITY IS NOT IN THE HASH TABLE SO ADD IT
      //===========================================
      AcmeDLLnode* link = new AcmeDLLnode(entity);
      hash* hlink   = new hash;
      hlink->id     = ID;
      hlink->entity = entity;
      hlink->link   = link;
      InsertBeforeHash(bin, hlink, ptr);
      Append(link);
    } else {
      //=========================================
      // ENTITY IS ALREADY IN THE HASH TABLE BUT 
      // CHECK THE OWNERSHIP OF THE TWO ENTITIES
      //=========================================
      PRECONDITION (ID == ptr->id);
      //=================================================
      // DON'T HAVE TO ADD THIS ENTITY SO JUST DELETE IT 
      //=================================================
      switch (entity->EntityType()) {
      case (AcmeEntity::CT_NODE):
        {
        AcmeNode* node = static_cast<AcmeNode*>(entity);
        node->~AcmeNode();
        }
        break;
      case (AcmeEntity::CT_FACE):
        {
        AcmeFace* face = static_cast<AcmeFace*>(entity);
        face->~AcmeFace();
        }
        break;
      case (AcmeEntity::CT_ELEM):
        {
        AcmeElem* elem = static_cast<AcmeElem*>(entity);
        elem->~AcmeElem();
        }
        break;
      default:
        POSTCONDITION(false);
        break;
      }
    }
  }
}

void AcmeBlockEntityList::Delete( AcmeEntity* entity )
{
  int ID = entity->Exodus_ID();
  
  int       ibin = ID % hash_nbins;
  hash_bin* bin  = &hash_bins[ibin];
  hash*     ptr  = bin->head;
  
  HashID_FindPtr(ptr,ID);
  if (HashID_InvalidPtr(ptr,ID)) return;
  
  // delete the linked list link
  AcmeDLL::Delete( ptr->link );
  
  // delete the hash table link
  if (ptr == bin->head && ptr == bin->tail) {
    // only one entry in the hash bin
    bin->head = NULL;
    bin->tail = NULL;
  } else if (ptr == bin->head) {
    // entry is at the beginning of the hash bin
    hash* next = ptr->next;
    next->prev = NULL;
    bin->head  = next;
  } else if (ptr == bin->tail) {
    // entry is at the end of the hash bin
    hash* prev = ptr->prev;
    prev->next = NULL;
    bin->tail  = prev;
  } else {
    // entry is in the middle of the hash bin
    hash* prev = ptr->prev;
    hash* next = ptr->next;
    prev->next = next;
    next->prev = prev;
  }
  delete ptr;
  hash_size--;
  
  // delete the entity itself
  switch( entity->EntityType() ){
  case AcmeEntity::CT_NODE: 
    {
      AcmeNode* node = static_cast<AcmeNode*>(entity);
      node->~AcmeNode();
    }
    break;
  case AcmeEntity::CT_FACE:
    {
      AcmeFace* face = static_cast<AcmeFace*>(entity);
      face->~AcmeFace();
    }
    break;
  case AcmeEntity::CT_ELEM:
    {
      AcmeElem* elem = static_cast<AcmeElem*>(entity);
      elem->~AcmeElem();
    }
    break;
  }
}

void
AcmeBlockEntityList::SetupHash(int n)
{
  if (hash_nbins>0) {
    for (int i=0; i<hash_nbins; i++) {
      hash* hlink = hash_bins[i].head;
      while (hlink) {
        hash* prev = hlink;
        hlink = hlink->next;
        delete prev;
      }
    }
    delete [] hash_bins;
  }
  
  // find the number of bins
  hash_size  = n;
  hash_nbins = (int) (BIN_FRACTION*hash_size);
  hash_nbins = MAX( hash_nbins, BIN_MINIMUM );

  // Round up nbins to make it more prime-like
  if ( !(hash_nbins % 2) ) hash_nbins += 1;
  if ( !(hash_nbins % 3) ) hash_nbins += 2;
  if ( !(hash_nbins % 6) ) hash_nbins += 6;

  // Allocate the bins
  hash_bins = new hash_bin[hash_nbins];

  // Initialize the bins to empty
  for( int i=0 ; i<hash_nbins ; i++ ) {
    hash_bins[i].head = NULL;
    hash_bins[i].tail = NULL;
  }
}

void
AcmeBlockEntityList::Rehash()
{
  SetupHash(NumEntities());
  
  IteratorStart();
  while (AcmeDLLnode* link=LinkIteratorForward()) {
    AcmeEntity* entity = link->Entity();
    int ID = entity->Exodus_ID();

    int       ibin = ID % hash_nbins;
    hash_bin* bin  = &hash_bins[ibin];
    hash*     ptr  = bin->head;
    
    if (Hash_AtTheEnd(ptr)) {
      hash* hlink   = new hash;
      hlink->id     = ID;
      hlink->entity = entity;
      hlink->link   = link;
      PrependHash(bin, hlink);
    } else {
      HashID_FindPtr(ptr,ID);
      if (Hash_AtTheEnd(ptr)) {
        // at the end of the bin, insert after ptr
        hash* hlink   = new hash;
        hlink->id     = ID;
        hlink->entity = entity;
        hlink->link   = link;
        AppendHash(bin, hlink);
      } else if (HashID_InTheMiddle(ptr,ID)) {
        // middle of the bin, insert before ptr
        hash* hlink   = new hash;
        hlink->id     = ID;
        hlink->entity = entity;
        hlink->link   = link;
        InsertBeforeHash(bin, hlink, ptr);
      }
    }
  }
}

void
AcmeBlockEntityList::PrependHash( hash_bin* bin, hash* hlink )
{
  if (bin->head == NULL) {
    // HASH BIN IS EMPTY
    hlink->prev = NULL;
    hlink->next = NULL;
    bin->head   = hlink;
    bin->tail   = hlink;
  } else {
    hash* next  = bin->head;
    next->prev  = hlink;
    hlink->prev = NULL;
    hlink->next = next;
    bin->head   = hlink;
  }
}

void
AcmeBlockEntityList::AppendHash( hash_bin* bin, hash* hlink )
{
  if (bin->tail == NULL) {
    // HASH BIN IS EMPTY
    hlink->prev = NULL;
    hlink->next = NULL;
    bin->head   = hlink;
    bin->tail   = hlink;
  } else {
    hash* prev  = bin->tail;
    prev->next  = hlink;
    hlink->prev = prev;
    hlink->next = NULL;
    bin->tail   = hlink;
  }
}

void
AcmeBlockEntityList::InsertBeforeHash( hash_bin* bin, hash* hlink, hash* ptr )
{
  PRECONDITION (ptr != NULL);
  if (ptr == bin->head) {
    PrependHash(bin, hlink);
  } else {
    hash* prev  = ptr->prev;
    hash* next  = ptr;
    hlink->prev = prev;
    hlink->next = next;
    prev->next  = hlink;
    next->prev  = hlink;
  }
}

AcmeEntity* 
AcmeBlockEntityList::Find( int ID )
{
  if (hash_nbins==0) return(NULL);
  int   ibin = ID % hash_nbins;
  hash* ptr  = hash_bins[ibin].head;
  
  HashID_FindPtr(ptr,ID);
  if( HashID_InvalidPtr(ptr,ID) ) return(NULL);
  POSTCONDITION(ID == ptr->id);
  return ptr->entity;
}

AcmeEntity* 
AcmeBlockEntityList::Find( AcmeEntity* entity  )
{
  if (hash_nbins==0) return(NULL);
  int   ID   = entity->Exodus_ID();
  int   ibin = ID % hash_nbins;
  hash* ptr  = hash_bins[ibin].head;
  
  HashID_FindPtr(ptr,ID);
  if( HashID_InvalidPtr(ptr,ID) ) return(NULL);
  POSTCONDITION(ID == ptr->id);
  return ptr->entity;
}


AcmeDLLnode* 
AcmeBlockEntityList::FindLink( int ID )
{
  if (hash_nbins==0) return(NULL);
  int   ibin = ID % hash_nbins;
  hash* ptr  = hash_bins[ibin].head;
  
  HashID_FindPtr(ptr,ID);
  if( HashID_InvalidPtr(ptr,ID) ) return(NULL);
  POSTCONDITION(ID == ptr->id);
  return ptr->link;
}

AcmeDLLnode* 
AcmeBlockEntityList::FindLink( AcmeEntity* entity )
{
  if (hash_nbins==0) return(NULL);
  int ID = entity->Exodus_ID();
  
  int   ibin = ID % hash_nbins;
  hash* ptr  = hash_bins[ibin].head;
  
  HashID_FindPtr(ptr,ID);
  if( HashID_InvalidPtr(ptr,ID) ) return(NULL);
  POSTCONDITION(ID == ptr->id);
  return ptr->link;
}

#ifndef CONTACT_NO_MPI

void AcmeBlockEntityList::Insert( char* buffer )
{
  int*  ibuffer = REINTERPRET_CAST(int*) (buffer);
  
  int ID = ibuffer[AcmeEntity::EXO_ID];

  int       ibin = ID % hash_nbins;
  hash_bin* bin  = &hash_bins[ibin];
  hash*     ptr  = bin->head;
  
  if (Hash_AtTheEnd(ptr)) {
    //=============================================
    // HASH BIN IS EMPTY SO WE KNOW THIS ENTITY 
    // IS NOT YET IN THE HASH TABLE SO JUST ADD IT
    //=============================================
    AcmeDLLnode* link = CreateLink(buffer);
    hash* hlink   = new hash;
    hlink->id     = ID;
    hlink->entity = link->Entity();
    hlink->link   = link;
    PrependHash(bin, hlink);
    Append(link);
  } else {
    //==================================================
    // HASH BIN IS NOT EMPTY SO SEE IF IT'S IN THIS BIN
    //==================================================
    HashID_FindPtr(ptr,ID);
    if (Hash_AtTheEnd(ptr)) {
      //===========================================
      // ENTITY IS NOT IN THE HASH TABLE SO ADD IT
      //===========================================
      AcmeDLLnode* link = CreateLink(buffer);
      hash* hlink   = new hash;
      hlink->id     = ID;
      hlink->entity = link->Entity();
      hlink->link   = link;
      AppendHash(bin, hlink);
      Append(link);
    } else if (HashID_InTheMiddle(ptr,ID)) {
      //===========================================
      // ENTITY IS NOT IN THE HASH TABLE SO ADD IT
      //===========================================
      AcmeDLLnode* link = CreateLink(buffer);
      hash* hlink   = new hash;
      hlink->id     = ID;
      hlink->entity = link->Entity();
      hlink->link   = link;
      InsertBeforeHash(bin, hlink, ptr);
      Append(link);
    } else {
      //=========================================
      // ENTITY IS ALREADY IN THE HASH TABLE
      //=========================================
      PRECONDITION (ID == ptr->id);
      if (ibuffer[AcmeEntity::OWNERSHIP]==AcmeEntity::OWNED &&
          ptr->link->Entity()->Ownership()!=AcmeEntity::OWNED) {
        //===============================================
        // IF THE NEW ENTITY IS OWNED AND THE OLD ENTITY
        // IS NOT_OWNED THEN REPLACE IT WITH THE NEW ONE
        //===============================================
        // delete old entity and replace it with this one
        AcmeEntity* old_entity = ptr->link->Entity();
        switch (old_entity->EntityType()) {
        case (AcmeEntity::CT_NODE):
          {
          AcmeNode* node = static_cast<AcmeNode*>(old_entity);
          node->~AcmeNode();
          }
          break;
        case (AcmeEntity::CT_FACE):
          {
          AcmeFace* face = static_cast<AcmeFace*>(old_entity);
          face->~AcmeFace();
          }
          break;
        case (AcmeEntity::CT_ELEM):
          {
          AcmeElem* elem = static_cast<AcmeElem*>(old_entity);
          elem->~AcmeElem();
          }
          break;
	default:
	  POSTCONDITION(false);
	  break;
        }
        AcmeEntity* entity = CreateEntity(buffer);
        ptr->link->Entity(entity);
      }
    }
  }
}

AcmeDLLnode* 
AcmeBlockEntityList::CreateLink(char* buffer)
{
  AcmeDLLnode* link = NULL;
  int*  ibuffer     = REINTERPRET_CAST(int*) (buffer);
  int   entity_type = ibuffer[AcmeEntity::BASE_TYPE];
  switch (entity_type) {
  case AcmeEntity::CT_NODE:
    {
      AcmeEntity::AcmeNodeType node_type = (AcmeEntity::AcmeNodeType)(ibuffer[AcmeEntity::DERIVED_TYPE]);
      AcmeNode* node = new AcmeNode(node_type);
      POSTCONDITION(node);
      node->Unpack(buffer);
      link = new AcmeDLLnode(node);
    }
  break;
  case AcmeEntity::CT_FACE:
    {
      AcmeFace* face = NULL;
      int face_type = ibuffer[AcmeEntity::DERIVED_TYPE];
      switch( face_type ){
      case AcmeEntity::QUADFACEL4 :
        face = new AcmeQuadFaceL4();
        break;
      case AcmeEntity::QUADFACEQ8 :
        face = new AcmeQuadFaceQ8();
        break;
      case AcmeEntity::TRIFACEL3 :
        face = new AcmeTriFaceL3();
        break;
      case AcmeEntity::TRIFACEQ6 :
        face = new AcmeTriFaceQ6();
        break;
#ifdef ACME_2D
      case AcmeEntity::LINEFACEL2 :
        face = new AcmeLineFaceL2();
        break;
      case AcmeEntity::LINEFACEL3 :
        face = new AcmeLineFaceL3();
        break;
#endif
      default:
        POSTCONDITION(false);
        break;
      }
      POSTCONDITION(face);
      face->Unpack(buffer);
      link = new AcmeDLLnode(face);
    }
  break;
  case AcmeEntity::CT_ELEM:
    {
      AcmeElem* elem = NULL;
      int elem_type = ibuffer[AcmeEntity::DERIVED_TYPE];
      switch( elem_type ){
      case AcmeEntity::HEXELEML8 :
        elem = new AcmeHexElemL8();
        break;
      default:
        POSTCONDITION(false);
        break;
      }
      POSTCONDITION(elem);
      elem->Unpack(buffer);
      link = new AcmeDLLnode(elem);
    }
  break;
  default:
    POSTCONDITION(false);
    break;
  }
  POSTCONDITION(link);
  return link;
}

AcmeEntity* 
AcmeBlockEntityList::CreateEntity(char* buffer)
{
  AcmeEntity* entity = NULL;
  int*  ibuffer     = REINTERPRET_CAST(int*) (buffer);
  int   entity_type = ibuffer[AcmeEntity::BASE_TYPE];
  switch (entity_type) {
  case AcmeEntity::CT_NODE:
    {
      AcmeEntity::AcmeNodeType node_type = (AcmeEntity::AcmeNodeType)(ibuffer[AcmeEntity::DERIVED_TYPE]);
      AcmeNode* node = new AcmeNode(node_type);
      POSTCONDITION(node);
      node->Unpack(buffer);
      entity = node;
    }
  break;
  case AcmeEntity::CT_FACE:
    {
      AcmeFace* face = NULL;
      int face_type = ibuffer[AcmeEntity::DERIVED_TYPE];
      switch( face_type ){
      case AcmeEntity::QUADFACEL4 :
        face = new AcmeQuadFaceL4();
        break;
      case AcmeEntity::QUADFACEQ8 :
        face = new AcmeQuadFaceQ8();
        break;
      case AcmeEntity::TRIFACEL3 :
        face = new AcmeTriFaceL3();
        break;
      case AcmeEntity::TRIFACEQ6 :
        face = new AcmeTriFaceQ6();
        break;
#ifdef ACME_2D
      case AcmeEntity::LINEFACEL2 :
        face = new AcmeLineFaceL2();
        break;
      case AcmeEntity::LINEFACEL3 :
        face = new AcmeLineFaceL3();
        break;
#endif
      default:
        POSTCONDITION(false);
        break;
      }
      POSTCONDITION(face);
      face->Unpack(buffer);
      entity = face;
    }
  break;
  case AcmeEntity::CT_ELEM:
    {
      AcmeElem* elem = NULL;
      int elem_type = ibuffer[AcmeEntity::DERIVED_TYPE];
      switch( elem_type ){
      case AcmeEntity::HEXELEML8 :
        elem = new AcmeHexElemL8();
        break;
      default:
        POSTCONDITION(false);
        break;
      }
      POSTCONDITION(elem);
      elem->Unpack(buffer);
      entity = elem;
    }
  break;
  default:
    POSTCONDITION(false);
    break;
  }
  POSTCONDITION(entity);
  return entity;
}

void 
AcmeBlockEntityList::DisplayHash()
{
  cout<<"    Displaying Hash Table..."<<endl;
  if (NumEntities()==0) {
    cout<<"      The table is empty"<<endl;
  } else {
    if (hash_nbins>0) {
      for (int i=0; i<hash_nbins; i++) {
        hash* hlink = hash_bins[i].head;
        if (hlink) cout<<"      bin "<<i<<endl;
        while (hlink) {
          cout<<"        entity = "<<hlink->entity<<endl;
          cout<<"        link   = "<<hlink->link<<endl;
          hlink = hlink->next;
        }
      }
    }
  }
  cout<<flush;
}
#endif
