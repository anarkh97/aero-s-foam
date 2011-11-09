// $Id: AcmeEntityHash.C,v 2002.1 2004/06/18 17:47:04 mwglass Exp $

#include "AcmeEntityHash.h"

#define BIN_FRACTION 0.25
#define BIN_MINIMUM 100


AcmeEntityHash::AcmeEntityHash( )
{
  nbins              = 0;
  number_of_entities = 0;
  bins               = NULL;
  hash_space         = NULL;
  hash_space_list    = NULL;
}

AcmeEntityHash::AcmeEntityHash(int number_entities, 
			       AcmeEntity** entities )
{
  number_of_entities = number_entities;
  nbins              = 0;
  bins               = NULL;
  hash_space         = NULL;
  hash_space_list    = NULL;


  // Create space for bins
  nbins = (int) (BIN_FRACTION*number_of_entities);
  nbins = MAX( nbins, BIN_MINIMUM );

  // Round up nbins to make it more prime-like
  if ( !(nbins % 2) ) nbins += 1;
  if ( !(nbins % 3) ) nbins += 2;
  if ( !(nbins % 6) ) nbins += 6;

  // Allocate the bins
  bins = new hash*[nbins];

  // Initialize the bins to empty
  int i;
  for( i=0 ; i<nbins ; i++ )
    bins[i] = NULL;

  // Create hash space
  create_space();

  // Add each entity to the hash table
  for( i=number_entities-1 ; i>=0 ; i-- )
    find( entities[i]->Exodus_ID(), 1, NULL, entities[i] );
}

AcmeEntityHash::AcmeEntityHash(AcmeLinkList* link_list )
{
  number_of_entities = link_list->NumEntities();
  nbins              = 0;
  bins               = NULL;
  hash_space         = NULL;
  hash_space_list    = NULL;


  // Create space for bins
  nbins = (int) (BIN_FRACTION*number_of_entities);
  nbins = MAX( nbins, BIN_MINIMUM );

  // Round up nbins to make it more prime-like
  if ( !(nbins % 2) ) nbins += 1;
  if ( !(nbins % 3) ) nbins += 2;
  if ( !(nbins % 6) ) nbins += 6;

  // Allocate the bins
  bins = new hash*[nbins];

  // Initialize the bins to empty
  int i;
  for( i=0 ; i<nbins ; i++ )
    bins[i] = NULL;

  // Create hash space
  create_space();

  // Add each entity to the hash table
  for(AcmeLLnode* llnode=link_list->Head(); llnode; llnode=llnode->Next())
    find( llnode->Entity()->Exodus_ID(), 1, NULL, llnode->Entity() );
}

AcmeEntityHash::~AcmeEntityHash()
{
  hash_list* next_list;
  while( hash_space_list != NULL ){
    delete [] hash_space_list->ptr;
    next_list = hash_space_list->next;
    delete hash_space_list;
    hash_space_list = next_list;
  }
  delete [] bins;
}

void
AcmeEntityHash::SetupHash(int number_entities, 
			  AcmeEntity** entities )
{
  number_of_entities = number_entities;
  nbins              = 0;
  bins               = NULL;
  hash_space         = NULL;
  hash_space_list    = NULL;

  // Create space for bins
  nbins = (int) (BIN_FRACTION*number_of_entities);
  nbins = MAX( nbins, BIN_MINIMUM );

  // Round up nbins to make it more prime-like
  if ( !(nbins % 2) ) nbins += 1;
  if ( !(nbins % 3) ) nbins += 2;
  if ( !(nbins % 6) ) nbins += 6;

  // Allocate the bins
  bins = new hash*[nbins];

  // Initialize the bins to empty
  int i;
  for( i=0 ; i<nbins ; i++ )
    bins[i] = NULL;

  // Create hash space
  create_space();

  // Add each entity to the hash table
  for( i=number_entities-1 ; i>=0 ; i-- )
    find( entities[i]->Exodus_ID(), 1, NULL, entities[i] );
}

void
AcmeEntityHash::ReHash( int number_entities, 
			AcmeEntity** entities )
{
  ClearHash();
  SetupHash(number_entities, entities);
}

void
AcmeEntityHash::ClearHash()
{
  hash_list* next_list;
  while( hash_space_list != NULL ){
    delete [] hash_space_list->ptr;
    next_list = hash_space_list->next;
    delete hash_space_list;
    hash_space_list = next_list;
  }
  delete [] bins;
  nbins              = 0;
  number_of_entities = 0;
  bins               = NULL;
  hash_space         = NULL;
  hash_space_list    = NULL;
}

void AcmeEntityHash::create_space()
{
  int n = number_of_entities;
  if( n<=0 ) n=1;

  hash_space          = new hash[n];
  hash_list* new_list = new hash_list;
  new_list->next      = hash_space_list;
  new_list->ptr       = hash_space;
  hash_space_list     = new_list;
}


AcmeEntity* AcmeEntityHash::find( int Global_ID, int flag, 
				  hash* newpt, AcmeEntity* Entity )
{
  hash* ptr;
  hash** previous;
  int ibin;

  ibin = Global_ID % nbins;
  previous = bins + ibin;
  ptr = bins[ibin];
  
  while( ptr != NULL && Global_ID < ptr->global_id ){
    previous = &(ptr->next);
    ptr = ptr->next;
  }
  
  if( ptr == NULL || Global_ID > ptr->global_id ){
    if( flag ){
      if( newpt == NULL ) newpt = hash_space++;
      *previous = newpt;
      newpt->global_id = Global_ID;
      newpt->entity = Entity;
      newpt->next = ptr;
    }
    return(NULL);
  }
  else if( flag == 2 ){
    ptr->entity = Entity;
  }
  return ptr->entity;
}

