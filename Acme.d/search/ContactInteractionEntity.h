// $Id$

#ifndef ContactInteractionEntity_H_
#define ContactInteractionEntity_H_

#include "ContactTopologyEntity.h"
#include "Contact_Defines.h"
#include "ContactParOStream.h"
#include "contact_assert.h"

class ContactInteractionEntity {

 public:
 
  struct entity_data{
    int type;                       //only used by UpdateTopology()
    int owner;
    int block_id;
    int index_in_host_array;        //only used to rtn to host-code
    int index_in_owner_proc_array;
    int host_gid[2];
  };
 
  struct connection_data{
    int type;
    int owner;
    int index_in_owner_proc_array;
    int host_gid[2];
  };
  
  enum Packed_Variables { BASE_TYPE = 0, 
                          NUMBER_PACKED_VARS };

  // Constructors/Destructors
  virtual ~ContactInteractionEntity();
  ContactInteractionEntity(Real *data_array_, ContactType base_type_);

  // Access Functions
  inline void Index(int i)          {index=i;};
  inline int  Index( ) const        {return index;};
  inline void ProcIndex(int i)      {proc_index=i;};
  inline int  ProcIndex( ) const    {return proc_index;};
  inline void EnfArrayIndex(int i)  {enf_index=i;};
  inline int  EnfArrayIndex() const {return enf_index;};  

  void Sort( int, ContactInteractionEntity** );
  
  // Access Functions
  inline ContactType Base_Type()              {return base_type;};
  inline Real* DataArray_Buffer()             {return data_array;};
  inline Real* Variable( VariableHandle & vh) {return DataArray_Buffer()+vh;};
  
  virtual void Display(ContactParOStream&);
  inline void SetEntityData( entity_data*, ContactTopologyEntity<Real>*);

 protected:
  Real *data_array;
  ContactType base_type;

  int index;
  int proc_index;
  int enf_index;
  
  // Parallel Packing/Unpacking Functions
  // These are protected because the derived class function generally needs
  // to add data and should call this function.  It shouldn't be called
  // directly
  inline int  Size();
  inline void Pack( char* );
  inline void Unpack( char* );
  inline void Copy( ContactInteractionEntity* );
  
  inline int  Size_ForSecondary();
  inline void Pack_ForSecondary( char* );
  
  inline int  PackEntityData(entity_data*, int*);
  inline int  UnPackEntityData(entity_data*, int*);
  inline int  PackEntityData(entity_data*, Real*);
  inline int  UnPackEntityData(entity_data*, Real*);

 private:
  // not defined; all ContactEntities are not copyable or assignable
  ContactInteractionEntity( const ContactInteractionEntity& );
  ContactInteractionEntity& operator=( const ContactInteractionEntity& );
  ContactInteractionEntity( ContactInteractionEntity&);

};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack functions that are to be used for DLB
//--------------------------------------------------------------------
inline int ContactInteractionEntity::Size()
{
  return( NUMBER_PACKED_VARS*sizeof(int) );
}

inline void ContactInteractionEntity::Pack( char* buffer )
{
  int* i_buffer = reinterpret_cast<int*> (buffer);
  i_buffer[BASE_TYPE] = Base_Type();
}

inline void ContactInteractionEntity::Unpack( char* buffer )
{
  REMEMBER( int* i_buffer = reinterpret_cast<int*> (buffer)) ;
  PRECONDITION( i_buffer[BASE_TYPE] == Base_Type() );
}

inline void ContactInteractionEntity::Copy( ContactInteractionEntity* src )
{
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Size/Pack/Unpack/Copy functions that are to be used for 
// transferring entities from the primary to secondary decomposition
//--------------------------------------------------------------------
inline int ContactInteractionEntity::Size_ForSecondary()
{
  return( sizeof(int) );
}

inline void ContactInteractionEntity::Pack_ForSecondary( char* buffer )
{
  int* i_buffer = reinterpret_cast<int*> (buffer);
  i_buffer[BASE_TYPE] = Base_Type();
}

inline int ContactInteractionEntity::PackEntityData(entity_data* data, 
                                                    int* buffer)
{
  int* buf = buffer;
  *buf++   = data->type;
  *buf++   = data->owner;
  *buf++   = data->block_id;
  *buf++   = data->index_in_host_array;
  *buf++   = data->index_in_owner_proc_array;
  *buf++   = data->host_gid[0];
  *buf++   = data->host_gid[1];
  POSTCONDITION(7*sizeof(int)==sizeof(entity_data));
  return 7;
}

inline int ContactInteractionEntity::PackEntityData(entity_data* data, 
                                                    Real* buffer)
{
  Real* buf = buffer;
  *buf++    = data->type;
  *buf++    = data->owner;
  *buf++    = data->block_id;
  *buf++    = data->index_in_host_array;
  *buf++    = data->index_in_owner_proc_array;
  *buf++    = data->host_gid[0];
  *buf++    = data->host_gid[1];
  POSTCONDITION(7*sizeof(int)==sizeof(entity_data));
  return 7;
}

inline int ContactInteractionEntity::UnPackEntityData(entity_data* data, 
                                                      int* buffer)
{
  int* buf                        = buffer;
  data->type                      = *buf++;
  data->owner                     = *buf++;
  data->block_id                  = *buf++;
  data->index_in_host_array       = *buf++;
  data->index_in_owner_proc_array = *buf++;
  data->host_gid[0]               = *buf++;
  data->host_gid[1]               = *buf++;
  POSTCONDITION(7*sizeof(int)==sizeof(entity_data));
  return 7;
}

inline int ContactInteractionEntity::UnPackEntityData(entity_data* data, 
                                                      Real* buffer)
{
  Real* buf                       = buffer;
  data->type                      = (int)*buf++;
  data->owner                     = (int)*buf++;
  data->block_id                  = (int)*buf++;
  data->index_in_host_array       = (int)*buf++;
  data->index_in_owner_proc_array = (int)*buf++;
  data->host_gid[0]               = (int)*buf++;
  data->host_gid[1]               = (int)*buf++;
  POSTCONDITION(7*sizeof(int)==sizeof(entity_data));
  return 7;
}

inline void ContactInteractionEntity::SetEntityData( entity_data* data, 
                                                     ContactTopologyEntity<Real>* entity)
{
  data->type                      = entity->Base_Type(); 
  data->owner                     = entity->Owner(); 
  data->block_id                  = entity->BlockID(); 
  data->index_in_host_array       = entity->HostArrayIndex();   
  data->index_in_owner_proc_array = entity->OwnerProcArrayIndex(); 
  data->host_gid[0]               = entity->Global_ID().HiInt(); 
  data->host_gid[1]               = entity->Global_ID().LoInt(); 
}

#endif  // ifdef ContactInteractionEntity_H_
