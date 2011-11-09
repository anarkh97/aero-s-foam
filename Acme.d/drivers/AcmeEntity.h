// $Id: AcmeEntity.h,v 2002.2 2004/06/24 14:51:35 mwglass Exp $

#ifndef _AcmeEntity_H_
#define _AcmeEntity_H_

#include "AcmeZoltanID.h"
#include "Contact_Defines.h"

#ifndef CONTACT_NO_MPI
#include "lbi_const.h"
#endif

#include <string.h>

class AcmeEntity {

 public:
 
  enum BaseEntityType { CT_UNKNOWN=0, CT_NODE, 
		        CT_FACE, CT_ELEM, 
		        CT_NUM_ENTITY_TYPES };
  enum DerivedEntityType{ ACME_UNKNOWN=0, 
                          ACME_NODE, ACME_POINT, 
                          ACME_QUADFACEL4, ACME_QUADFACEQ8, 
                          ACME_TRIFACEL3, ACME_TRIFACEQ6,
                          ACME_LINEFACEL2, ACME_LINEFACEQ3, 
                          ACME_HEXELEML8,
                          NDERIVED_TYPES };
  enum AcmeNodeType{ NODE, POINT,
                     NNODE_TYPES };
  enum AcmeFaceType{ QUADFACEL4, QUADFACEQ8, 
                     TRIFACEL3, TRIFACEQ6,
                     LINEFACEL2, LINEFACEQ3, 
                     NFACE_TYPES };
  enum AcmeElemType{ HEXELEML8, 
                     NELEM_TYPES };
                    
  enum EntityState  {INACTIVE=0, ACTIVE};
  
  enum EntityAction {NONE=0, MARK_FOR_BIRTH, MARK_FOR_DEATH, MARK_FOR_DEATH1};
                  
  enum ProcessorOwnership { NOT_OWNED=0, OWNED };
  
  enum Packed_Variables   { BASE_TYPE = 0, DERIVED_TYPE, 
			    EXO_ID, OWNERSHIP, OWNER, IS_ACTIVE, TAG,
			    BLOCK_INDEX, PROC_INDEX, HOST_INDEX, BLOCK_ID,
                            NUM_DATA, NUM_ATTR,
			    NUMBER_PACKED_VARS };

  int temp_tag;    

  // Constructors/Destructors
  AcmeEntity( AcmeEntity::BaseEntityType, 
              int, int, int );
  virtual ~AcmeEntity();

  inline void  NumData(int n) {
    num_data=n;
    if (data!=NULL) delete [] data;
    data = new Real[num_data];
  };
  inline int   NumData() { return num_data; };
  inline void  Data(Real* a) { memcpy(data, a, num_data*sizeof(Real)); };
  inline Real* Data() { return data; };
  inline Real  Data(int i) { return data[i]; };
  
  inline void  NumAttributes(int n) {
    num_attributes=n;
    if (attributes!=NULL) delete [] attributes;
    attributes = new Real[num_attributes];
  };
  inline int   NumAttributes() { return num_attributes; };
  inline void  Attributes(Real* a) { memcpy(attributes, a, num_attributes*sizeof(Real)); };
  inline Real* Attributes() { return attributes; };
  inline Real  Attributes(int i) { return attributes[i]; };
  
  // Access Functions
  inline ProcessorOwnership Ownership() { return ownership; };
  inline void Ownership( ProcessorOwnership  PO) { ownership = PO;};
  inline int  Owner() { return owner; };
  inline void Owner(int proc) { owner = proc; };
  inline int  Exodus_ID() { return exodus_id; };
  inline void Exodus_ID(int id) { exodus_id=id; };
  inline int  BlockID() { return block_id; };
  inline void BlockID(int i) { block_id=i; };
  inline int  ProcArrayIndex() { return proc_array_index; };
  inline void ProcArrayIndex(int i) { proc_array_index = i; };
  inline int  HostArrayIndex() { return host_array_index; };
  inline void HostArrayIndex(int i) { host_array_index = i; };
  inline int  BlockArrayIndex() { return block_array_index; };
  inline void BlockArrayIndex(int i) { block_array_index = i; };

  // Access Functions
  inline AcmeEntity::BaseEntityType EntityType( ) {return entity_type;};
  inline void State(AcmeEntity::EntityState s ) { state=s; };
  inline AcmeEntity::EntityState State( ) {return state;};
  inline void Action(AcmeEntity::EntityAction s ) { action=s; };
  inline AcmeEntity::EntityAction Action( ) {return action;};
  
#ifndef CONTACT_NO_MPI
  inline AcmeZoltanLID& Zoltan_LID() { return zoltan_lid; };
  inline AcmeZoltanGID& Zoltan_GID() { return zoltan_gid; };
  inline void ZoltanLID(LB_ID_PTR lid) 
    {zoltan_lid.ZoltanLID((int)entity_type, proc_array_index, lid);};
  inline void ZoltanGID(LB_ID_PTR gid) 
    {zoltan_gid.ZoltanGID((int)entity_type, exodus_id, gid);};
#endif

 protected:
 
  AcmeEntity::BaseEntityType    entity_type;
  AcmeEntity::EntityState       state;
  AcmeEntity::EntityAction      action;
  
  // Parallel Packing/Unpacking Functions
  // These are protected because the derived class function generally needs
  // to add data and should call this function.  It shouldn't be called
  // directly
  int  Size();
  void Pack( char* );
  void Unpack( char* );

  int block_array_index;        // This is an object's 0..N-1 index into
                                // the host array (on a per block basis)
  int proc_array_index;         // This is an object's 0..N-1 index 
                                // on a processor
  int host_array_index;         // This is an object's 0..N-1 index 
                                // on a processor
  int block_id;                 // This is the index of the block 
                                // this entity belongs to
  int exodus_id;
  int owner;
  
#ifndef CONTACT_NO_MPI
  AcmeZoltanLID zoltan_lid;
  AcmeZoltanGID zoltan_gid;
#endif
  
  ProcessorOwnership ownership;
#ifndef CONTACT_NO_MPI
  //ContactZoltanLID zoltan_lid;
  //ContactZoltanGID zoltan_gid;
#endif
  int   num_data;
  Real* data;
  int   num_attributes;
  Real* attributes;

 private:
  // not defined; all ContactEntities are not copyable or assignable
  AcmeEntity( const AcmeEntity& );
  AcmeEntity& operator=( const AcmeEntity& );

};

#endif  // ifdef _AcmeEntity_H_
