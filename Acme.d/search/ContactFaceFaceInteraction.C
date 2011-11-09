// $Id$

#include "allocators.h"
#include "ContactFaceFaceInteraction.h"
#include "CString.h"
#include "ContactTopologyEntityList.h"
#include "ContactTopologyEntityHash.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactTopology.h"
#include "ContactFaceBlock.h"
#include <cstddef>
#include <cstring>
#include <new>

ContactFaceFaceInteraction::ContactFaceFaceInteraction( )
   : ContactInteractionEntity(DataArray, CT_FFI)
{
  slave_face  = NULL;
  master_face = NULL;
  num_edges   = 0;
  num_edges   = 0;  // why twice REJ?
  vertices    = NULL;
}                                                                               

ContactFaceFaceInteraction::ContactFaceFaceInteraction( ContactFace* Sface,
							ContactFace* Mface,
							int Nedges, 
                                                        int* FaceEdge,
                                                        int* EdgeMaster,
                                                        Real* Sarea, 
                                                        Real* Marea )
: ContactInteractionEntity(DataArray, CT_FFI)
{
  PRECONDITION( Sface && Mface );
  num_edges   = Nedges;
  slave_face  = Sface;
  master_face = Mface;
  Set_SlaveFaceEntityData();
  Set_MasterFaceEntityData();
  if (num_edges>0) {
    int i;
    vertices = new ContactFaceFaceVertex[num_edges+1];
    for (i=0; i<num_edges; ++i) {
      vertices[i].slave_x          = Sarea[2*i];
      vertices[i].slave_y          = Sarea[2*i+1];
      vertices[i].master_x         = Marea[2*i];
      vertices[i].master_y         = Marea[2*i+1];
      vertices[i].slave_edge_id    = FaceEdge[i];
      vertices[i].master_edge_flag = EdgeMaster[i];
    }
    i = num_edges;
    vertices[i].slave_x          = Sarea[2*i];
    vertices[i].slave_y          = Sarea[2*i+1];
    vertices[i].master_x         = Marea[2*i];
    vertices[i].master_y         = Marea[2*i+1];
    vertices[i].slave_edge_id    = 0;
    vertices[i].master_edge_flag = 0;
  }
}

ContactFaceFaceInteraction::ContactFaceFaceInteraction( 
                                       ContactFaceFaceInteraction& ffi )
: ContactInteractionEntity(DataArray, CT_FFI)
{
  slave_face              = ffi.slave_face;
  master_face             = ffi.master_face;
  slave_face_entity_data  = ffi.slave_face_entity_data;
  master_face_entity_data = ffi.master_face_entity_data;
  num_edges               = ffi.num_edges;
  if (num_edges>0) {
    vertices = new ContactFaceFaceVertex[num_edges+1];
    for (int i=0; i<num_edges+1; ++i) {
      vertices[i] = ffi.vertices[i];
    }
  }
}

ContactFaceFaceInteraction* 
ContactFaceFaceInteraction::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactFace* Sface,
				     ContactFace* Mface,
				     int Nedges, 
                                     int* FaceEdge,
                                     int* EdgeMaster,
                                     Real* Sarea, Real* Marea )
{
  return new (alloc.New_Frag())
    ContactFaceFaceInteraction( Sface, Mface, Nedges, 
                                FaceEdge, EdgeMaster,Sarea, Marea );
}

ContactFaceFaceInteraction* 
ContactFaceFaceInteraction::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc )
{
  return new (alloc.New_Frag())
    ContactFaceFaceInteraction( );
}


ContactFaceFaceInteraction* 
ContactFaceFaceInteraction::new_ContactFaceFaceInteraction(
				     ContactFixedSizeAllocator& alloc,
				     ContactFaceFaceInteraction& cffi )
{
  return new (alloc.New_Frag())
    ContactFaceFaceInteraction( cffi );
}


void ContactFaceFaceInteraction_SizeAllocator(ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactFaceFaceInteraction),
                100,  // block size
                0);  // initial block size
  alloc.Set_Name( "ContactQFaceFaceInteraction allocator" );
}


ContactFaceFaceInteraction::~ContactFaceFaceInteraction()
{
  if (vertices) delete [] vertices;
}

int ContactFaceFaceInteraction::Size()
{
  return(ContactInteractionEntity::Size() + 
         2*sizeof(entity_data)+
         1*sizeof(int) + 
         DataArray_Length()*sizeof(Real) +
         (num_edges+1)*sizeof(ContactFaceFaceVertex));
}

void ContactFaceFaceInteraction::Pack( char* buffer )
{
  int cnt=0;
  ContactInteractionEntity::Pack( buffer );
  int* i_buf = reinterpret_cast<int*>(buffer + ContactInteractionEntity::Size());
  cnt += PackEntityData(&slave_face_entity_data, &i_buf[cnt]);
  cnt += PackEntityData(&master_face_entity_data, &i_buf[cnt]);
  i_buf[cnt++] = num_edges;
  
  char* buf = buffer+ContactInteractionEntity::Size()+cnt*sizeof(int);
  std::memcpy( buf, DataArray, DataArray_Length()*sizeof(Real));

  buf += DataArray_Length()*sizeof(Real);
  std::memcpy( buf, vertices,(num_edges+1)*sizeof(ContactFaceFaceVertex));
}

void ContactFaceFaceInteraction::Unpack( char* buffer )
{
  int cnt=0;
  ContactInteractionEntity::Unpack( buffer );
  int* i_buf = reinterpret_cast<int*>(buffer + ContactInteractionEntity::Size());
  cnt += UnPackEntityData(&slave_face_entity_data, &i_buf[cnt]);
  cnt += UnPackEntityData(&master_face_entity_data, &i_buf[cnt]);
  num_edges = i_buf[cnt++];
  
  char* buf = buffer+ContactInteractionEntity::Size()+cnt*sizeof(int);
  std::memcpy( DataArray, buf, DataArray_Length()*sizeof(Real));

  buf += DataArray_Length()*sizeof(Real);
  vertices  = new ContactFaceFaceVertex[num_edges+1];
  std::memcpy( vertices,buf,(num_edges+1)*sizeof(ContactFaceFaceVertex));
}

void ContactFaceFaceInteraction::Copy( ContactFaceFaceInteraction* src )
{
  ContactInteractionEntity::Copy( src );
  slave_face_entity_data  = src->slave_face_entity_data;
  master_face_entity_data = src->master_face_entity_data;
  num_edges               = src->num_edges;
  std::memcpy( DataArray,   src->DataArray, DataArray_Length()*sizeof(Real));
  vertices                = new ContactFaceFaceVertex[num_edges+1];
  std::memcpy( vertices,    src->vertices,(num_edges+1)*sizeof(ContactFaceFaceVertex));
}

void ContactFaceFaceInteraction::Connect_SlaveFace( ContactTopologyEntityList& hash_table )
{   
  slave_face = static_cast<ContactFace *>(hash_table.Find( &slave_face_entity_data ));
  POSTCONDITION( slave_face );
}

void ContactFaceFaceInteraction::Connect_MasterFace( ContactTopologyEntityList& hash_table )
{
  master_face = static_cast<ContactFace *>(hash_table.Find( &master_face_entity_data ));
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( master_face );
}

void ContactFaceFaceInteraction::Connect_SlaveFace( ContactTopologyEntityHash& hash_table )
{   
  slave_face = static_cast<ContactFace *>(hash_table.find( &slave_face_entity_data ));
  POSTCONDITION( slave_face );
}

void ContactFaceFaceInteraction::Connect_MasterFace( ContactTopologyEntityHash& hash_table )
{
  master_face = static_cast<ContactFace *>(hash_table.find( &master_face_entity_data ));
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( master_face );
}

void ContactFaceFaceInteraction::Connect_SlaveFace( ContactTopology* topology )
{   
  int block = slave_face_entity_data.block_id;
  slave_face = static_cast<ContactFace *>
    (topology->Face_Block(block)->FaceList()->Find( &slave_face_entity_data ));
  POSTCONDITION( slave_face );
}
void ContactFaceFaceInteraction::Connect_MasterFace( ContactTopology* topology )
{
  int block = master_face_entity_data.block_id;
  master_face = static_cast<ContactFace *>
    (topology->Face_Block(block)->FaceList()->Find( &master_face_entity_data ));
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( slave_face );
}

void ContactFaceFaceInteraction::Connect_SlaveFace( ContactFace* Face )
{   
  slave_face = Face;
  POSTCONDITION( slave_face );
}
void ContactFaceFaceInteraction::Connect_MasterFace( ContactFace* Face )
{
  master_face = Face;
  //Can't have this condition always be satisfied in parallel
  //POSTCONDITION( master_face );
}

int ContactFaceFaceInteraction::Data_Size()
{
  return 2+num_edges+num_edges+4*num_edges+num_edges;
}

int ContactFaceFaceInteraction::Restart_Size()
{
  return 2*sizeof(entity_data)/sizeof(int)+1+6*(num_edges+1);
}

void ContactFaceFaceInteraction::Restart_Pack( Real* buffer )
{
  int i;
  int cnt=0;
  Real* buf_loc = buffer;

  cnt += PackEntityData(&slave_face_entity_data, &buf_loc[cnt]);
  cnt += PackEntityData(&master_face_entity_data, &buf_loc[cnt]);

  // Pack the data
  buf_loc[cnt++] = num_edges;
  for (i=0; i<num_edges+1; ++i) {
    buf_loc[cnt++] = vertices[i].slave_x;
    buf_loc[cnt++] = vertices[i].slave_y;
    buf_loc[cnt++] = vertices[i].master_x;
    buf_loc[cnt++] = vertices[i].master_y;
    buf_loc[cnt++] = vertices[i].slave_edge_id;
    buf_loc[cnt++] = vertices[i].master_edge_flag;
  }
}

void ContactFaceFaceInteraction::Restart_Unpack( Real* buffer )
{
  int i;
  int cnt=0;
  Real* buf_loc = buffer;

  cnt += UnPackEntityData(&slave_face_entity_data, &buf_loc[cnt]);
  cnt += UnPackEntityData(&master_face_entity_data, &buf_loc[cnt]);
  
  // Unpack the Data Array
  num_edges = (int) *buf_loc++;
  vertices = new ContactFaceFaceVertex[num_edges+1];
  for (i=0; i<num_edges+1; ++i) {
    vertices[i].slave_x          = buf_loc[cnt++];
    vertices[i].slave_y          = buf_loc[cnt++];
    vertices[i].master_x         = buf_loc[cnt++];
    vertices[i].master_y         = buf_loc[cnt++];
    vertices[i].slave_edge_id    = (int) buf_loc[cnt++];
    vertices[i].master_edge_flag = (int) buf_loc[cnt++];
  }
}

int ContactFaceFaceInteraction::Set_SlaveFaceEntityData() 
{
  if (slave_face) {
    SetEntityData(&slave_face_entity_data, slave_face);
    return 1;
  }else{      
    return 0;
  }
}

int ContactFaceFaceInteraction::Set_MasterFaceEntityData() 
{
  if (master_face) {
    SetEntityData(&master_face_entity_data, master_face);
    return 1;
  }else{      
    return 0;
  }
}

