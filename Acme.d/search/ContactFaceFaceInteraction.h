// $Id$

#ifndef ContactFaceFaceInteraction_h_
#define ContactFaceFaceInteraction_h_

#include "Contact_Defines.h"
#include "ContactInteractionEntity.h"
#include "ContactFace.h"
#include "ContactSearch.h"
 
struct ContactFaceFaceVertex {
  Real slave_x;
  Real slave_y;
  Real master_x;
  Real master_y;
  int  slave_edge_id;
  int  master_edge_flag;
#if (MAX_FFI_DERIVATIVES > 0)
  Real slave_x_derivatives[MAX_FFI_DERIVATIVES];
  Real slave_y_derivatives[MAX_FFI_DERIVATIVES];
  Real master_x_derivatives[MAX_FFI_DERIVATIVES];
  Real master_y_derivatives[MAX_FFI_DERIVATIVES];
#ifdef COMPUTE_FFI_SECOND_DERIVATIVES
  Real slave_x_second_derivatives[MAX_FFI_SECOND_DERIVATIVES];
  Real slave_y_second_derivatives[MAX_FFI_SECOND_DERIVATIVES];
  Real master_x_second_derivatives[MAX_FFI_SECOND_DERIVATIVES];
  Real master_y_second_derivatives[MAX_FFI_SECOND_DERIVATIVES];
#endif
#endif
};

class CString;
class ContactTopologyEntityList;
class ContactTopologyEntityHash;
class ContactHostGlobalID;
class ContactFixedSizeAllocator;

class ContactFaceFaceInteraction : public ContactInteractionEntity {
  
 public:
  
  enum InteractionSource { UNKNOWN_SOURCE=-1,CLOSEST_POINT_PROJECTION_1=1, 
                           CLOSEST_POINT_PROJECTION_2, MOVING_INTERSECTION };

  enum Scalar_Vars { UNKNOWN_SCALAR_VAR = -1,
#include "contact_variables.define"
#include "contact_variables.def"
#include "contact_variables.undefine"
                     NUMBER_SCALAR_VARS };

  enum Vector_Vars { UNKNOWN_VECTOR_VAR = -1,
#include "contact_variables.define"
#include "contact_variables.def"
#include "contact_variables.undefine"
		     NUMBER_VECTOR_VARS };

  ContactFaceFaceInteraction();
  ContactFaceFaceInteraction( ContactFace<Real>*, ContactFace<Real>*, 
			      int, int*, int*, Real*, Real* );
  ContactFaceFaceInteraction( ContactFaceFaceInteraction& );
  static ContactFaceFaceInteraction* new_ContactFaceFaceInteraction(
            ContactFixedSizeAllocator&, ContactFace<Real>*, ContactFace<Real>*, 
	    int, int*, int*, Real*, Real* );
#if (MAX_FFI_DERIVATIVES > 0)
  ContactFaceFaceInteraction( ContactFace<ActiveScalar>*, ContactFace<ActiveScalar>*,
                              int, int*, int*, ActiveScalar*, ActiveScalar* );
  static ContactFaceFaceInteraction* new_ContactFaceFaceInteraction(
            ContactFixedSizeAllocator&, ContactFace<ActiveScalar>*, ContactFace<ActiveScalar>*,
            int, int*, int*, ActiveScalar*, ActiveScalar* );
//  {
//    std::cerr << "ContactFaceFaceInteraction::new_ContactFaceFaceInteraction is not implemented for type ActiveScalar\n";
//    return NULL;
//  };
#endif
  static ContactFaceFaceInteraction* new_ContactFaceFaceInteraction(
	     ContactFixedSizeAllocator& );
  static ContactFaceFaceInteraction* new_ContactFaceFaceInteraction(
	     ContactFixedSizeAllocator&, ContactFaceFaceInteraction& );
  ~ContactFaceFaceInteraction();
  
#ifndef CONTACT_NO_MPI
  inline ContactZoltanLID& Zoltan_LID() { return zoltan_lid; };
  inline ContactZoltanGID& Zoltan_GID() { return zoltan_gid; };
  inline void ZoltanFaceLID(LB_ID_PTR lid, int flag=0) 
    {//if (flag) zoltan_lid.ZoltanLID(CT_FACE, 
     //                               master_face_entity_data.index_in_owner_proc_array, 
     //                               lid);
     //else      zoltan_lid.ZoltanLID(CT_FACE, 
     //                               master_face_entity_data.index_in_proc_array, 
     //                               lid);
     zoltan_lid.ZoltanLID(CT_FACE, 
                                    master_face_entity_data.index_in_owner_proc_array, 
                                    lid);};
  inline void ZoltanFaceGID(LB_ID_PTR gid) 
    {zoltan_gid.ZoltanGID(CT_FACE, 
                          master_face_entity_data.host_gid[0],
                          master_face_entity_data.host_gid[1],
                          gid);};
#endif

  inline ContactFace<Real>* SlaveFace() {return slave_face;};
  inline entity_data* SlaveFaceEntityData()  {return &slave_face_entity_data;};
  inline ContactFace<Real>* MasterFace() {return master_face;};
  inline entity_data* MasterFaceEntityData() {return &master_face_entity_data;};
  inline void Set_SlaveFace(ContactFace<Real>* cf) {slave_face=cf;};
  inline void Set_MasterFace(ContactFace<Real>* cf) {master_face=cf;};
  int Set_SlaveFaceEntityData();
  int Set_MasterFaceEntityData();
  inline int NumEdges() {return num_edges;};
  inline void NumEdges(int n) {num_edges=n;vertices=new ContactFaceFaceVertex[n+1];};
  inline ContactFaceFaceVertex* Get_Vertices() { return vertices; };
  inline ContactFaceFaceVertex* Get_Vertex(int n) { return &vertices[n]; };

  inline Real& Scalar_Var( VariableHandle vh ) {return DataArray[vh];};
  inline Real* Vector_Var( VariableHandle vh ) 
    { return (DataArray+NUMBER_SCALAR_VARS+3*vh); };

  inline int DataArray_Length() {return NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS;};
  inline void Initialize_Memory() {std::memset(DataArray_Buffer(), 0, DataArray_Length()*sizeof(Real));};
  void Connect_SlaveFace ( ContactTopologyEntityList& );
  void Connect_MasterFace( ContactTopologyEntityList& );
  void Connect_SlaveFace ( ContactTopologyEntityHash& );
  void Connect_MasterFace( ContactTopologyEntityHash& );
  void Connect_SlaveFace ( ContactTopology* );
  void Connect_MasterFace( ContactTopology* );
  void Connect_SlaveFace ( ContactFace<Real>* );
  void Connect_MasterFace( ContactFace<Real>* );

  // Parallel packing/unpacking functions
  int  Size();
  void Pack( char* buffer );
  void Unpack( char* buffer );
  void Copy( ContactFaceFaceInteraction* src );

  // Restart Pack/Unpack functions
  int  Restart_Size();
  void Restart_Pack( Real* buffer );
  void Restart_Unpack( Real* buffer );

  int Data_Size();
  
 protected:

 private:
  
  Real DataArray[NUMBER_SCALAR_VARS+3*NUMBER_VECTOR_VARS+1];
  ContactFace<Real>* slave_face;
  entity_data slave_face_entity_data;
  ContactFace<Real>* master_face;
  entity_data master_face_entity_data;
  int num_edges;
  ContactFaceFaceVertex* vertices;
#ifndef CONTACT_NO_MPI
  ContactZoltanLID zoltan_lid;
  ContactZoltanGID zoltan_gid;
#endif
};

#endif // ContactFaceFaceInteraction_h_
