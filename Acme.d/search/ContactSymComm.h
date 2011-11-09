// $Id$

#ifndef ContactSymComm_h_
#define ContactSymComm_h_

#include "Contact_Defines.h"
#include "contact_assert.h"

#ifndef CONTACT_NO_MPI
#include "mpi.h"
class ContactTopologyEntity;
class ContactTopologyEntityHash;
class ContactTopologyEntityList;
class ContactHostGlobalID;

class ContactSymComm {
  
 public:
  ContactSymComm( int, const int*, const int*, ContactTopologyEntity** );
  ContactSymComm();
  ~ContactSymComm();
  
  inline int Size() {return num_entities;};
  inline int Num_Comm_Partners() {return num_comm_partners;};
  inline int Comm_Proc_ID( int i ){
    PRECONDITION( i<num_comm_partners || num_comm_partners==0 );
    return comm_proc_ids[i];
  };

  inline int Num_to_Proc( int i ){
    PRECONDITION( i<num_comm_partners || num_comm_partners==0 );
    return num_to_proc[i];
  }
  inline ContactTopologyEntity** Entity_List( int i ){
    PRECONDITION( i<num_comm_partners || num_comm_partners==0 );
    if (num_comm_partners==0)
      return NULL;
    else
      return entity_list+offset[i];
  }

  void Build_Comm_Plan( ContactSymComm*, int, int, int*, ContactHostGlobalID*,
			ContactTopologyEntityList*, MPI_Comm& );
  void Build_Comm_Plan( ContactSymComm*, int, int, int*, ContactHostGlobalID*,
			ContactTopologyEntityHash*, MPI_Comm& );
  void Build_Subset_Comm_Plan_Using_Temp_Tag( ContactSymComm& );

  void Check(MPI_Comm Comm);
  void Sort();

 private:
  
  ContactSymComm(ContactSymComm&);
  ContactSymComm& operator=(ContactSymComm&);

  int num_comm_partners;    // How many procs I communicate to
  int num_entities;         // Number of entities in entity_list
  int* comm_proc_ids;       // Processor ids (num_com_partners long)
  int* num_to_proc;         // Number to each proc (num_com_partners long)
  int* offset;              // Offset in entity_list for comm_partner i 
  ContactTopologyEntity** entity_list;

  int allocated_procs;
  int allocated_entities;

};

#endif
#endif
