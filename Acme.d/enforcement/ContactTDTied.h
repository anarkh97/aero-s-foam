// $Id$

#ifndef ContactTDTied_h_
#define ContactTDTied_h_

#include "ContactTDEnfModel.h"
#include "ContactSearch.h"
#include "contact_assert.h"

class ContactTDTied : public ContactTDEnfModel {

 public:
  ContactTDTied( int ID, int* int_data, Real* real_data, 
		 ContactTopology* Topology );
  ContactTDTied( ContactTopology* Topology );
  ~ContactTDTied();

  virtual int Num_Interaction_State_Variables() { return 0; };
  virtual int Num_Node_State_Variables() { return 0; };
  virtual ContactSearch::ContactErrorCode 
    Initialize_Model( int num_models, ContactEnfModel** models,
		      int num_tables, ContactTable** tables )
    { return ContactSearch::NO_ERROR; };
  // virtual void Initialize_Interaction_State_Data(Real*) {};
  virtual void Initialize_Node_State_Data(Real*) {};
  virtual int Restart_Size() { return 1; };
  virtual int Extract_Restart_Data(Real* restart_data);
  virtual int Implant_Restart_Data(const Real* restart_data);
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable(int n, Real* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable(int n, const Real* data );

  virtual int  Interaction_Type  (ContactNodeEntityInteraction*)       {return TDEM_TIED;}
  virtual bool Active_Interaction(ContactNodeEntityInteraction*, Real) {return true;}
  virtual bool Limit_Force(ContactNodeEntityInteraction*, Real, Real*, Real*,Real*,Real, Real, Real*);
};

#endif
