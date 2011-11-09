// $Id$

#ifndef ContactTDJunction_h_
#define ContactTDJunction_h_

#include "ContactTDEnfModel.h"
#include "ContactTable.h"
#include "ContactSearch.h"

class ContactTDJunction : public ContactTDEnfModel {

 public:
  ContactTDJunction( int ID, int* int_data, Real* real_data,
			 ContactTopology* topology );
  ContactTDJunction( ContactTopology* topology );
  ~ContactTDJunction();

  virtual int Num_Interaction_State_Variables() { return 0; };
  virtual int Num_Node_State_Variables() { return 0; };
  virtual ContactSearch::ContactErrorCode 
    Initialize_Model( int num_models, ContactEnfModel** models,
		      int num_tables, ContactTable** tables );
  // virtual void Initialize_Interaction_State_Data(Real*) {};
  virtual void Initialize_Node_State_Data(Real* state_data) {};

  virtual int Restart_Size() { return 4; };
  virtual int Extract_Restart_Data(double* restart_data);
  virtual int Implant_Restart_Data(const double* restart_data);
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable(int n, Real* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable(int n, const Real* data );

  virtual int  Interaction_Type(ContactNodeEntityInteraction*)
  {return TDEM_ADHESIVE_FRICTION;};
  virtual bool Limit_Force(ContactNodeEntityInteraction*,
        Real, Real*, Real*,Real*,Real, Real, Real*);
  virtual bool Active_Interaction(ContactNodeEntityInteraction*,Real);

 private:
  int Ntable_id; // ContactTable id for F_n
  int Ttable_id; // ContactTable id for F_t
  Real NCutoff; // adhesion cutoff [length]
  Real TCutoff; // shear junction normal distance cutoff [length]
  ContactTable* F_n; // (normal) force-displ. law
  ContactTable* F_t; // (tangential) force-displ. law


};
#endif
