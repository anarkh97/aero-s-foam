// $Id$

#ifndef ContactTDFrictionless_h_
#define ContactTDFrictionless_h_

#include "ContactTDEnfModel.h"
#include "ContactSearch.h"

class ContactTDFrictionless : public ContactTDEnfModel {

 public:
  ContactTDFrictionless( int ID, int* int_data, Real* real_data,
			 ContactTopology* topology );
  ContactTDFrictionless( ContactTopology* topology );
  ~ContactTDFrictionless();

  virtual int Num_Interaction_State_Variables() { return 0; };
  virtual int Num_Node_State_Variables() { return 0; };
  virtual ContactSearch::ContactErrorCode 
    Initialize_Model( int num_models, ContactEnfModel** models,
		      int num_tables, ContactTable** tables )
    { return ContactSearch::NO_ERROR; };
  // virtual void Initialize_Interaction_State_Data(Real*) {};
  virtual void Initialize_Node_State_Data(Real* state_data) {};
  virtual int Restart_Size() { return 1; };
  virtual int Extract_Restart_Data(double* restart_data);
  virtual int Implant_Restart_Data(const double* restart_data);
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable(int n, Real* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable(int n, const Real* data );
  virtual int  Interaction_Type(ContactNodeEntityInteraction*)
    {return TDEM_FRICTIONLESS;}
  virtual bool Limit_Force(ContactNodeEntityInteraction*, Real, Real*, Real*,Real*,Real, Real, Real*);
};
#endif
