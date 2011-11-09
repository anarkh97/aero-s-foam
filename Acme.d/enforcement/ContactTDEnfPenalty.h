// $Id$

#ifndef ContactTDEnfPenalty_h_
#define ContactTDEnfPenalty_h_

#include "ContactSearch.h"
#include "ContactEnforcement.h"
#include "ContactTDEnforcement.h"
#include "ContactTDEnfModel.h"

class ContactNode;
class ContactNodeFaceInteraction;
class ContactNodeSurfaceInteraction;
class ContactNodeEntityInteraction;
class ContactTopology;

class ContactTDEnfPenalty : public ContactTDEnforcement {
  
 public:
  ContactTDEnfPenalty( const double*, ContactSearch*,
 			     ContactSearch::ContactErrorCode& error,
                             bool get_cvars, 
                             bool get_plot_force );
  ContactTDEnfPenalty( ContactSearch*, const double* restart_data,
			     ContactSearch::ContactErrorCode& error,
                             bool get_cvars,
                             bool get_plot_force );
  ~ContactTDEnfPenalty();

  ContactSearch::ContactErrorCode Compute_Contact_Force( double dt_old, 
							 double dt,
							 double* mass,
							 double* density,
							 double* wavespeed,
							 double* force );
   ContactSearch::ContactErrorCode Set_Penalty_Scale( double scale );

 private:
  void TDPenalty_Release_Scratch(void);

  // Nodal Scratch Variable Handles
  ScratchVariable NODAL_MASS;
  ScratchVariable INC_FORCE;
  ScratchVariable TOTAL_FORCE;
//  ScratchVariable NEW_POSITION;

  ContactNodeFaceInteraction* cnfi;
  ContactNodeEntityInteraction* cnei;
  double penalty_scale;

};

#endif

