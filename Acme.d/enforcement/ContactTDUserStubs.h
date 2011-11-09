// $Id$


#include "Contact_Defines.h"
#include "ContactTDUserSubTypes.h"

extern "C" {
  
  CONTACT_INIT_MODEL_FN                 FORTRAN(user_initialize_model);
  CONTACT_INIT_TIME_STEP_FN             FORTRAN(user_initialize_time_step);
  CONTACT_INIT_NODE_STATE_DATA_FN       FORTRAN(user_init_node_state_data);
  CONTACT_LIMIT_FORCE_FN                FORTRAN(user_limit_force);
  CONTACT_INTERACTION_ACTIVE_FN         FORTRAN(user_is_active);
  CONTACT_INTERACTION_TYPE_FN           FORTRAN(user_interaction_type);
  
}

