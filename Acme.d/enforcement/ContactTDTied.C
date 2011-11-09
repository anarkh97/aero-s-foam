// $Id$


#include "ContactTDTied.h"
#include "contact_assert.h"
#if (CONTACT_DEBUG_PRINT_LEVEL>=2)
#include <iostream>
#include <cstdio>
#include "ContactParOStream.h"
#include "Contact_Communication.h"
#include "ContactTopology.h"
#endif

ContactTDTied::ContactTDTied( int ID, int* /*int_data*/, Real* real_data,
			      ContactTopology* topo)
  : ContactTDEnfModel( ID, ContactEnforcement::TD_TIED, topo )
{
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  MPI_Comm comm = topology->Search()->Get_Comm();
  if (contact_processor_number(comm)==0) {
    std::cout<<"Adding tied enforcement model, ID ="<<ID<<std::endl<<std::flush;
  }
#endif
}


ContactTDTied::ContactTDTied(ContactTopology* topo)
  : ContactTDEnfModel(ContactEnforcement::TD_TIED, topo)
{
}

ContactTDTied::~ContactTDTied()
{
}

int ContactTDTied::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  restart_data[words_added++] = id;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTDTied::Implant_Restart_Data( const Real* restart_data )
{
  int words_read = 0;
  id = (int) restart_data[words_read++];
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}

ContactSearch::ContactErrorCode
ContactTDTied::Extract_Nodal_Restart_Variable(int n, Real* data )
{ // we should never call this -- no restart variables in this model
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

ContactSearch::ContactErrorCode
ContactTDTied::Implant_Nodal_Restart_Variable(int n, const Real* data )
{ // we should never call this -- no restart variables in this model
  PRECONDITION( false);
  return ContactSearch::INVALID_ID; 
}

bool ContactTDTied::Limit_Force(ContactNodeEntityInteraction*,
Real gap, Real* rel_disp, Real* slip, Real* normal, Real dt, Real area, Real* force) 
{
    return false;
}
