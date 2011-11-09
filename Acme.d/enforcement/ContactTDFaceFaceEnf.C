// $Id$

#ifdef CONTACT_TD_FACE_FACE_ENF

#include "ContactTDFaceFaceEnf.h"
#include "ContactErrors.h"
#include "ContactFace.h"
#include "ContactNode.h"
#include "ContactFaceFaceInteraction.h"

#include <cstring>

ContactTDFaceFaceEnf::ContactTDFaceFaceEnf( const Real* enf_data,
					    ContactSearch* Search,
			     ContactSearch::ContactErrorCode& error )
  : ContactEnforcement( Search, ContactEnforcement::TDFaceFaceEnf,
			NSIZED, enf_data, false )
{
}

ContactTDFaceFaceEnf::ContactTDFaceFaceEnf( ContactSearch* Search,
					    const Real* restart_data,
			     ContactSearch::ContactErrorCode& error)
  : ContactEnforcement( Search, ContactEnforcement::TDFaceFaceEnf,
			restart_data )
{  
}

ContactTDFaceFaceEnf::~ContactTDFaceFaceEnf()
{
}


ContactSearch::ContactErrorCode 
ContactTDFaceFaceEnf::Compute_Forces( const Real* mass, Real* force )
{
  int i;

  ContactSearch::ContactErrorCode error_code = ContactSearch::NO_ERROR;

  // Initialize the force
  std::memset( force, 0, number_of_nodes*3*sizeof(Real) );

  // Call base class Set_Up function to prepare for enforcement
  error_code = Set_Up();
  if( error_code ) return error_code;

  // Loop over the face-face interactions
  for( i=0 ; i<number_face_face_interactions ; ++i ){
    ContactFaceFaceInteraction* ffi = face_face_interaction_list[i];
    ContactFace* slave_face = ffi->SlaveFace();
    ContactFace* master_face = ffi->MasterFace();
    POSTCONDITION( slave_face && master_face );
    std::cout << "Slave  Face = " << slave_face->Global_ID() << std::endl;
    std::cout << "Master Face = " << master_face->Global_ID() << std::endl;
  }
  errors->Add_Error_Message("ContactTDFaceFaceEnf::Compute_Forces() is not implemented");
  error_code = ContactSearch::UNIMPLEMENTED_FUNCTION;
  return error_code;
}

#endif
