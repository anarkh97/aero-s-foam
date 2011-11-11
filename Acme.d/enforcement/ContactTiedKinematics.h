// $Id$

#ifndef ContactTiedKinematics_h_
#define ContactTiedKinetmaics_h_

#include "ContactSearch.h"
#include "ContactEnforcement.h"

template<typename DataType> class ContactNode;
class ContactTopology;

class ContactTiedKinematics : public ContactEnforcement {

 public:

  enum Enforcement_Data_Index { KINEMATIC_PARTITION=0, NSIZED };

  ContactTiedKinematics( const double*, ContactSearch*,
			 ContactSearch::ContactErrorCode& );

  ContactTiedKinematics( ContactSearch*, const double* restart_data,
			 ContactSearch::ContactErrorCode& );

  ~ContactTiedKinematics();

  ContactSearch::ContactErrorCode Compute_Position( double* position );

  // regression test output functions
  virtual Real Get_Global_Plot_Variable( int ) { return -1.0; };
  virtual void Get_Nodal_Plot_Variable( int, Real* );
  virtual void Get_Element_Plot_Variable( int, Real* ) { return; };

  // Restart Functions
  virtual int Number_General_Restart_Variables() {return 0;};
  virtual int Number_Nodal_Restart_Variables() {return 0;};
  virtual int Number_Edge_Restart_Variables() {return 0;};
  virtual int Number_Face_Restart_Variables() {return 0;};
  virtual int Number_Element_Restart_Variables() {return 0;};
  virtual ContactSearch::ContactErrorCode
    Extract_General_Restart_Variable( double* data );
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable( int n, double* data, int* node_ids );
  virtual ContactSearch::ContactErrorCode
    Extract_Edge_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Extract_Face_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Extract_Element_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Implant_General_Restart_Variable( double* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Edge_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Face_Restart_Variable( int n, double* data );
  virtual ContactSearch::ContactErrorCode
    Implant_Element_Restart_Variable( int n, double* data );

  virtual ContactSearch::ContactErrorCode
    Update_For_Topology_Change( int new_number_of_nodes, int* old_to_new_map );

 private:
  
  ContactTiedKinematics(ContactTiedKinematics&);
  ContactTiedKinematics& operator=(ContactTiedKinematics&);
  
  Real* final_position;
};

#endif
