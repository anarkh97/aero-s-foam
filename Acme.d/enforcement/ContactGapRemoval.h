// $Id$

#ifndef ContactGapRemoval_h_
#define ContactGapRemoval_h_

#include "ContactSearch.h"
#include "ContactEnforcement.h"
#include "ContactScratchManager.h"

template<typename DataType> class ContactNode;
class ContactTopology;
class ContactNodeFaceInteraction;

class ContactGapRemoval : public ContactEnforcement {
  
 public:

  enum Enforcement_Data_Index { KINEMATIC_PARTITION=0, NSIZED };
  enum Contact_Status { INACTIVE=0, ACTIVE };

  ContactGapRemoval( const double*, ContactSearch*,
		     ContactSearch::ContactErrorCode& error );
  ContactGapRemoval( ContactSearch*, const double* restart_data,
		     ContactSearch::ContactErrorCode& error );
  ~ContactGapRemoval();

  ContactSearch::ContactErrorCode Compute_Gap_Removal( int max_iterations,
						       double trivial_gap,
						       double* displc );

  // regression test output functions
  virtual Real Get_Global_Plot_Variable( int ) { return -1.0; };
  virtual void Get_Nodal_Plot_Variable( int, Real* );
  virtual void Get_Element_Plot_Variable( int, Real* ) { return; };

  // Restart Functions
  virtual int Number_General_Restart_Variables() 
    {return Number_Base_General_Restart_Variables();};
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
  
  ContactGapRemoval(ContactGapRemoval&);
  ContactGapRemoval& operator=(ContactGapRemoval&);

  void Display_Enforcement_Data();
  void Gapremoval_Set_Up();
  void Compute_Displacement_Correction();
  void Compute_Displacement_Correction_Correction();
  void Assemble_Corrections( );
  void Update_Gaps( double, bool& );
  void Make_Independent(int, ContactNodeFaceInteraction**, int* );
  double Kinematic_Partition( ContactNodeFaceInteraction* );
  int  Number_of_Contact_Constraints( ContactNodeFaceInteraction** );
  double Enforcement_Data( Enforcement_Data_Index, int, int );
  void Unify_gaps(ContactNodeFaceInteraction**, double*, double*);
  void Partition_gap(ContactNodeFaceInteraction**, double*, double*);
  double Displacement_Correction( ContactNodeFaceInteraction* );
  void Invert_3x3Matrix(double m[3][3], double i[3][3]);

  // Internal memory
  double* displ_correction;
  int max_interactions;
  int num_w_sing_gap;
  int num_w_mult_gap;
  ContactNodeFaceInteraction** single_contact_constraints;
  ContactNodeFaceInteraction** multiple_contact_constraints;

  // Scratch_Variables
  ScratchVariable D_COR;
  ScratchVariable D_CORCOR;
  ScratchVariable TOTAL_COR;
  ScratchVariable UPDATED_POS;  
  ScratchVariable GAP_TO_ENFORCE;
};

#endif
