// $Id$

#ifndef ContactAnalyticPlane_h_
#define ContactAnalyticPlane_h_

#include "ContactAnalyticSurface.h"
#include "ContactSearch.h"
#include "Contact_Defines.h"

class ContactParOStream;

class ContactAnalyticPlane : public ContactAnalyticSurface {

 public:

  ContactAnalyticPlane( int id, int key, const Real* data );
  ContactAnalyticPlane( const ContactAnalyticPlane& );
  ContactAnalyticPlane( int id );
  ~ContactAnalyticPlane();

  ContactSearch::ContactErrorCode Check_for_Errors( ContactErrors* );
  ContactSearch::ContactErrorCode Set_Configuration( const Real* data );

  void Bounding_Box( Real* min, Real* max );
  void ComputeBoundingBox(ContactBoundingBox* bb);
  virtual bool Process( Real* node_position,
                        Real& penetration_mag,
                        Real* contact_point,
                        Real* surface_normal,
                        Real* pushback_dir,
                        Real& time_to_contact,
                        int&  location,
                        bool PRINT_THIS_NODE,
                        ContactParOStream* postream );
  virtual bool Process( Real* node_position,    //1st configuration
                        Real* node_position_2,  //2nd configuration (aug or predicted); 
                        Real& penetration_mag, 
			Real* contact_point,
                        Real* surface_normal,
                        Real* pushback_dir,
                        Real& time_to_contact,
                        int&  location,
                        bool PRINT_THIS_NODE,
                        ContactParOStream* postream );
  void Display( ContactParOStream& );

  virtual int Restart_Size();
  virtual int Extract_Restart_Data( Real* );
  virtual int Implant_Restart_Data( Real* );

 private:
  
  Real point[3];
  Real normal_vector[3];

};

#endif // ContactAnalyticPlane_h_
