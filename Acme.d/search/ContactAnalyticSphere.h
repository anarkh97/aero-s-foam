// $Id$

#ifndef ContactAnalyticSphere_h_
#define ContactAnalyticSphere_h_

#include "ContactAnalyticSurface.h"
#include "ContactSearch.h"
#include "Contact_Defines.h"

class ContactParOStream;

class ContactAnalyticSphere : public ContactAnalyticSurface {

 public:
  ContactAnalyticSphere( int id, int key, const Real* data );
  ContactAnalyticSphere( const ContactAnalyticSphere& );
  ContactAnalyticSphere( int id );
  ~ContactAnalyticSphere();

  ContactSearch::ContactErrorCode Check_for_Errors(ContactErrors*);
  ContactSearch::ContactErrorCode Set_Configuration( const Real* );

  void Bounding_Box( Real* min, Real* max );
  void ComputeBoundingBox(ContactBoundingBox*);
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

  Real center[3];
  Real radius;

};

#endif // ContactAnalyticSphere_h_
