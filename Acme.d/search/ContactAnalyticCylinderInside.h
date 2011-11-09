// $Id$

#ifndef ContactAnalyticCylinderInside_h_
#define ContactAnalyticCylinderInside_h_

#include "ContactAnalyticSurface.h"
#include "ContactSearch.h"
#include "Contact_Defines.h"

class ContactParOStream;

class ContactAnalyticCylinderInside : public ContactAnalyticSurface {

 public:

  ContactAnalyticCylinderInside( int id, int key, const Real* data );
  ContactAnalyticCylinderInside( const ContactAnalyticCylinderInside& );
  ContactAnalyticCylinderInside( int ID );
  ~ContactAnalyticCylinderInside();

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
  
  Real center[3];
  Real axial_dir[3];
  Real radius;
  Real length;

};

#endif // ContactAnalyticCylinderInside_h_
