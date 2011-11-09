// $Id$

#ifndef ContactAnalyticSurface_h_
#define ContactAnalyticSurface_h_ 

#include "Contact_Defines.h"
#include "ContactTopologyEntity.h"
#include "ContactSearch.h"
#include "ContactBoundingBox.h"

class ContactParOStream;

// The enum for the AnalyticSurface types lives in ContactSearch so that
// host codes only need to include that file to be able to use the
// search package.

class ContactAnalyticSurface : public ContactTopologyEntity {

 public:

  ContactAnalyticSurface( ContactSearch::AnalyticSurface_Type, int, int );
  ContactAnalyticSurface( ContactSearch::AnalyticSurface_Type, 
                          ContactAnalyticSurface* );
  ContactAnalyticSurface( ContactSearch::AnalyticSurface_Type, int );
  virtual ~ContactAnalyticSurface();

  virtual ContactSearch::ContactErrorCode Check_for_Errors(ContactErrors*) = 0;
  virtual ContactSearch::ContactErrorCode 
    Set_Configuration( const Real* data ) = 0;

  inline ContactSearch::AnalyticSurface_Type Surface_Type() {return type;};
  
  inline int DataArray_Length() {return 0;};
  
  // This function allows the search to eliminate nodes that are not near
  // the analytic surface to be culled out
  virtual void Bounding_Box( Real* min, Real* max ) = 0;

  // This function process a node in proximity to the analytic surface and
  // determines if an interaction should be defined.
  virtual bool Process( Real* node_position,
                        Real& penetration_mag, 
			Real* contact_point,
                        Real* surface_normal,
                        Real* pushback_dir,
                        Real& time_to_contact,
                        int&  location,
                        bool PRINT_THIS_NODE = false,
                        ContactParOStream* postream = NULL ) = 0;

  // This function process a node in proximity to the analytic surface and
  // determines if an interaction should be defined.
  virtual bool Process( Real* node_position,    //1st configuration
                        Real* node_position_2,  //2nd configuration (aug or predicted); 
                        Real& penetration_mag, 
			Real* contact_point,
                        Real* surface_normal,
                        Real* pushback_dir,
                        Real& time_to_contact,
                        int&  location,
                        bool PRINT_THIS_NODE,
                        ContactParOStream* postream ) = 0;
                        
  // This function is for output information to stdout

  virtual void Display( ContactParOStream& ) = 0;

  virtual int Restart_Size() = 0;
  virtual int Extract_Restart_Data( Real* restart_data ) = 0;
  virtual int Implant_Restart_Data( Real* restart_data ) = 0;
  virtual void ComputeBoundingBox(ContactBoundingBox*) = 0;
  ContactBoundingBox* BoundingBox() { return &bounding_box; };

 protected:
  ContactBoundingBox bounding_box;
 
 private:
  ContactSearch::AnalyticSurface_Type type;

};

#endif // ContactAnalyticalSurface_h_
