// $Id$

#ifndef ContactEnforcementData_h_
#define ContactEnforcementData_h_

#include "Contact_Defines.h"
#include "contact_assert.h"

class ContactTopology;

class ContactEnforcementData {

 public:
  ContactEnforcementData( int number_entity_keys, int size_data_per_pair,
			  const Real* Data );
  ContactEnforcementData();
  ~ContactEnforcementData();

  Real* Data_Array() { return data; };
  Real Get_Data( int index, int f_key, int n_key  )
    { PRECONDITION (f_key >= 0);
      PRECONDITION (n_key >= 0);
      PRECONDITION (f_key < number_entity_keys);
      PRECONDITION (n_key < number_entity_keys);
      PRECONDITION (index < size_data_per_pair );
      return data[(f_key*number_entity_keys+n_key)*size_data_per_pair+index]; 
    };
  void Set_Data( int index, int f_key, int n_key , Real value )
    { PRECONDITION (f_key >= 0);
      PRECONDITION (n_key >= 0);
      PRECONDITION (f_key < number_entity_keys);
      PRECONDITION (n_key < number_entity_keys);
      PRECONDITION (index < size_data_per_pair );
      data[(f_key*number_entity_keys+n_key)*size_data_per_pair+index]=value; 
    };
  int Number_Entity_Keys() { return number_entity_keys; };
  int Size_Data_Per_Pair() { return size_data_per_pair; };
  int Restart_Size() { return 2+number_entity_keys*number_entity_keys*
			      size_data_per_pair; };
  int Extract_Restart_Data( Real* restart_data );
  int Implant_Restart_Data( const Real* restart_data );

 private:
  
  ContactEnforcementData(ContactEnforcementData&);
  ContactEnforcementData& operator=(ContactEnforcementData&);
  
  int number_entity_keys;
  int size_data_per_pair;
  Real* data;
  
};

#endif
