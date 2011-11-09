// $Id$

#include "ContactTDEnfModel.h"
#include "ContactEnforcement.h"
#include "ContactEnfModel.h"

ContactTDEnfModel::ContactTDEnfModel( int ID,
         ContactEnforcement::Enforcement_Model_Types Type,
				      ContactTopology* Topology )
  : ContactEnfModel( ID, Type, Topology )
{
}

ContactTDEnfModel::ContactTDEnfModel( 
	 ContactEnforcement::Enforcement_Model_Types Type,
	 ContactTopology* Topology )
  : ContactEnfModel( Type, Topology )
{
}

ContactTDEnfModel::~ContactTDEnfModel()
{}

bool 
ContactTDEnfModel::Active_Interaction
(ContactNodeEntityInteraction* cnfi,Real gap)
{
  return (gap < 0 ? true : false);
}

void 
ContactTDEnfModel::Set_Node_State_Data(Real value, int node, 
                                       int offset, int state)
{
  Node_State_Data(node,state)[offset] = value;
}

Real 
ContactTDEnfModel::Get_Node_State_Data(int node, int offset, int state)
{
  return Node_State_Data(node,state)[offset];
}

#if 0
ContactParOStream& 
ContactTDEnfModel::ParOStream()
{
	  return search->postream;
}
#endif

