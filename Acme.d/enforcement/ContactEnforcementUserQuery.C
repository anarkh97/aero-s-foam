
#include "Contact_Defines.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactEnforcement.h"
#include "ContactTable.h"

void FORTRAN(UserQuery_Table_Last_Abscissa)(void *data, int *id, Real *value)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)data;
  ContactSearch*      search      = enforcement->Search();
  int num_tables = search->Num_Tables();
  ContactTable** tables = search->Tables();
  for( int i=0 ; i<num_tables ; ++i ){
    ContactTable* table = tables[i];
    if( *id == table->ID() ){
      *value = table->Last_Abscissa();
      break;
    }
  }
}
