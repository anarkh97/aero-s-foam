// $Id$

#include "ContactTable.h"
#include "contact_assert.h"
#include <cstring>

ContactTable::ContactTable( int ID, int Num_Points, Real* Abscissa, 
			    Real* Ordinate )
{
  PRECONDITION( ID > 0 );
  PRECONDITION( Num_Points > 0 );
  PRECONDITION( Abscissa );
  PRECONDITION( Ordinate );
  id = ID;
  num_points = Num_Points;
  abscissa = new Real[num_points];
  ordinate = new Real[num_points];
  std::memcpy( abscissa, Abscissa, num_points*sizeof(Real) );
  std::memcpy( ordinate, Ordinate, num_points*sizeof(Real) );
}


ContactTable::~ContactTable()
{
  delete [] abscissa;
  delete [] ordinate;
}

Real ContactTable::Interpolate_Value( Real Abscissa )
{
  for( int i=0 ; i<num_points-1 ; ++i){
    if( Abscissa >= abscissa[i] && Abscissa <= abscissa[i+1] ){
      // linearly interpolate
      return( ordinate[i] +
	      ((Abscissa      - abscissa[i] )/
	       (abscissa[i+1] - abscissa[i] )) *
	       (ordinate[i+1] - ordinate[i]) );
    }
  }

  // Out of the bounds.  Return the value at the bound
  if( Abscissa < abscissa[0] )
    return ordinate[0];
  else
    return ordinate[num_points-1];
}


int ContactTable::Extract_Restart_Data( Real* restart_data )
{
  int words_added = 0;
  restart_data[words_added++] = id;
  restart_data[words_added++] = num_points;
  std::memcpy( &restart_data[words_added], abscissa, num_points*sizeof(Real) );
  words_added += num_points;
  std::memcpy( &restart_data[words_added], ordinate, num_points*sizeof(Real) );
  words_added += num_points;
  POSTCONDITION( words_added == Restart_Size() );
  return words_added;
}

int ContactTable::Implant_Restart_Data( Real* restart_data )
{
  int words_read = 0;
  id = (int) restart_data[words_read++];
  num_points = (int) restart_data[words_read++];
  abscissa = new Real[num_points];
  ordinate = new Real[num_points];
  std::memcpy( abscissa, &restart_data[words_read], num_points*sizeof(Real) );
  words_read += num_points;
  std::memcpy( ordinate, &restart_data[words_read], num_points*sizeof(Real) );
  words_read += num_points;
  POSTCONDITION( words_read == Restart_Size() );
  return words_read;
}
