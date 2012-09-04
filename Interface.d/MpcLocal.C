#include <MpcLocal.H>
#include <zero.H>
#include <echo.h>
// $Date: 2006/04/13 22:49:39 $
// $Revision: 1.2 $
#include <maps.H>
#include <salinas.h>
#include <code_types.h>

/**********************************************************************/
// resizes the object. keeps current data
void MpcLocal::NumEntries( int n )
{
  Int *tmpId=0;
  Real* tmpC=0;
  short* tmpDof=0;
  int copyn=0;

  if (n>0 ) {
    if ( mNumEntries ) {
      tmpId=mLocalId;
      tmpC =mCoef;
      tmpDof=mNodalDof;
      copyn=mNumEntries;
    }
    mLocalId=NEW Int[n];
    mCoef=NEW Real[n];
    mNodalDof=NEW short[n];
    mNumEntries=n;
    zero(mCoef,n);
    zero(mLocalId,n);
    zero(mNodalDof,n);
    if ( copyn ){
      if ( copyn > n ) copyn=n;
      memcpy(mLocalId,tmpId,copyn*sizeof(Int));
      memcpy(mCoef,tmpC,copyn*sizeof(Real));
      memcpy(mNodalDof,tmpDof,copyn*sizeof(short));
      delete[] tmpId; tmpId = 0;
      delete[] tmpC; tmpC = 0;
      delete[] tmpDof; tmpDof = 0;
    }
  }
  else { // n=0
    delete[] mLocalId; mLocalId = 0;
    delete[] mCoef; mCoef = 0;
    delete[] mNodalDof; mNodalDof = 0;
    mLocalId=NULL;
    mCoef=NULL;
    mNodalDof=NULL;
    mNumEntries=n;
  }
}
/**********************************************************************/
MpcLocal::MpcLocal( const MpcLocal& old )        // copy constructor
{
  mNumEntries = old.NumEntries();
  mNumEntriesGlobal = old.NumEntriesGlobal();
  if ( mNumEntries ) 
  {
    mLocalId   = NEW Int[mNumEntries];    
    mNodalDof   = NEW short[mNumEntries]; 
    mCoef       = NEW Real[mNumEntries];  
  }
  for (unsigned int i=0; i<mNumEntries; i++) {
    mLocalId[i] = old.LocalId(i);
    mNodalDof[i] = old.NodalDof(i);
    mCoef[i]     = old.Coef(i);
  }
  G( old.G() );
}
/**********************************************************************/
MpcLocal& MpcLocal::operator=(const MpcLocal& old)
  // overloaded assignment operator
{
  if (&old == this) // Check for assignment to self in operator =
    return *this;
  // first destroy what's in this
  if (mNumEntries){
    delete[] mLocalId; mLocalId = 0;
    delete[] mCoef; mCoef = 0;
    delete[] mNodalDof; mNodalDof = 0;
  }
  // now, allocate for new object of size of rhs
  mNumEntries = old.NumEntries();
  mNumEntriesGlobal = old.NumEntriesGlobal();
  if ( mNumEntries )
  {
    mLocalId   = NEW Int[mNumEntries];     
    mNodalDof   = NEW short[mNumEntries]; 
    mCoef       = NEW Real[mNumEntries]; 
  }
  // finally, copy rhs data into this
  for (unsigned int i=0; i<mNumEntries; i++) {
    mLocalId[i] = old.LocalId(i);
    mNodalDof[i] = old.NodalDof(i);
    mCoef[i]     = old.Coef(i);
  }
  G( old.G() );

  return *this;
}
