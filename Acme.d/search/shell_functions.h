// $Id$

// These functions override base class virtual functions for shells.
// I chose to do this with an include file instead of having a 
// diamond multiple inheritance.


 public: 

#ifndef CONTACT_NO_MPI
  virtual int Size() { return (ContactFace<DataType>::Size() + sizeof(Real)); };
  virtual void Pack( char* buffer ){ 
    ContactFace<DataType>::Pack(buffer);
    std::memcpy( buffer+ContactFace<DataType>::Size(), &thickness, sizeof(Real) );
  };
  virtual void Unpack( char* buffer ){
    ContactFace<DataType>::Unpack(buffer);
    std::memcpy( &thickness, buffer+ContactFace<DataType>::Size(), sizeof(Real) );
  };
#endif

  Real Thickness() { return thickness; };
  void Thickness(Real t) {thickness = t;};

  Real Lofting_Factor() { return lofting_factor; };
  void Lofting_Factor( Real lf ) { lofting_factor = lf; };

 private:

  Real thickness;
  Real lofting_factor;

