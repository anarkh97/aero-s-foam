#ifndef FWK_MACROS_H
#define FWK_MACROS_H

// A macro to disallow the copy constructor and operator= functions
#define DISALLOW_COPY_AND_ASSIGN(TypeName)   \
  TypeName(const TypeName &);                 \
  const TypeName & operator=(const TypeName &)

#define EXPORT_PTRINTERFACE_TYPES(Typename) \
  typedef Fwk::Ptr< Typename > Ptr;           \
  typedef Fwk::Ptr< const Typename > PtrConst  


#endif /* FWK_MACROS_H */
