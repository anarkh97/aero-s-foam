// $Id$

#ifndef contact_assert_h_
#define contact_assert_h_

#ifdef CONTACT_DEBUG

#include <cstdio>
#include <cstdlib>

#ifndef PRECONDITION
#define PRECONDITION(e) \
  (!(e)? std::fprintf(stderr,"Precondition not satisfied: file \"%s\",\
  line %d\n%s\n", __FILE__, __LINE__, #e), (void)std::abort() : \
  (void)0)
#endif

#ifndef POSTCONDITION
#define POSTCONDITION(e) \
  (!(e)? std::fprintf(stderr,"Postcondition not satisfied: file \"%s\",\
  line %d\n%s\n", __FILE__, __LINE__, #e), (void)std::abort() : \
  (void)0)
#endif

#ifndef INVARIANT
#define INVARIANT(e) \
  (!(e)? std::fprintf(stderr,"Postcondition not satisfied: file \"%s\",\
  line %d\n%s\n", __FILE__, __LINE__, #e), (void)std::abort() : \
  (void)0)
#endif

#ifndef REMEMBER
#define REMEMBER(e) e
#endif

#else

#ifndef PRECONDITION
#define PRECONDITION(e) 
#endif

#ifndef POSTCONDITION
#define POSTCONDITION(e) 
#endif

#ifndef INVARIANT
#define INVARIANT(e) 
#endif

#ifndef REMEMBER
#define REMEMBER(e)
#endif

#endif // #ifdef CONTACT_DEBUG

#endif // #ifndef contact_asserts_h_


