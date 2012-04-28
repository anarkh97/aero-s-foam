// $Id$

#include "ContactFaceFaceSearch.C"

#define INSTANTIATION_HELPER(T) \
template ContactFaceFaceInteraction* ContactSearch::Face_Face_Search<T>(ContactFace<T>* slave_face, \
                                                 ContactFace<T>* master_face, \
                                                 ContactElem<T>* element, \
                                                 VariableHandle POSITION); \
template bool ContactSearch::Face_Face_Search_Step1<T>(ContactFace<T>* slave_face, \
                                                 ContactFace<T>* master_face, \
                                                 ContactElem<T>* element, \
                                                 VariableHandle POSITION, int&, int&); \
template ContactFaceFaceInteraction* ContactSearch::Face_Face_Search_Step2<T>(ContactFace<T>* slave_face, \
                                                 ContactFace<T>* master_face, \
                                                 ContactElem<T>* element, \
                                                 VariableHandle POSITION, int, int);

INSTANTIATION_HELPER(Real);
#if (MAX_FFI_DERIVATIVES > 0)
INSTANTIATION_HELPER(ActiveScalar);
#endif
