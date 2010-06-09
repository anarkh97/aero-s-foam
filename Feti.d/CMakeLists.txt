ADD_LIBRARY(Feti
CoarseSetCore.C
DistrVectorCore.C
FetiCore.C
FetiOpCore.C
NLFetiCore.C)

SET(INC_DIR ..)

INCLUDE_DIRECTORIES( ${INC_DIR} )
